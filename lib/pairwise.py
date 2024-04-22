import numpy as np
import scipy.special
import concurrent.futures
import tqdm
import itertools

from common import Models, NUM_MODELS, ALL_MODELS
import common
import lh
import mutrel

def swap_A_B(arr):
  swapped = np.zeros(len(arr)) + np.nan
  for midx, M in enumerate(ALL_MODELS):
    if midx in (Models.garbage, Models.cocluster, Models.diff_branches):
      swapped[midx] = arr[midx]
    elif midx == Models.A_B:
      swapped[midx] = arr[Models.B_A]
    elif midx == Models.B_A:
      swapped[midx] = arr[Models.A_B]
    else:
      raise Exception('Unknown model')
  assert not np.any(np.isnan(swapped))
  return swapped

def _calc_posterior(evidence, logprior):
  # If both variants have zero variant reads, the evidence will be all zeroes.
  # Ensure we don't hit this case when working with coclustering.
  if evidence[Models.cocluster] == 0 and np.all(np.delete(evidence, Models.cocluster) == -np.inf):
    # Regardless of what the prior is on coclustering (even if it's zero), if
    # we're dealing with a pair consisting of variant `i` and variant `i`, we
    # want the posterior to indicate coclustering certainty.
    posterior = np.zeros(len(evidence))
    posterior[Models.cocluster] = 1
    return posterior

  joint = evidence + logprior
  B = np.max(joint)
  joint -= B
  posterior = np.exp(joint) / np.sum(np.exp(joint))
  return posterior

def _calc_posterior_full(evidence, logprior):
  # This function is currently used only to double-check the results of
  # `_calc_posterior`.
  joint = evidence + logprior[None,None,:]
  diag = range(len(joint))
  joint[diag,diag,:] = -np.inf
  joint[diag,diag,Models.cocluster] = 0

  B = np.max(joint, axis=2)
  joint -= B[:,:,None]
  expjoint = np.exp(joint)
  posterior = expjoint / np.sum(expjoint, axis=2)[:,:,None]

  mutrel.check_posterior_sanity(posterior)
  return posterior

def make_full_posterior(evidence, logprior):
  logprior = _complete_logprior(logprior)
  posterior = mutrel.Mutrel(
    vids = evidence.vids,
    rels = _calc_posterior_full(evidence.rels, logprior),
  )
  return posterior


def _complete_logprior(logprior):
  '''
  If you don't want to include a certain model in the posterior (e.g., perhaps
  you're computing pairwise probabilities for supervariants, and so don't want
  to allow them to cocluster unless they're the same supervariant), set the
  corresponding entry in prior to 0.
  '''
  if logprior is None:
    logprior = {}
    logtotal = -np.inf
  else:
    for K in logprior.keys():
      assert K in ALL_MODELS
      assert logprior[K] <= 0
    logtotal = scipy.special.logsumexp(list(logprior.values()))

  assert logtotal <= 0
  remaining = set(ALL_MODELS) - set(logprior.keys())
  if np.isclose(0, logtotal):
    remaining_val = -np.inf
  else:
    remaining_val = np.log(1 - np.exp(logtotal)) - np.log(len(remaining))
  for K in remaining:
      logprior[K] = remaining_val

  logprior_vals = -np.inf * np.ones(len(logprior))
  for K in logprior:
    logprior_vals[getattr(Models, K)] = logprior[K]
  assert np.isclose(0, scipy.special.logsumexp(logprior_vals))
  return logprior_vals

def _compute_pairs(pairs, variants, logprior, posterior, evidence, pbar=None, parallel=1):
  logprior = _complete_logprior(logprior)
  # TODO: change ordering of pairs based on what will provide optimal
  # integration accuracy according to Quaid's advice.
  pairs = list(pairs)
  # Don't bother starting more workers than jobs.
  parallel = min(parallel, len(pairs))

  # If you set parallel = 0, we don't invoke the parallelism machinery. This
  # makes debugging easier.
  if parallel > 0:
    futures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=parallel) as ex:
      for A, B in pairs:
        futures.append(ex.submit(_calc_lh_and_posterior, variants[A], variants[B], logprior))
      if pbar is not None:
        for F in concurrent.futures.as_completed(futures):
          pbar.update()
    for (A, B), F in zip(pairs, futures):
      evidence.rels[A,B], posterior.rels[A,B] = F.result()
  else:
    for A, B in pairs:
      evidence.rels[A,B], posterior.rels[A,B] = _calc_lh_and_posterior(variants[A], variants[B], logprior)

  # Duplicate evidence's keys, since we'll be modifying dictionary.
  for A, B in pairs:
    if A == B:
      continue
    evidence.rels[B,A] = swap_A_B(evidence.rels[A,B])
    posterior.rels[B,A] = swap_A_B(posterior.rels[A,B])

  mutrel.check_mutrel_sanity(evidence.rels)
  mutrel.check_posterior_sanity(posterior.rels)
  assert np.all(np.isclose(1, np.sum(posterior.rels, axis=2)))

  # TODO: only calculate posterior once here, instead of computing it within
  # each worker separately for a given variant pair.
  other = _calc_posterior_full(evidence.rels, logprior)
  assert np.allclose(posterior.rels, other)
  return (posterior, evidence)

def calc_posterior(variants, logprior, rel_type, parallel=1):
  M = len(variants)
  # Allow Numba use by converting to namedtuple.
  vids = common.extract_vids(variants)
  variants = [common.convert_variant_dict_to_tuple(variants[V]) for V in vids]

  posterior = mutrel.init_mutrel(vids)
  evidence = mutrel.init_mutrel(vids)
  pairs = list(itertools.combinations(range(M), 2)) + [(V, V) for V in range(M)]

  _compute = lambda pbar: _compute_pairs(
     pairs,
     variants,
     logprior,
     posterior,
     evidence,
     pbar,
     parallel,
  )

  if parallel > 0:
    with tqdm.tqdm(total=len(pairs), desc='Computing %s relations' % rel_type, unit='pair', dynamic_ncols=True) as pbar:
      posterior, evidence =_compute(pbar)
  else:
    posterior, evidence =_compute(None)
  return (posterior, evidence)

def merge_variants(to_merge, evidence, logprior):
  assert np.all(np.array([V for group in to_merge for V in group]) < len(evidence.vids))
  already_merged = set()

  for vidxs in to_merge:
    vidxs = set(vidxs)
    assert len(vidxs & already_merged) == 0

    M_old = len(evidence.vids)
    merged_vid = ','.join([evidence.vids[V] for V in vidxs])
    new_vids = evidence.vids + [merged_vid]

    new_evidence = mutrel.init_mutrel(new_vids)
    new_evidence.rels[:-1,:-1] = evidence.rels

    merged_row = np.sum(np.array([evidence.rels[V] for V in vidxs]), axis=0)
    assert merged_row.shape == (M_old, NUM_MODELS)
    merged_col = np.copy(merged_row)
    merged_col[:,Models.A_B] = merged_row[:,Models.B_A]
    merged_col[:,Models.B_A] = merged_row[:,Models.A_B]

    new_evidence.rels[-1,:-1] = merged_row
    new_evidence.rels[:-1,-1] = merged_col
    new_evidence.rels[-1,-1,:] = -np.inf
    new_evidence.rels[-1,-1,Models.cocluster] = 0

    already_merged |= vidxs
    evidence = new_evidence

  evidence = mutrel.remove_variants_by_vidx(evidence, already_merged)
  posterior = make_full_posterior(evidence, logprior)
  return (posterior, evidence)

def add_variants(vids_to_add, variants, mutrel_posterior, mutrel_evidence, logprior, pbar, parallel):
  for vid in vids_to_add:
    assert vid in variants

  new_vids = mutrel_posterior.vids + vids_to_add
  # M: number of variants we have now
  # A: number of variants we added
  M = len(new_vids)
  A = len(vids_to_add)
  assert len(mutrel_posterior.vids) == len(mutrel_evidence.vids) == M - A
  variants = [common.convert_variant_dict_to_tuple(variants[V]) for V in new_vids]

  new_posterior = mutrel.init_mutrel(new_vids)
  new_evidence = mutrel.init_mutrel(new_vids)
  new_posterior.rels[:-A,:-A] = mutrel_posterior.rels
  new_evidence.rels[:-A,:-A] = mutrel_evidence.rels

  pairs = [(I, J) for I in range(M) for J in range(M - A, M)]

  return _compute_pairs(
    pairs,
    variants,
    logprior,
    new_posterior,
    new_evidence,
    pbar,
    parallel,
  )

def _calc_lh_and_posterior(V1, V2, logprior):
  evidence, evidence_per_sample = lh.calc_lh(V1, V2)
  posterior = _calc_posterior(evidence, logprior)
  return (evidence, posterior)

def _examine(V1, V2, variants, logprior=None, _calc_lh=None):
  E, Es = lh.calc_lh(*[common.convert_variant_dict_to_tuple(V) for V in (variants[V1], variants[V2])], _calc_lh)
  Es -= np.max(Es, axis=1)[:,None]
  sep = np.nan * np.ones(len(variants[V1]['var_reads']))[:,None]
  persamp = np.hstack((
    variants[V1]['var_reads'][:,None],
    variants[V1]['total_reads'][:,None],
    variants[V1]['vaf'][:,None],
    variants[V1]['omega_v'][:,None],
    variants[V1]['vaf'][:,None] / variants[V1]['omega_v'][:,None],
    sep,
    variants[V2]['var_reads'][:,None],
    variants[V2]['total_reads'][:,None],
    variants[V2]['vaf'][:,None],
    variants[V2]['omega_v'][:,None],
    variants[V2]['vaf'][:,None] / variants[V2]['omega_v'][:,None],
    sep,
    Es,
  ))

  post = _calc_posterior(E, _complete_logprior(logprior))
  return (persamp, E, post)
