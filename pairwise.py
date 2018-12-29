import numpy as np
import scipy.special
import concurrent.futures
from progressbar import progressbar
import itertools

from common import Models
import common
import lh
import mutrel

def swap_A_B(arr):
  swapped = np.zeros(len(arr)) + np.nan
  for midx, M in enumerate(Models._all):
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

def _calc_posterior_full(evidence, prior):
  logprior = _complete_prior(prior)
  # This function is currently unused.
  joint = evidence + logprior[None,None,:]
  B = np.max(joint, axis=2)
  joint -= B[:,:,None]
  expjoint = np.exp(joint)
  posterior = expjoint / np.sum(expjoint, axis=2)[:,:,None]

  # Regardless of what the prior is on coclustering (even if it's zero), if
  # we're dealing with a pair consisting of variant `i` and variant `i`, we
  # want the posterior to indicate coclustering certainty.
  M = range(len(posterior))
  posterior[M,M,:] = 0
  posterior[M,M,Models.cocluster] = 1
  mutrel.check_posterior_sanity(posterior)

  return posterior

def _complete_prior(prior):
  '''
  If you don't want to include a certain model in the posterior (e.g., perhaps
  you're computing pairwise probabilities for supervariants, and so don't want
  to allow them to cocluster unless they're the same supervariant), set the
  corresponding entry in prior to 0.
  '''
  if prior is None:
    prior = {}
  for K in prior.keys():
    assert K in Models._all
    assert 0 <= prior[K] <= 1

  total = np.sum(list(prior.values()))
  assert 0 <= total <= 1
  remaining = set(Models._all) - set(prior.keys())
  for K in remaining:
    prior[K] = (1 - total) / len(remaining)

  logprior = -np.inf * np.ones(len(prior))
  for K in prior:
    if prior[K] == 0:
      continue
    logprior[getattr(Models, K)] = np.log(prior[K])

  assert np.isclose(0, scipy.special.logsumexp(logprior))
  return logprior

def _compute_pairs(pairs, variants, prior, posterior, evidence, pbar=None, parallel=1):
  logprior = _complete_prior(prior)
  # TODO: change ordering of pairs based on what will provide optimal
  # integration accuracy according to Quaid's advice.
  pairs = [sorted(C) for C in pairs]
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
  assert np.allclose(posterior.rels, _calc_posterior_full(evidence.rels, prior))
  return (posterior, evidence)

def calc_posterior(variants, prior, rel_type, parallel=1):
  M = len(variants)
  # Allow Numba use by converting to namedtuple.
  vids = common.extract_vids(variants)
  variants = [common.convert_variant_dict_to_tuple(variants[V]) for V in vids]

  posterior = mutrel.init_mutrel(vids)
  evidence = mutrel.init_mutrel(vids)
  pairs = list(itertools.combinations(range(M), 2)) + [(V, V) for V in range(M)]

  with progressbar(total=len(pairs), desc='Computing %s relations' % rel_type, unit='pair', dynamic_ncols=True) as pbar:
    return _compute_pairs(
      pairs,
      variants,
      prior,
      posterior,
      evidence,
      pbar = pbar,
      parallel = parallel,
    )

def merge_variants(to_merge, evidence, prior):
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
    assert merged_row.shape == (M_old, len(Models._all))
    merged_col = np.copy(merged_row)
    merged_col[:,Models.A_B] = merged_row[:,Models.B_A]
    merged_col[:,Models.B_A] = merged_row[:,Models.A_B]

    new_evidence.rels[-1,:-1] = merged_row
    new_evidence.rels[:-1,-1] = merged_col
    new_evidence.rels[-1,-1,:] = -np.inf
    new_evidence.rels[-1,-1,Models.cocluster] = 0

    already_merged |= vidxs
    evidence = new_evidence

  evidence = mutrel.remove_variants(evidence, already_merged)
  posterior = mutrel.Mutrel(
    vids = evidence.vids,
    rels = _calc_posterior_full(evidence.rels, prior),
  )
  return (posterior, evidence)

def add_variants(vids_to_add, variants, mutrel_posterior, mutrel_evidence, prior, pbar, parallel):
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
    prior,
    new_posterior,
    new_evidence,
    pbar,
    parallel,
  )

def _calc_lh_and_posterior(V1, V2, logprior):
  evidence, evidence_per_sample = lh.calc_lh(V1, V2)
  posterior = _calc_posterior(evidence, logprior)
  return (evidence, posterior)

def _examine(V1, V2, variants, _calc_lh=None):
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

  prior = {'garbage': 0.001}
  post1 = _calc_posterior(E, _complete_prior(None))
  post2 = _calc_posterior(E, _complete_prior({'garbage': 0.001}))
  return (persamp, E, post1, post2)
