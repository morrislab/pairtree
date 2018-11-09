import numpy as np
import scipy.special
import concurrent.futures
from tqdm import tqdm
import itertools

from common import Models
import common
import lh

def swap_evidence(evidence):
  swapped = np.zeros(len(evidence)) + np.nan
  for midx, M in enumerate(Models._all):
    if midx in (Models.garbage, Models.cocluster, Models.diff_branches):
      swapped[midx] = evidence[midx]
    elif midx == Models.A_B:
      swapped[midx] = evidence[Models.B_A]
    elif midx == Models.B_A:
      swapped[midx] = evidence[Models.A_B]
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

def complete_prior(prior):
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

def _sanity_check_tensor(tensor):
  assert not np.any(np.isnan(tensor))
  for model in ('garbage', 'cocluster', 'diff_branches'):
    # These should be symmetric.
    midx = getattr(common.Models, model)
    mat = tensor[:,:,midx]
    assert np.allclose(mat, mat.T)

def calc_posterior(variants, prior=None, parallel=1):
  logprior = complete_prior(prior)
  make_tqdm = lambda I, I_len: tqdm(
    I,
    total=I_len,
    desc='Computing relations',
    unit='pair',
    dynamic_ncols=True,
    disable=None # Disable on non-TTY
  )

  M = len(variants)
  # Allow Numba use by converting to namedtuple.
  variants = [common.convert_variant_dict_to_tuple(V) for V in variants]

  evidence = np.nan * np.ones((M, M, len(Models._all)))
  posterior = np.copy(evidence)

  pairs = list(itertools.combinations(range(M), 2)) + [(V, V) for V in range(M)]
  # TODO: change ordering of pairs based on what will provide optimal
  # integration accuracy according to Quaid's advice.
  pairs = [sorted(C) for C in pairs]
  # Don't bother starting more workers than jobs.
  parallel = min(parallel, len(pairs))
  # Don't bother invoking all the parallelism machinery if we're only fitting a
  # single sample at a time.
  if parallel > 1:
    futures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=parallel) as ex:
      for A, B in pairs:
        futures.append(ex.submit(calc_lh, variants[A], variants[B]))
      for F in make_tqdm(concurrent.futures.as_completed(futures), len(futures)):
        pass
    for (A, B), F in zip(pairs, futures):
      evidence[A,B], _ = F.result()
  else:
    for A, B in make_tqdm(pairs, len(pairs)):
      evidence[A,B], _ = calc_lh(variants[A], variants[B])

  # Duplicate evidence's keys, since we'll be modifying dictionary.
  for A, B in pairs:
    if A == B:
      continue
    evidence[B,A] = swap_evidence(evidence[A,B])
  for A in range(M):
    for B in range(M):
      # TODO: move this into the parallel calcs.
      posterior[A,B] = _calc_posterior(evidence[A,B], logprior)

  _sanity_check_tensor(evidence)
  _sanity_check_tensor(posterior)
  return (posterior, evidence)

# In cases where we have zero variant reads for both variants, the garbage
# models ends up taking considerable posterior mass, but we have no information
# to judge whether the pair should be garbage. The contribution from lots of
# samples with zero variant reads on both variants can end up overwhelming
# informative samples with non-zero variant reads, such that the pair across
# samples jointly is deemed garbage. This is undesireable behaviour, clearly,
# and so ignore such (0, 0) cases by zeroing out their evidence. This
# effectively states that we have a uniform posterior over the relations of
# those mutations in that sample, which is the sensible thing to do.
#
# See this discussion on Slack for more details:
# https://morrislab.slack.com/archives/CCE5HNVSP/p1540392641000100
def _fix_zero_var_read_samples(V1, V2, evidence):
  zero_indices = np.logical_and(V1.var_reads == 0, V2.var_reads == 0)
  evidence[zero_indices,:] = 0

def calc_lh(V1, V2, _calc_lh=lh.calc_lh_quad):
  if V1.id == V2.id:
    # If they're the same variant, they should cocluster with certainty.
    evidence = -np.inf * np.ones(len(Models._all))
    evidence[Models.cocluster] = 0
    return (evidence, evidence)

  evidence_per_sample = _calc_lh(V1, V2) # SxM
  evidence_per_sample[:,Models.garbage] = lh.calc_garbage(V1, V2)
  _fix_zero_var_read_samples(V1, V2, evidence_per_sample)
  evidence = np.sum(evidence_per_sample, axis=0)
  return (evidence, evidence_per_sample)
