import numpy as np
import scipy.stats
import scipy.special
import concurrent.futures
from tqdm import tqdm
import itertools

from common import Models
import common
import lh
import util

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

def calc_posterior(variants, prior=None, parallel=1):
  logprior = complete_prior(prior)
  N = len(variants)
  posterior = {}
  evidence = {}

  # Allow Numba use by converting to namedtuple.
  variants = {vid: common.convert_variant_dict_to_tuple(V) for vid, V in variants.items()}

  combos = list(itertools.combinations(variants.keys(), 2)) + [(K,K) for K in variants.keys()]
  combos = [sorted(C) for C in combos]
  # Don't bother starting more workers than jobs.
  parallel = min(parallel, len(combos))
  # Don't bother invoking all the parallelism machinery if we're only fitting a
  # single sample at a time.
  if parallel > 1:
    pairs = []
    futures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=parallel) as ex:
      for vidx1, vidx2 in combos:
        pairs.append((vidx1, vidx2))
        futures.append(ex.submit(calc_lh, variants[vidx1], variants[vidx2]))
      for F in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc='Computing relations', unit=' pairs', dynamic_ncols=True):
        pass
    for pair, F in zip(pairs, futures):
      evidence[pair], _ = F.result()
  else:
    for vidx1, vidx2 in combos:
      pair = (vidx1, vidx2)
      evidence[pair], _ = calc_lh(variants[vidx1], variants[vidx2])

  # Duplicate evidence's keys, since we'll be modifying dictionary.
  for A, B in list(evidence.keys()):
    if A == B:
      continue
    evidence[(B,A)] = swap_evidence(evidence[(A,B)])
  for pair in evidence.keys():
    posterior[pair] = _calc_posterior(evidence[pair], logprior)
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
