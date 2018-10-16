import numpy as np
import scipy.stats
import concurrent.futures
from tqdm import tqdm
import itertools

from common import Models
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

def _calc_posterior(evidence, include_garbage, include_cocluster):
  evidence = np.copy(evidence)
  if not include_garbage:
    evidence[Models.garbage] = -np.inf
  if not include_cocluster:
    # We still want the posterior for coclustering of node `i` with `i` to be
    # certain. To signal that we're computing the (i, i) pairwise posterior,
    # the evidence for coclustering will be set to 0, and all other entries
    # will be -inf. In this instance, don't set the coclustering evidence to
    # -inf -- instead, leave it at 0 so that the posterior will indicate
    # certainty about coclustering.
    if evidence[Models.cocluster] < 0:
      evidence[Models.cocluster] = -np.inf
  B = np.max(evidence)
  evidence -= B
  posterior = np.exp(evidence) / np.sum(np.exp(evidence))
  return posterior

def calc_posterior(variants, parallel=1, include_garbage_in_posterior=False, include_cocluster_in_posterior=False):
  N = len(variants)
  posterior = {}
  evidence = {}

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
    posterior[pair] = _calc_posterior(evidence[pair], include_garbage_in_posterior, include_cocluster_in_posterior)
  return (posterior, evidence)

def calc_lh(V1, V2, _calc_lh=lh.calc_lh_quad):
  if V1['id'] == V2['id']:
    # If they're the same variant, they should cocluster with certainty.
    evidence = -np.inf * np.ones(len(Models._all))
    evidence[Models.cocluster] = 0
    return (evidence, evidence)

  evidence_per_sample = _calc_lh(V1, V2)
  garb1, garb2 = [scipy.stats.randint.logpmf(V['var_reads'], 0, V['total_reads'] + 1) for V in (V1, V2)]
  evidence_per_sample[:,Models.garbage] = garb1 + garb2
  evidence = np.sum(evidence_per_sample, axis=0)
  return (evidence, evidence_per_sample)

def test_calc_lh(V1, V2):
  inputs = {
    'vaf': np.array([V['vaf'] for V in (V1, V2)]),
    'var_reads': np.array([V['var_reads'] for V in (V1, V2)]),
    'total_reads': np.array([V['total_reads'] for V in (V1, V2)]),
  }
  records = {'evidence': {}, 'norm_evidence': {}, 'posterior': {}}

  for M in (
    lh.calc_lh_quad,
    lh.calc_lh_mc_1D,
    lh.calc_lh_mc_2D,
    lh.calc_lh_grid,
  ):
    M_name = M.__name__
    M = util.time_exec(M)
    evidence, evidence_per_sample = calc_lh(V1, V2, _calc_lh=M)
    records['evidence'][M_name] = evidence
    records['norm_evidence'][M_name] = evidence - np.max(evidence)
    records['posterior'][M_name] = _calc_posterior(evidence, include_garbage=True, include_cocluster=False)

    sep = np.zeros(len(evidence_per_sample))[:,None] + np.nan
    combined = np.hstack((
      inputs['vaf'].T,
      sep,
      inputs['var_reads'].T,
      sep,
      inputs['total_reads'].T,
      sep,
      evidence_per_sample,
      sep,
      evidence_per_sample - np.max(evidence_per_sample, axis=1)[:,None],
    ))
    print(M_name, '(%.3f ms)' % util.time_exec._ms)
    print(combined)
    for J in sorted(records.keys()):
      print(J, records[J][M_name], sep='\t')
    print()

  for K in inputs.keys():
    print(K)
    print(inputs[K])
  print()

  for J in sorted(records.keys()):
    print(J)
    for K in sorted(records[J].keys()):
      print(K, records[J][K], sep='\t')
    print()

  import sys
  sys.exit()

def generate_results(posterior, evidence, variants):
  model_probs = {}
  model_evidence = {}
  for midx, model in enumerate(Models._all):
    model_probs[model]    = {'%s,%s' % (vid1, vid2): P[midx] for (vid1, vid2), P in posterior.items() }
    model_evidence[model] = {'%s,%s' % (vid1, vid2): P[midx] for (vid1, vid2), P in evidence.items()  }

  results = {
    'models': Models._all,
    'model_probs': model_probs,
    'model_evidence': model_evidence,
    'variants': { V: {'name': variants[V]['name']} for V in variants.keys() },
  }
  return results
