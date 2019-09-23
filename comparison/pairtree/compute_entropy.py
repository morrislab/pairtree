#!/usr/bin/env python
import argparse
import numpy as np
import numpy.ma as ma
import scipy.special

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import resultserializer
import util
import pardist

from collections import defaultdict

def resolve_unique(adjms, llhs):
  unique = {}

  for adjm, llh in zip(adjms, llhs):
    key = hash(adjm.tobytes())
    if key in unique:
      assert np.isclose(unique[key]['llh'], llh)
      assert np.array_equal(unique[key]['adjm'], adjm)
      unique[key]['count'] += 1
    else:
      unique[key] = {
        'adjm': adjm,
        'llh': llh,
        'count': 1,
      }

  merged = sorted(unique.values(), key = lambda R: -R['llh'])
  return merged

def calc_entropy(A):
  A = ma.masked_equal(A, 0)
  ent = -ma.sum(A * ma.log2(A))
  return np.abs(np.array(ent))

def compute_tree_kld(P, Q):
  '''Compute KL(P || Q).'''
  assert set(P.keys()).issubset(Q.keys())
  keys = list(P.keys())
  Parr = np.array([P[K] for K in keys])
  Qarr = np.array([Q[K] for K in keys])
  P_nonzero = Parr != 0
  assert np.all(Qarr[P_nonzero] > 0)

  logP = np.log2(Parr[P_nonzero])
  logQ = np.log2(Qarr[P_nonzero])
  kld = np.sum(Parr[P_nonzero] * (logP - logQ))
  assert kld >= 0
  return kld

def compute_tree_jsd(adjms1, probs1, adjms2, probs2):
  '''Compute Jensen-Shannon divergence.'''
  M = defaultdict(lambda: 0)

  def _dostuff(A, P):
    assert len(A) == len(P)
    A_probs = {}
    for _A, _P in zip(A, P):
      # Sometimes LLH is so negative that the softmax of it ends up being
      # actually zero. Pretend that such trees weren't sampled.
      # Actually, since we use the halved-weight to create the mixture, it's this
      # value that needs to be nonzero.
      if 0.5*_P == 0:
        continue
      key = hash(util.find_parents(_A).tobytes())
      M[key] += 0.5 * _P
      assert key not in A_probs
      A_probs[key] = _P
    return A_probs

  A_probs1 = _dostuff(adjms1, probs1)
  A_probs2 = _dostuff(adjms2, probs2)
  for arr in (M, A_probs1, A_probs2):
    V = np.array(list(arr.values()))
    assert np.isclose(1, np.sum(V)) and np.all(V >= 0)

  kld_1 = compute_tree_kld(A_probs1, M)
  kld_2 = compute_tree_kld(A_probs2, M)
  jsd = 0.5*(kld_1 + kld_2)
  # These bounds are guaranteed by JSD (when it's computed in bits).
  assert 0 <= jsd <= 1
  return jsd

def compute_parents_kld(P, Q):
  for A in (P, Q):
    assert np.all(A >= 0)
    assert np.allclose(1, np.sum(A, axis=0))
  logP = ma.log2(ma.masked_equal(P, 0))
  logQ = ma.log2(ma.masked_equal(Q, 0))
  kld = np.sum(P * (logP - logQ), axis=0)

  assert np.allclose(0, kld[kld < 0])
  kld = np.abs(kld)
  assert np.all(kld >= 0)
  return kld

def compute_parents_jsd(pardist1, pardist2):
  M = 0.5*(pardist1 + pardist2)
  kld1 = compute_parents_kld(pardist1, M)
  kld2 = compute_parents_kld(pardist2, M)
  jsd = 0.5*(kld1 + kld2)
  assert np.allclose(1, jsd[jsd > 1])
  jsd = np.minimum(1, jsd)
  assert np.all(jsd >= 0) and np.all(jsd <= 1)
  return jsd

def compute_indices(sampled, truth):
  sampled = set([util.find_parents(A).tobytes() for A in sampled])
  truth = set([util.find_parents(A).tobytes() for A in truth])
  truth_recovered = len(sampled & truth) / len(truth)
  jaccard = len(sampled & truth) / len(sampled | truth)
  return (truth_recovered, jaccard)

def make_sorted(arr, N=10):
  arr = np.array(sorted(arr, reverse=True))[:N]
  return np.array2string(arr, suppress_small=True, separator=' ')

def process(results, truth):
  merged = resolve_unique(results['adjm'], results['llh'])
  unique_adjm = [M['adjm'] for M in merged]
  K = len(unique_adjm[0]) - 1

  pard = pardist.compute_parent_dist(results['adjm'], results['llh'])
  parentropy, _, _ = pardist.compute_entropy(pard)

  llhs   = np.array([M['llh']   for M in merged])
  counts = np.array([M['count'] for M in merged])
  probs1 = counts / np.sum(counts)
  probs2 = util.softmax(llhs)
  probs3 = util.softmax(llhs + np.log(counts))

  truth_num = len(truth['adjm'])
  truth_llh = np.zeros(truth_num)
  truth_probs = np.ones(truth_num) / truth_num
  truth_pard = pardist.compute_parent_dist(truth['adjm'], truth_llh)
  truth_parentropy, _, _ = pardist.compute_entropy(truth_pard)

  top_probs = 10
  good_thresh = 1e-3
  jsd_parents = compute_parents_jsd(pard, truth_pard)
  jsd_parents_phi_mean = jsd_parents * np.mean(truth['phi'][0,1:], axis=1)
  jsd_parents_phi_max  = jsd_parents * np.max(truth['phi'][0,1:], axis=1)

  results = {}
  results['true_trees'] = truth_num
  results['sampled_unique_trees'] = len(probs3)
  results['num_good'] = np.sum(probs3 >= good_thresh)
  results['prop_good'] = '%.3f' % (np.sum(probs3 >= good_thresh) / len(probs3))
  results['H_trees_truth'] = calc_entropy(truth_probs)
  results['H_trees_pairtree_1'] = calc_entropy(probs1)
  results['H_trees_pairtree_2'] = calc_entropy(probs2)
  results['H_trees_pairtree_3'] = calc_entropy(probs3)
  results['H_parents_truth'] = truth_parentropy
  results['H_parents_pairtree'] = parentropy
  results['prop_truth_recovered'], results['jaccard'] = compute_indices(unique_adjm, truth['adjm'])
  results['jsd_trees'] = compute_tree_jsd(unique_adjm, probs3, truth['adjm'], truth_probs)
  results['jsd_parents_sum'] = np.sum(jsd_parents)
  results['jsd_parents_mean'] = np.sum(jsd_parents) / K
  results['jsd_parents_max'] = np.max(jsd_parents)
  results['jsd_parents'] = jsd_parents
  results['jsd_parents_phi_mean'] = np.max(jsd_parents_phi_mean)
  results['jsd_parents_phi_max'] = np.max(jsd_parents_phi_max)
  results['jsd_parents_phi_mean_top10'] = make_sorted(jsd_parents_phi_mean)
  results['jsd_parents_phi_max_top10']  = make_sorted(jsd_parents_phi_max)
  results['top_probs_1_top10'] = make_sorted(probs1)
  results['top_probs_2_top10'] = make_sorted(probs2)
  results['top_probs_3_top10'] = make_sorted(probs3)

  keys = list(results.keys())
  vals = [results[key] for key in keys]
  for A in (keys, vals):
    print(*A, sep=',')

def main():
  np.set_printoptions(linewidth=400, precision=3, threshold=sys.maxsize, suppress=True)
  np.seterr(divide='raise', invalid='raise')

  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--truth', dest='truth_fn')
  parser.add_argument('results_fn')
  args = parser.parse_args()

  results = np.load(args.results_fn)
  if args.truth_fn is not None:
    truth = np.load(args.truth_fn)
  else:
    truth = None
  process(results, truth)

if __name__ == '__main__':
  main()
