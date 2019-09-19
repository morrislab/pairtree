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
      key = hash(util.find_parents(_A).tobytes())
      M[key] += 0.5 * _P
      assert key not in A_probs
      A_probs[key] = _P
    return A_probs

  A_probs1 = _dostuff(adjms1, probs1)
  A_probs2 = _dostuff(adjms2, probs2)
  a1 = set(A_probs1)
  a2 = set(A_probs2)

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
  assert np.all(jsd >= 0) and np.all(jsd <= 1)
  joint_jsd = np.sum(jsd)
  return joint_jsd

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
  results['jsd_trees'] = compute_tree_jsd(unique_adjm, probs3, truth['adjm'], truth_probs)
  results['jsd_parents'] = compute_parents_jsd(pard, truth_pard)
  results['jsd_parents_per_node'] = compute_parents_jsd(pard, truth_pard) / K
  results['top_probs_1'] = np.array(sorted(probs1, reverse=True)[:top_probs])
  results['top_probs_2'] = np.array(sorted(probs2, reverse=True)[:top_probs])
  results['top_probs_3'] = np.array(sorted(probs3, reverse=True)[:top_probs])

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

  results = resultserializer.load(args.results_fn)
  if args.truth_fn is not None:
    truth = resultserializer.load(args.truth_fn)
  else:
    truth = None
  process(results, truth)

if __name__ == '__main__':
  main()
