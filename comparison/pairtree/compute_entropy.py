#!/usr/bin/env python
import argparse
import numpy as np
import numpy.ma as ma

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import util
import evalutil

from collections import defaultdict

def compute_parent_dist(structs, weights):
  K = len(structs[0]) + 1
  parent_dist = np.zeros((K, K))
  assert np.all(weights >= 0) and np.isclose(1, np.sum(weights))

  for struct, weight in zip(structs, weights):
    adjm = util.convert_parents_to_adjmatrix(struct)
    np.fill_diagonal(adjm, 0)
    parent_dist += weight * adjm

  assert np.all(0 == parent_dist[:,0])
  parent_dist = parent_dist[:,1:]
  parent_dist = evalutil.fix_rounding_errors(parent_dist)

  assert np.all(0 <= parent_dist) and np.all(parent_dist <= 1)
  assert np.allclose(1, np.sum(parent_dist, axis=0))
  return parent_dist

def compute_parentropy(parent_dist):
  K = len(parent_dist)
  assert parent_dist.shape == (K, K-1)
  parent_dist = ma.masked_equal(parent_dist, 0)
  parent_entropy = -ma.sum(parent_dist * ma.log2(parent_dist), axis=0)
  total_entropy = np.sum(parent_entropy)
  entropy_per_node = total_entropy / (K-1)
  return (total_entropy, entropy_per_node, parent_entropy)

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

def compute_tree_jsd(structs1, probs1, structs2, probs2):
  '''Compute Jensen-Shannon divergence.'''
  M = defaultdict(lambda: 0)

  def _dostuff(structs, probs):
    assert len(structs) == len(probs)
    T_probs = {}
    for struct, prob in zip(structs, probs):
      # Sometimes LLH is so negative that the softmax of it ends up being
      # actually zero. Pretend that such trees weren't sampled.
      # Actually, since we use the halved-weight to create the mixture, it's this
      # value that needs to be nonzero.
      if 0.5*prob == 0:
        continue
      key = hash(struct.tobytes())
      M[key] += 0.5 * prob
      assert key not in T_probs
      T_probs[key] = prob
    return T_probs

  T_probs1 = _dostuff(structs1, probs1)
  T_probs2 = _dostuff(structs2, probs2)
  for arr in (M, T_probs1, T_probs2):
    V = np.array(list(arr.values()))
    assert np.isclose(1, np.sum(V)) and np.all(V >= 0)

  kld_1 = compute_tree_kld(T_probs1, M)
  kld_2 = compute_tree_kld(T_probs2, M)
  jsd = 0.5*(kld_1 + kld_2)
  if jsd > 1:
    assert np.isclose(1, jsd)
    jsd = 1
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
  sampled = set([T.tobytes() for T in sampled])
  truth = set([T.tobytes() for T in truth])
  truth_recovered = len(sampled & truth) / len(truth)
  jaccard = len(sampled & truth) / len(sampled | truth)
  return (truth_recovered, jaccard)

def make_sorted(arr, N=10):
  arr = np.array(sorted(arr, reverse=True))[:N]
  return np.array2string(arr, suppress_small=True, separator=' ')

def process(results, truth):
  K = len(results['struct'][0])

  probs1 = results['count'] / np.sum(results['count'])
  probs2 = util.softmax(results['llh'])
  probs3 = util.softmax(results['llh'] + np.log(results['count']))
  # Verify that the means by which we compute posterior probabilities in the
  # results files hasn't changed. (Even if it has, we don't use
  # `results['prob']` in this file, so it should be fine.)
  assert np.allclose(probs3, results['prob'])

  pard = compute_parent_dist(results['struct'], probs3)
  parentropy, _, _ = compute_parentropy(pard)

  truth_num = len(truth['struct'])
  truth_llh = np.zeros(truth_num)
  truth_probs = np.ones(truth_num) / truth_num
  truth_pard = compute_parent_dist(truth['struct'], truth_probs)
  truth_parentropy, _, _ = compute_parentropy(truth_pard)
  assert np.allclose(truth_probs, truth['prob'])

  top_probs = 10
  good_thresh = 1e-3
  jsd_parents = compute_parents_jsd(pard, truth_pard)
  jsd_parents_phi_mean = jsd_parents * np.mean(truth['phi'][0,1:], axis=1)
  jsd_parents_phi_max  = jsd_parents * np.max(truth['phi'][0,1:], axis=1)

  stats = {}
  stats['true_trees'] = truth_num
  stats['sampled_unique_trees'] = len(probs3)
  stats['num_good'] = np.sum(probs3 >= good_thresh)
  stats['prop_good'] = '%.3f' % (np.sum(probs3 >= good_thresh) / len(probs3))
  stats['H_trees_truth'] = calc_entropy(truth_probs)
  stats['H_trees_pairtree_1'] = calc_entropy(probs1)
  stats['H_trees_pairtree_2'] = calc_entropy(probs2)
  stats['H_trees_pairtree_3'] = calc_entropy(probs3)
  stats['H_parents_truth'] = truth_parentropy
  stats['H_parents_pairtree'] = parentropy
  stats['prop_truth_recovered'], stats['jaccard'] = compute_indices(results['struct'], truth['struct'])
  stats['jsd_trees'] = compute_tree_jsd(results['struct'], probs3, truth['struct'], truth_probs)
  stats['jsd_parents_sum'] = np.sum(jsd_parents)
  stats['jsd_parents_mean'] = np.sum(jsd_parents) / K
  stats['jsd_parents_max'] = np.max(jsd_parents)
  stats['jsd_parents'] = jsd_parents
  stats['jsd_parents_phi_mean'] = np.max(jsd_parents_phi_mean)
  stats['jsd_parents_phi_max'] = np.max(jsd_parents_phi_max)
  stats['jsd_parents_phi_mean_top10'] = make_sorted(jsd_parents_phi_mean)
  stats['jsd_parents_phi_max_top10']  = make_sorted(jsd_parents_phi_max)
  stats['top_probs_1_top10'] = make_sorted(probs1)
  stats['top_probs_2_top10'] = make_sorted(probs2)
  stats['top_probs_3_top10'] = make_sorted(probs3)

  keys = list(stats.keys())
  vals = [stats[key] for key in keys]
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
