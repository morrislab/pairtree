import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from common import Models
import mutrel
import resultserializer
import evalutil

def calc_mutrel_from_trees(adjms, llhs, clusters, tree_weights):
  clusterings = [clusters for idx in range(len(adjms))]
  return evalutil.calc_mutrel_from_trees(adjms, llhs, clusterings, tree_weights)

def calc_mutrel_from_clustrel(clustrel, clusters):
  mutrel.check_posterior_sanity(clustrel.rels)
  assert len(clusters[0]) == 0
  vids, membership = evalutil.make_membership_mat(clusters[1:])
  # K: number of non-empty clusters
  M, K = membership.shape

  num_models = len(Models._all)
  mrel = np.zeros((M, M, num_models))
  assert clustrel.rels.shape == (K, K, num_models)

  for modelidx in range(num_models):
    mut_vs_cluster = np.dot(membership, clustrel.rels[:,:,modelidx]) # MxK
    mrel[:,:,modelidx] = np.dot(mut_vs_cluster, membership.T)
  mutrel.check_posterior_sanity(mrel)

  return mutrel.Mutrel(
    vids = vids,
    rels = mrel,
  )

def choose_subset(results, subset_size):
  keys = ('adjm', 'llh', 'phi')
  lengths = set([len(results[K]) for K in keys])
  assert len(lengths) == 1
  total = lengths.pop()
  assert 0 < subset_size <= total
  assert total % subset_size == 0

  partitions = int(total / subset_size)
  part_idx = np.random.randint(partitions)
  start  = part_idx * subset_size
  end = start + subset_size
  for K in keys:
    results[K] = results[K][start:end]

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--use-subset', type=int)
  parser.add_argument('--weight-trees-by', choices=('llh', 'uniform'), default='uniform')
  parser.add_argument('--clustrel-mutrel')
  parser.add_argument('--trees-mutrel')
  parser.add_argument('--pure-mutrel')
  parser.add_argument('--phi', dest='mutphifn')
  parser.add_argument('pairtree_results_fn')
  args = parser.parse_args()

  results = resultserializer.load(args.pairtree_results_fn)
  clusters = [[]] + list(results['clusters'])
  garbage = list(results['garbage'])
  all_vids = set([V for C in results['clusters'] for V in C] + garbage)

  if args.use_subset is not None:
    choose_subset(results, args.use_subset)

  if args.pure_mutrel is not None:
    pure_mutrel = results['mutrel_posterior']
    assert set(pure_mutrel.vids) == all_vids
    evalutil.save_sorted_mutrel(pure_mutrel, args.pure_mutrel)

  if args.trees_mutrel is not None:
    tree_mutrel = calc_mutrel_from_trees(results['adjm'], results['llh'], clusters, args.weight_trees_by)
    tree_mutrel = evalutil.add_garbage(tree_mutrel, garbage)
    assert set(tree_mutrel.vids) == all_vids
    evalutil.save_sorted_mutrel(tree_mutrel, args.trees_mutrel)

  if args.clustrel_mutrel is not None:
    clustrel_mutrel = calc_mutrel_from_clustrel(results['clustrel_posterior'], clusters)
    clustrel_mutrel = evalutil.add_garbage(clustrel_mutrel, garbage)
    assert set(clustrel_mutrel.vids) == all_vids
    evalutil.save_sorted_mutrel(clustrel_mutrel, args.clustrel_mutrel)

  if args.mutphifn is not None:
    clusterings = [clusters for idx in range(len(results['adjm']))]
    mutphi = evalutil.calc_mutphi(results['phi'], results['llh'], clusterings, args.weight_trees_by)
    evalutil.save_mutphi(mutphi, args.mutphifn)

if __name__ == '__main__':
  main()
