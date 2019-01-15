import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import common
import mutrel
import resultserializer
import evalutil

def calc_mutrel_from_trees(adjms, llhs, clusters, tree_weights):
  if tree_weights == 'llh':
    weights = evalutil.softmax(llhs)
  elif tree_weights == 'uniform':
    weights = np.ones(len(adjms)) / len(adjms)
  else:
    raise Exception('Unknown tree_weights=%s' % tree_weights)
  vids = None

  for adjm, weight in zip(adjms, weights):
    mrel = mutrel.make_mutrel_tensor_from_cluster_adj(adjm, clusters)
    if vids is None:
      vids = mrel.vids
      soft_mutrel = np.zeros(mrel.rels.shape)
    else:
      assert mrel.vids == vids
    soft_mutrel += weight * mrel.rels

  soft_mutrel = evalutil.fix_rounding_errors(soft_mutrel)
  return mutrel.Mutrel(
    vids = vids,
    rels = soft_mutrel,
  )

def make_membership_mat(clusters):
  vids = common.sort_vids([vid for C in clusters for vid in C])
  vidmap = {vid: vidx for vidx, vid in enumerate(vids)}
  N = len(vids)
  K = len(clusters)

  # membership[i,j] = 1 iff mutation `i` is in cluster `j`
  membership = np.zeros((N, K))
  for cidx, C in enumerate(clusters):
    members = [vidmap[vid] for vid in C]
    membership[members,cidx] = 1
  return (vids, membership)

def calc_mutrel_from_clustrel(clustrel, clusters):
  mutrel.check_posterior_sanity(clustrel.rels)
  assert len(clusters[0]) == 0
  vids, membership = make_membership_mat(clusters[1:])
  # K: number of non-empty clusters
  M, K = membership.shape

  num_models = len(common.Models._all)
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

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--weight-trees-by', choices=('llh', 'uniform'), default='uniform')
  parser.add_argument('--clustrel-mutrel')
  parser.add_argument('--trees-mutrel')
  parser.add_argument('--pure-mutrel')
  parser.add_argument('--phi')
  parser.add_argument('pairtree_results_fn')
  args = parser.parse_args()

  results = resultserializer.load(args.pairtree_results_fn)
  clusters = [[]] + list(results['clusters'])
  garbage = list(results['garbage'])
  all_vids = set([V for C in results['clusters'] for V in C] + garbage)

  if args.pure_mutrel is not None:
    pure_mutrel = results['mutrel_posterior']
    assert set(pure_mutrel.vids) == all_vids
    evalutil.save_sorted_mutrel(pure_mutrel, args.pure_mutrel)

  if args.trees_mutrel is not None:
    tree_mutrel = calc_mutrel_from_trees(results['adjm'], results['llh'], clusters, args.weight_trees_by)
    tree_mutrel = mutrel.add_garbage(tree_mutrel, garbage)
    assert set(tree_mutrel.vids) == all_vids
    evalutil.save_sorted_mutrel(tree_mutrel, args.trees_mutrel)

  if args.clustrel_mutrel is not None:
    clustrel_mutrel = calc_mutrel_from_clustrel(results['clustrel_posterior'], clusters)
    clustrel_mutrel = mutrel.add_garbage(clustrel_mutrel, garbage)
    assert set(clustrel_mutrel.vids) == all_vids
    evalutil.save_sorted_mutrel(clustrel_mutrel, args.clustrel_mutrel)

  if args.phi is not None:
    np.savez_compressed(args.phi, phi=results['phi'])

if __name__ == '__main__':
  main()
