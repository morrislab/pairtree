import argparse
import os
import numpy as np

import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import common
import resultserializer

def softmax(V):
  V = np.copy(V) - np.max(V)
  return np.exp(V) / np.sum(np.exp(V))

def calc_mutrel_from_trees(adjms, llhs, clusters):
  #weights = softmax(llhs)
  weights = np.eye(len(adjms)) / len(adjms)

  soft_mutrel = None
  for adjm, weight in zip(adjms, weights):
    mutrel = common.make_mutrel_tensor_from_cluster_adj(adjm, clusters)
    if soft_mutrel is None:
      soft_mutrel = np.zeros(mutrel.rels.shape)
    soft_mutrel += weight * mutrel.rels
  return soft_mutrel

def make_membership_mat(clusters):
  children = [int(vid[1:]) for C in clusters for vid in C]
  N = len(children)
  K = len(clusters)
  max_child = max(children)
  children = set(children)
  assert N == len(children)
  assert children == set(range(N))

  # membership[i,j] = 1 iff mutation `i` is in cluster `j`
  membership = np.zeros((N, K))
  for cidx, C in enumerate(clusters):
    members = [int(vid[1:]) for vid in C]
    membership[members,cidx] = 1
  return membership

def calc_mutrel_from_clustrel(clustrel, clusters):
  assert np.all(0 <= clustrel.rels) and np.all(clustrel.rels <= 1)
  assert len(clusters[0]) == 0
  membership = make_membership_mat(clusters[1:])
  # K: number of non-empty clusters
  M, K = membership.shape

  num_models = len(common.Models._all)
  mutrel = np.zeros((M, M, num_models))
  assert clustrel.rels.shape == (K, K, num_models)

  for modelidx in range(num_models):
    mut_vs_cluster = np.dot(membership, clustrel.rels[:,:,modelidx]) # MxK
    mutrel[:,:,modelidx] = np.dot(mut_vs_cluster, membership.T)
  assert np.all(0 <= mutrel) and np.all(mutrel <= 1)
  return mutrel

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--clustrel-mutrel')
  parser.add_argument('--trees-mutrel')
  parser.add_argument('--pure-mutrel')
  parser.add_argument('--phi')
  parser.add_argument('pairtree_results_fn')
  args = parser.parse_args()

  results = resultserializer.load(args.pairtree_results_fn)
  clusters = [[]] + list(results['clusters'])

  if args.pure_mutrel is not None:
    np.savez_compressed(args.pure_mutrel, mutrel=results['mutrel_posterior'].rels)
  if args.trees_mutrel is not None:
    tree_mutrel = calc_mutrel_from_trees(results['adjm'], results['llh'], clusters)
    np.savez_compressed(args.trees_mutrel, mutrel=tree_mutrel)
  if args.clustrel_mutrel is not None:
    clustrel_mutrel = calc_mutrel_from_clustrel(results['clustrel_posterior'], clusters)
    np.savez_compressed(args.clustrel_mutrel, mutrel=tree_mutrel)
  if args.phi is not None:
    np.savez_compressed(args.phi, phi=results['phi'])

if __name__ == '__main__':
  main()
