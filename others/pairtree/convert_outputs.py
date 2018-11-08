import argparse
import sys
import os
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import common
import inputparser

def softmax(V):
  V = np.copy(V) - np.max(V)
  return np.exp(V) / np.sum(np.exp(V))

def calc_tree_mutrel(adjms, llhs, clusters):
  weights = softmax(llhs)
  soft_mutrel = None
  for adjm, weight in zip(adjms, weights):
    mutrel = common.make_mutrel_tensor_from_cluster_adj(adjm, clusters)
    if soft_mutrel is None:
      soft_mutrel = np.zeros(mutrel.shape)
    soft_mutrel += weight * mutrel
  return soft_mutrel

def make_membership_mat(clusters):
  children = [vid for C in clusters for vid in C]
  N = len(children)
  K = len(clusters)
  max_child = max(children)
  children = set(children)
  assert N == len(children)
  assert children == set(range(N))

  # membership[i,j] = 1 iff mutation `i` is in cluster `j`
  membership = np.zeros((N, K))
  for cidx, C in enumerate(clusters):
    membership[C,cidx] = 1
  return membership

def calc_pairs_mutrel(posterior, clusters):
  assert len(clusters[0]) == 0
  membership = make_membership_mat(clusters[1:])
  # K: number of non-empty clusters
  M, K = membership.shape

  num_models = len(common.Models._all)
  mutrel = np.zeros((M, M, num_models))
  assert posterior.shape == (K, K, num_models)

  for modelidx in range(num_models):
    mut_vs_cluster = np.dot(membership, posterior[:,:,modelidx]) # MxK
    mutrel[:,:,modelidx] = np.dot(mut_vs_cluster, mut_vs_cluster.T) # KxK
  assert np.all(0 <= posterior) and np.all(posterior <= 1)
  assert np.all(0 <= mutrel) and np.all(mutrel <= 1)
  return mutrel

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('pairtree_results_fn')
  parser.add_argument('pairs_mutrel_fn')
  parser.add_argument('trees_mutrel_fn')
  args = parser.parse_args()

  params = inputparser.load_params(args.pairtree_params_fn)
  clusters = params['clusters']
  results = np.load(args.pairtree_results_fn)

  tree_mutrel = calc_tree_mutrel(results['adjm'], results['llh'], clusters)
  np.savez_compressed(args.trees_mutrel_fn, soft_mutrel=tree_mutrel)

  pairs_mutrel = calc_pairs_mutrel(results['posterior'], clusters)
  np.savez_compressed(args.pairs_mutrel_fn, soft_mutrel=pairs_mutrel)

if __name__ == '__main__':
  main()
