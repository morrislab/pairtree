import argparse
import pickle
import numpy as np
from numba import njit

import evalutil

import sys
import os

@njit
def _find_parents(adj):
  K = len(adj)
  assert adj.shape == (K, K)
  parents = np.full(K-1, np.nan)
  for j in range(1, K):
    for i in range(K):
      if i == j:
        continue
      if adj[i,j]:
        parents[j-1] = i
        break
  assert not np.any(np.isnan(parents))
  return parents.astype(np.uint8)

@njit
def enum_trees(tau, phi, order, traversal):
  assert traversal == 'dfs' or traversal == 'bfs'
  epsilon = 1e-10
  K = len(tau)
  expected_colsum = np.ones(K)
  expected_colsum[0] = 0

  tau_copy = np.copy(tau)
  first_partial = tau_copy.astype(np.int8)
  np.fill_diagonal(first_partial, 0)
  first_childsum = np.zeros(phi.shape)
  partial_trees = [(1, first_partial, first_childsum)]
  completed_trees = []

  while len(partial_trees) > 0:
    if traversal == 'dfs':
      to_resolve, P, childsum = partial_trees.pop()
    else:
      to_resolve, P, childsum = partial_trees.pop(0)
    if to_resolve == K:
      diff = phi - childsum
      assert np.all(expected_colsum == np.sum(P, axis=0))
      assert np.all(0 <= childsum + epsilon) and np.all(childsum <= 1 + epsilon)
      assert np.all(childsum <= phi + epsilon)
      np.fill_diagonal(P, 1)
      struct = _find_parents(P)
      completed_trees.append(struct)
      continue

    R = order[to_resolve]
    parents = np.nonzero(P[:,R])[0]
    for parent in parents:
      childsum_prime = np.copy(childsum[parent])
      childsum_prime += phi[R]
      if np.any(childsum_prime > phi[parent] + epsilon):
        continue
      new_childsum = np.copy(childsum)
      new_childsum[parent] = childsum_prime
      P_prime = np.copy(P)
      P_prime[:,R] = 0
      P_prime[parent,R] = 1
      partial_trees.append((to_resolve + 1, P_prime, new_childsum))

  return completed_trees

@njit
def make_order(phi):
  phisum = np.sum(phi, axis=1)
  order = np.argsort(-phisum)
  assert order[0] == 0
  return order

@njit
def make_tau(phi, order):
  K, S = phi.shape
  tau = np.eye(K, dtype=np.int8)

  for I in range(K):
    for J in range(I + 1, K):
      I_prime = order[I]
      J_prime = order[J]
      assert not np.all(phi[I_prime] == phi[J_prime])
      if np.all(phi[I_prime] >= phi[J_prime]):
        tau[I_prime,J_prime] = 1

  return tau

def ensure_truth_found(truth, trees):
  for T in trees:
    if np.all(T == truth):
      return
  raise Exception('No truth found')

def write_truth(structs, phi, clusters, garbage, results_fn):
  N = len(structs)
  llhs = np.zeros(N)
  probs = np.ones(N) / N
  phis = np.array([phi for _ in range(N)])
  counts = np.ones(N)

  np.savez_compressed(
    results_fn,
    struct = structs,
    count = counts,
    phi = phis,
    llh = llhs,
    prob = probs,
    clusters = clusters,
    garbage = garbage,
  )

def check_true_delta(struct, phi):
  import accupy
  delta = np.zeros(phi.shape)
  queue = [0]
  while len(queue) > 0:
    parent = queue.pop()
    children = np.flatnonzero(struct == parent) + 1
    delta[parent] = accupy.fsum(np.vstack((phi[parent], -phi[children])))
    queue += children.tolist()
  # Because of numerical errors some elements may be less than zero.
  # This may occur for any one of three reasons:
  #   1. Numerical errors in sampling the `eta` values, such that a sample's etas sum to >1
  #   2. Numerical error in summing phis as "sum of eta of all descendants"
  #   3. Numerical error in the computation of the delta within this function
  # We use `accupy.fsum()` to correct reason 3 to the greatest extent possible,
  # but, of course, this does nothing to help us if the errors result from
  # reason 1 or reason 2.
  assert np.all(np.isclose(0, delta[delta < 0]))
  assert np.all(delta <= 1)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--check-delta', action='store_true')
  parser.add_argument('sim_data_fn')
  parser.add_argument('results_fn')
  args = parser.parse_args()

  with open(args.sim_data_fn, 'rb') as dataf:
    simdata = pickle.load(dataf)
  if args.check_delta:
    check_true_delta(simdata['structure'], simdata['phi'])

  phi = simdata['phi']
  order = make_order(phi)
  tau = make_tau(phi, order)
  structs = enum_trees(tau, phi, order, 'dfs')
  ensure_truth_found(simdata['structure'], structs)
  write_truth(
    structs,
    simdata['phi'],
    simdata['clusters'],
    simdata['vids_garbage'],
    args.results_fn,
  )

if __name__ == '__main__':
  main()
