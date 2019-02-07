import argparse
import pickle
import numpy as np
from numba import njit

@njit
def count_trees(tau, phi, order, traversal):
  assert traversal == 'dfs' or traversal == 'bfs'
  K = len(tau)
  expected_colsum = np.ones(K)
  expected_colsum[0] = 0

  first_partial = np.copy(tau)
  np.fill_diagonal(first_partial, 0)
  first_delta = np.copy(phi)
  partial_trees = [(1, first_partial, first_delta)]
  completed_trees = 0

  while len(partial_trees) > 0:
    if traversal == 'dfs':
      to_resolve, P, delta = partial_trees.pop()
    else:
      to_resolve, P, delta = partial_trees.pop(0)
      #to_resolve, P, delta = partial_trees[0]
      #partial_trees = partial_trees[1:]
    if to_resolve == K:
      assert np.all(expected_colsum == np.sum(P, axis=0))
      assert np.all(0 <= delta) and np.all(delta <= 1)
      np.fill_diagonal(P, 1)
      completed_trees += 1
      continue

    R = order[to_resolve]
    parents = np.nonzero(P[:,R])[0]
    for parent in parents:
      P_prime = np.copy(P)
      P_prime[:,R] = 0
      P_prime[parent,R] = 1
      if np.any(delta[parent] - phi[R] < 0):
        continue
      delta_prime = np.copy(delta)
      delta_prime[parent] -= phi[R]
      partial_trees.append((to_resolve + 1, P_prime, delta_prime))
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
  tau = np.eye(K)

  for I in range(K):
    for J in range(I + 1, K):
      I_prime = order[I]
      J_prime = order[J]
      assert not np.all(phi[I_prime] == phi[J_prime])
      if np.all(phi[I_prime] >= phi[J_prime]):
        tau[I_prime,J_prime] = 1

  return tau

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sim_data_fn')
  args = parser.parse_args()

  with open(args.sim_data_fn, 'rb') as dataf:
    simdata = pickle.load(dataf)

  phi, true_tree = simdata['phi'], simdata['adjm']
  order = make_order(phi)
  tau = make_tau(phi, order)
  num_trees = count_trees(tau, phi, order, 'dfs')
  print(args.sim_data_fn, num_trees)

main()
