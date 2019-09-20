import argparse
import pickle
import numpy as np
from numba import njit

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import mutrel
import evalutil

@njit
def enum_trees(tau, phi, order, traversal):
  assert traversal == 'dfs' or traversal == 'bfs'
  K = len(tau)
  expected_colsum = np.ones(K)
  expected_colsum[0] = 0

  first_partial = np.copy(tau)
  np.fill_diagonal(first_partial, 0)
  first_delta = np.copy(phi)
  partial_trees = [(1, first_partial, first_delta)]
  completed_trees = []

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
      completed_trees.append(P)
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

def make_mutrel(trees, clusters):
  soft_mutrel = None

  for adjm in trees:
    mrel = mutrel.make_mutrel_tensor_from_cluster_adj(adjm, clusters)
    if soft_mutrel is None:
      soft_mutrel = mrel
    else:
      assert np.all(soft_mutrel.vids == mrel.vids)
      soft_mutrel = soft_mutrel._replace(rels = soft_mutrel.rels + mrel.rels)

  soft_mutrel = soft_mutrel._replace(rels = evalutil._fix_rounding_errors(soft_mutrel.rels / len(trees)))
  return soft_mutrel

def write_mutrel(trees, clusters, simdata, mrelfn):
  mrel = make_mutrel(trees, clusters)
  mrel = evalutil.add_garbage(mrel, simdata['vids_garbage'])

  assert set(mrel.vids) == set(simdata['vids_good'] + simdata['vids_garbage'])
  mrel = mutrel.sort_mutrel_by_vids(mrel)
  mutrel.check_posterior_sanity(mrel.rels)
  np.savez_compressed(
    mrelfn,
    rels = mrel.rels,
    vids = mrel.vids,
    num_trees = len(trees),
  )

def ensure_truth_found(truth, trees):
  for T in trees:
    if np.all(T == truth):
      return
  raise Exception('No truth found')

def write_truth(adjms, phi, truth_fn):
  N = len(adjms)
  llhs = np.zeros(N)
  phis = np.array([phi for _ in range(N)])
  np.savez_compressed(truth_fn, adjm=adjms, phi=phis, llh=llhs)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--enumerate-trees', action='store_true')
  parser.add_argument('--mutrel', dest='mutrel_fn')
  parser.add_argument('--truth', dest='truth_fn')
  parser.add_argument('sim_data_fn')
  args = parser.parse_args()

  with open(args.sim_data_fn, 'rb') as dataf:
    simdata = pickle.load(dataf)
  clusters = [[]] + simdata['clusters']

  if args.enumerate_trees:
    phi = simdata['phi']
    order = make_order(phi)
    tau = make_tau(phi, order)
    trees = enum_trees(tau, phi, order, 'dfs')
    ensure_truth_found(simdata['adjm'], trees)
  else:
    raise Exception('I have decided that making mutrels without enumerating trees is stupid')
    trees = [simdata['adjm']]

  if args.mutrel_fn is not None:
    write_mutrel(trees, clusters, simdata, args.mutrel_fn)
  if args.truth_fn is not None:
    write_truth(trees, simdata['phi'], args.truth_fn)

if __name__ == '__main__':
  main()
