import argparse
import pickle
import numpy as np
from numba import njit

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import resultserializer

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
def enum_trees(tau, phi, order, traversal, store_trees=True, print_progress=False):
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
  # Explicitly make this a Numpy int. If I left this as a native Python int,
  # Numba would automatically have converted it to a Numpy type to allow a
  # fixed-length encoding using machine-native datatypes. If integer overflow
  # occurs, it will automatically convert to a float:
  # https://numba.pydata.org/numba-doc/dev/proposals/integer-typing.html. To
  # avoid converting to a float (and losing precision), we make this an
  # explicit int, then check for overflow in the loop.
  num_trees = np.uint64(0)
  max_count = np.iinfo(np.uint64).max

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
      # Ensure we won't overflow.
      assert num_trees < max_count
      # Adding a Python int to `num_trees`, which is a Numpy uint64, results in
      # `num_trees` becoming a float. This is because the Numpy casting rules
      # say that a Python int can be negative, since it's inherently signed, so
      # the result must be a sufficiently general type to accommodate this. To
      # avoid this behaviour, explicitly add a Numpy unsigned integer. (Adding
      # a Numpy signed integer would result in a float again for the same
      # reason.)
      num_trees += np.uint8(1)
      # Print every hundred million trees.
      if print_progress and num_trees % int(1e8) == 0:
        print('partial', num_trees)
      if store_trees:
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

  return (num_trees, completed_trees)

@njit
def make_order(phi):
  phisum = np.sum(phi, axis=1)
  order = np.argsort(-phisum)
  assert order[0] == 0
  return order

# TODO: propagate "definite parents" constraints through the tau graph to
# reduce the number of edges. This is an-front optimization that would improve
# performance.
#
# As the algorithm currently stands, if a high-phi node near the top of the
# tree has only one possible parent, there will still be a bunch of spurious
# edges between that parent and low-phi nodes near the bottom of the tree.
# These will be elimianted through enumeration, but result in additional
# iterations of the loop. I think the performance difference will be
# inconsequential, though -- for each of those low-phi nodes, there will be at
# most `K` invalid parents for it, so I think that I can upper-bound the number
# of "bad" trees in that case at `K^2`, which is small enough to not matter.
@njit
def make_tau(phi, order):
  K, S = phi.shape
  tau = np.eye(K, dtype=np.int8)

  for I in range(K):
    for J in range(I + 1, K):
      I_prime = order[I]
      J_prime = order[J]
      # Because we can now have `eta = 0` for some nodes, phis will sometimes
      # be identical across nodes. (This is particularly common with only one
      # sample, since it only "gets one chance" to have a different phi.) This
      # implies that there's no longer a meaningful ordering between those two
      # nodes, which in turn means that we may not receover the true tree
      # through enumeration. We may need to alter the simulations to have a
      # tiny (non-zero) minimum `eta` so that this situation doesn't occur.
      #assert not np.all(phi[I_prime] == phi[J_prime])
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

  results = resultserializer.Results(results_fn)
  results.add('struct', structs)
  results.add('count', counts)
  results.add('phi', phis)
  results.add('llh', llhs)
  results.add('prob', probs)
  results.add('clusters', clusters)
  results.add('garbage', garbage)
  results.save()

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
  parser.add_argument('--only-count', action='store_true')
  parser.add_argument('--only-write-truth', action='store_true', help='Write only the single true tree structure')
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

  if args.only_count:
    num_trees, _ = enum_trees(tau, phi, order, 'dfs', store_trees=False)
    print(num_trees)
    return
  elif args.only_write_truth:
    num_trees = 1
    structs = np.array(simdata['structure'])
  else:
    num_trees, structs = enum_trees(tau, phi, order, 'dfs')
    structs = np.array(structs)

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
