from common import Models
import numpy as np
import scipy.stats

# Matrices: M mutations, K clusters
#   adj: KxK, adjacency matrix
#   Z: KxK, Z[a,b]=1 iff cluster a is ancestral to cluster b
#   A: MxK, A[m,k]=1 iff mut m is in cluster k
#   phi: Kx1, per-cluster phis
#   mut_phi: Mx1, per-mutation phis

def calc_Z(adj):
  K = len(adj)
  Z = np.zeros((K,K))

  def _find_desc(vec):
    # Given parent->child vector (i.e., row from `adj`), compute
    # parent->child->grandchild vector. Applying this recursively transforms
    # parent->child vector to parent->descendant vector.

    # Base case: if I have no children, my ancestor vec is just myself.
    leaf_vec = np.eye(K)[K - 1]
    if np.array_equal(vec, leaf_vec):
      return vec
    else:
      vec = vec[np.newaxis,:] # K -> 1xK
      self_and_child = np.dot(vec, adj) # Adjacency vector for self and child
      self_and_child[self_and_child > 1] = 1
      return self_and_child

  for k in range(K):
    # If we know `adj` is topologically sorted, we can reduce the complexity of
    # this -- we would start at leaves and work our way upward, eliminating
    # need for recursive DFS. But since we don't expect `K` to be large, we can
    # write more general version that works for non-toposorted trees.
    Z[k] = _find_desc(adj[k])

  return Z

def calc_mut_p(A, Z, psi):
  eta = softmax(psi) # Kx1
  phi = np.dot(Z, eta) # Kx1
  mut_phi = np.dot(A, phi) # Mx1
  mut_p = 0.5 * mut_phi
  return mut_p

def calc_llh(var_reads, ref_reads, A, Z, psi):
  total_reads = var_reads + ref_reads
  mut_p = calc_mut_p(A, Z, psi)
  mut_probs = scipy.stats.binom.logpmf(var_reads, total_reads, mut_p)
  return np.sum(mut_probs)

def extract_mut_info(clusters, cidxs, variants):
  S = len(list(variants.values())[0]['total_reads']) # Number of samples
  M = len(variants)
  K = len(cidxs)
  assert set(cidxs) == set(range(K))

  A = np.zeros((M, K))
  ref_reads = np.zeros((M, S))
  var_reads = np.zeros((M, S))

  for cidx, C in zip(cidxs, clusters):
    for V in C:
      # Copy variant so we don't modify original dict.
      variant = variants['s%s' % V]
      A[V,cidx] = 1
      ref_reads[V] = variant['ref_reads']
      var_reads[V] = variant['var_reads']

  return (A, ref_reads, var_reads)

def fit_all_phis(adj, clusters, cidxs, variants):
  Z = calc_Z(adj)
  A, ref_reads, var_reads = extract_mut_info(clusters, cidxs, variants)
  M, K = A.shape
  _, S = ref_reads.shape
  psi = np.zeros((S, K))

  for s in range(S):
    psi[s] = np.random.normal(size=K)
    llh = calc_llh(var_reads[:,s], ref_reads[:,s], A, Z, psi[s])
    print('sample', s, llh)

def _calc_grad(var_reads, ref_reads, A, Z, psi):
  mut_p = calc_mut_p(A, Z, psi)
  grad = (var_reads / mut_p) - (ref_reads / (1 - mut_p))

  eta = softmax(psi)


def grad_desc(var_reads, ref_reads, Z):
  learn_rate = 0.5
  iterations = 1000

  K = len(Z)
  M = len(var_reads)

  last_llh = -float('inf')
  last_psi = None

  for I in range(iterations):
    grad = 



def softmax(X):
  b = np.max(X)
  logS = X - (b + np.log(np.sum(np.exp(X - b))))
  return np.exp(logS)

def make_adj(relations):
  N = len(relations)
  assert relations.shape == (N, N)
  assert relations[1,0] == Models.B_A

  adj = np.eye(N)
  adj[0,1] = 1

  for I in range(2, N):
    I_placed = False

    for J in reversed(range(I)):
      if relations[I,J] == Models.B_A:
        if I_placed is False:
          adj[J,I] = 1
          I_placed = True
      elif relations[I,J] == Models.diff_branches:
        pass
      else:
        raise Exception('Unexpected relation for (%s,%s): %s' % (I, J, Models._all[relations[I,J]]))

  return adj

def build_tree(relations, clusters, cidxs, variants):
  adj = make_adj(relations)
  fit_all_phis(adj, clusters, cidxs, variants)
