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
  # Renumber variants so their indices are contiguous. They may not be
  # contiguous, e.g., when singleton clusters are removed.
  used_vars = set([V for C in clusters for V in C])
  var_map = {old: new for (new, old) in enumerate(sorted(used_vars))}

  S = len(list(variants.values())[0]['total_reads']) # Number of samples
  M = len(used_vars)
  K = len(cidxs)
  assert set(cidxs) == set(range(K))

  A = np.zeros((M, K))
  ref_reads = np.zeros((M, S))
  var_reads = np.zeros((M, S))

  for cidx, C in zip(cidxs, clusters):
    for V in C:
      # Copy variant so we don't modify original dict.
      variant = variants['s%s' % V]
      varidx = var_map[V]
      A[varidx,cidx] = 1
      ref_reads[varidx] = variant['ref_reads']
      var_reads[varidx] = variant['var_reads']

  return (A, ref_reads, var_reads)

def fit_all_phis(adj, A, ref_reads, var_reads):
  Z = calc_Z(adj)
  M, K = A.shape
  _, S = ref_reads.shape
  psi = np.zeros((S, K))

  for s in range(S):
    psi[s] = grad_desc(var_reads[:,s], ref_reads[:,s], A, Z)

def calc_grad(var_reads, ref_reads, A, Z, psi):
  M, K = A.shape
  AZ = np.dot(A, Z) # MxK

  mut_p = calc_mut_p(A, Z, psi)
  binom_grad = (var_reads / mut_p) - (ref_reads / (1 - mut_p)) # Mx1
  assert len(binom_grad) == M

  eta = softmax(psi)
  eta_rows = np.tile(eta, (K, 1))
  softmax_grad = eta_rows * (np.eye(K) - eta_rows.T) # KxK

  grad_elem = np.zeros((M, K))
  for m in range(M):
    active_b = softmax_grad * np.tile(AZ[m], (K, 1))
    grad_elem[m] = binom_grad[m] * np.sum(active_b, axis=1)

  grad = 0.5 * np.sum(grad_elem, axis=0)
  return grad

def grad_desc(var_reads, ref_reads, A, Z):
  learn_rate = 0.5
  iterations = 1000

  K = len(Z)
  M = len(var_reads)

  last_llh = -float('inf')
  last_psi = None
  psi = np.random.normal(size=K)

  for I in range(iterations):
    grad = calc_grad(var_reads, ref_reads, A, Z, psi)
    psi += learn_rate * grad
    llh = calc_llh(var_reads, ref_reads, A, Z, psi)

    if llh > last_llh:
      print('%s\tkeep\t%.2e\t%.2e' % (I, learn_rate, llh))
      last_llh = llh
      last_psi = psi
      learn_rate *= 1.1
    else:
      print('%s\trevert\t%.2e\t%.2e' % (I, learn_rate, llh))
      llh = last_llh
      psi = last_psi
      learn_rate *= 0.5

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
  A, ref_reads, var_reads = extract_mut_info(clusters, cidxs, variants)
  fit_all_phis(adj, A, ref_reads, var_reads)
