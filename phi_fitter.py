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

  def _find_desc(I, vec):
    # Base case: if I have no children, my ancestor vec is just myself.
    if np.sum(vec) == 0:
      return vec
    else:
      children = np.array([_find_desc(idx, adj[idx]) for (idx, val) in enumerate(vec) if idx != I and val == 1])
      self_and_child = vec + np.sum(children, axis=0)
      self_and_child[self_and_child > 1] = 1
      return self_and_child

  for k in range(K):
    # If we know `adj` is topologically sorted, we can reduce the complexity of
    # this -- we would start at leaves and work our way upward, eliminating
    # need for recursive DFS. But since we don't expect `K` to be large, we can
    # write more general version that works for non-toposorted trees.
    Z[k] = _find_desc(k, adj[k])

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

def fit_all_phis(adj, A, ref_reads, var_reads):
  Z = calc_Z(adj)
  M, K = A.shape
  _, S = ref_reads.shape
  psi = np.zeros((S, K))
  phi = np.zeros((S, K))

  for s in range(S):
    psi[s] = grad_desc(var_reads[:,s], ref_reads[:,s], A, Z)
    phi[s] = np.dot(Z, softmax(psi[s]))

  return phi

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

  grad = 0.5 * np.sum(grad_elem, axis=0) / M
  return grad

def calc_grad_numerical(var_reads, ref_reads, A, Z, psi):
  _, K = A.shape
  grad = np.zeros(K)
  for k in range(K):
    delta = 1e-20
    P = np.copy(psi)
    P[k] += delta
    grad[k] = (calc_llh(var_reads, ref_reads, A, Z, P) - calc_llh(var_reads, ref_reads, A, Z, psi)) / delta
  return grad

def grad_desc(var_reads, ref_reads, A, Z):
  learn_rate = 1e-4
  iterations = 1000

  K = len(Z)
  M = len(var_reads)

  last_llh = -float('inf')
  last_psi = None
  psi = np.random.normal(size=K)

  for I in range(iterations):
    grad = calc_grad(var_reads, ref_reads, A, Z, psi)
    #grad_num = calc_grad_numerical(var_reads, ref_reads, A, Z, psi)
    #print(grad - grad_num)

    psi += learn_rate * grad
    llh = calc_llh(var_reads, ref_reads, A, Z, psi)

    if llh > last_llh:
      print('%s\tkeep\t%.2e\t%.2e' % (I, learn_rate, llh))
      last_llh = llh
      last_psi = psi
      learn_rate *= 1.1
    else:
      #print('%s\trevert\t%.2e\t%.2e' % (I, learn_rate, llh))
      llh = last_llh
      psi = last_psi
      learn_rate *= 0.5

  return psi

def softmax(X):
  b = np.max(X)
  logS = X - (b + np.log(np.sum(np.exp(X - b))))
  return np.exp(logS)
