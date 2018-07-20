import numpy as np
import scipy.stats
import common
import concurrent.futures
from tqdm import tqdm

# Matrices: M mutations, K clusters, S samples
#   adj: KxK, adjacency matrix -- adj[a,b]=1 iff a is parent of b (with 1 on diagonal)
#   Z: KxK, Z[a,b]=1 iff cluster a is ancestral to cluster b (with 1 on diagonal)
#   A: MxK, A[m,k]=1 iff mut m is in cluster k
#   phi: Kx1, per-cluster phis
#   mut_phi: Mx1, per-mutation phis

def calc_mut_p(A, Z, psi):
  eta = softmax(psi) # Kx1
  phi = np.dot(Z, eta) # Kx1
  mut_phi = np.dot(A, phi) # Mx1
  mut_p = 0.5 * mut_phi

  # Avoid numerical issues.
  delta = 1e-30
  mut_p[np.isclose(mut_p, 0)]   += delta
  mut_p[np.isclose(mut_p, 0.5)] -= delta

  return mut_p

def calc_llh(var_reads, ref_reads, A, Z, psi):
  total_reads = var_reads + ref_reads
  mut_p = calc_mut_p(A, Z, psi)
  mut_probs = scipy.stats.binom.logpmf(var_reads, total_reads, mut_p)
  return np.sum(mut_probs)

def fit_phis(adj, superclusters, supervars, iterations, parallel):
  supervar_ids = sorted(supervars.keys(), key = lambda C: int(C[1:]))
  ref_reads = np.array([supervars[cid]['ref_reads'] for cid in supervar_ids])
  var_reads = np.array([supervars[cid]['var_reads'] for cid in supervar_ids])
  assert len(supervars) == len(adj) - 1
  # Supervar `i` is in cluster `i`. Cluster 0 is empty.
  A = np.insert(np.eye(len(supervars)), 0, 0, axis=1)
  return fit_all_phis(adj, A, ref_reads, var_reads, iterations, parallel)

def fit_phi_S(eta_S, var_reads_S, ref_reads_S, A, Z, iterations):
  eta_S = np.maximum(1e-10, eta_S)
  psi_S = np.log(eta_S)
  psi_S = grad_desc(var_reads_S, ref_reads_S, A, Z, psi_S, iterations)
  eta_S = softmax(psi_S)
  return eta_S

def fit_all_phis(adj, A, ref_reads, var_reads, iterations, parallel):
  Z = common.make_ancestral_from_adj(adj)
  M, K = A.shape
  _, S = ref_reads.shape

  phi_implied = 2*(var_reads / (ref_reads + var_reads))
  # Make first row 1 to account for normal root.
  phi_implied = np.insert(phi_implied, 0, 1, axis=0)
  # Since phi = Z.eta, we have eta = (Z^-1)phi.
  eta = np.dot(np.linalg.inv(Z), phi_implied).T
  assert eta.shape == (S, K)

  modified_samples = []
  futures = []
  with concurrent.futures.ProcessPoolExecutor(max_workers=parallel) as ex:
    for s in range(S):
      # If no negative elements are in eta, we have the analytic solution.
      # Otherwise, set the negative elements to zero to provide an initial
      # starting point, then run the optimizer.
      if np.any(eta[s] < 0):
        modified_samples.append(s)
        futures.append(ex.submit(fit_phi_S, eta[s], var_reads[:,s], ref_reads[:,s], A, Z, iterations))
    for F in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc='Fitting phis', unit=' samples', dynamic_ncols=True):
      pass
  for s, F in zip(modified_samples, futures):
    eta[s] = F.result()

  phi = np.dot(Z, eta.T).T
  return (phi, eta)

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

def grad_desc(var_reads, ref_reads, A, Z, psi, iterations):
  learn_rate = 1e-4

  K = len(Z)
  M = len(var_reads)

  last_llh = -float('inf')
  last_psi = None

  for I in range(iterations):
    grad = calc_grad(var_reads, ref_reads, A, Z, psi)
    #grad_num = calc_grad_numerical(var_reads, ref_reads, A, Z, psi)
    #print(grad - grad_num)

    psi += learn_rate * grad
    llh = calc_llh(var_reads, ref_reads, A, Z, psi)

    if llh > last_llh:
      #print('%s\tkeep\t%.2e\t%.2e' % (I, learn_rate, llh))
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

def extract_mut_info(clusters, variants):
  # Renumber variants so their indices are contiguous. They may not be
  # contiguous, e.g., when singleton clusters are removed.
  used_vars = set([V for C in clusters for V in C])
  var_map = {old: new for (new, old) in enumerate(sorted(used_vars))}

  S = len(list(variants.values())[0]['total_reads']) # Number of samples
  M = len(used_vars)
  K = len(clusters)

  A = np.zeros((M, K))
  ref_reads = np.zeros((M, S))
  var_reads = np.zeros((M, S))

  for cidx, C in enumerate(clusters):
    for V in C:
      # Copy variant so we don't modify original dict.
      variant = variants['s%s' % V]
      varidx = var_map[V]
      A[varidx,cidx] = 1
      ref_reads[varidx] = variant['ref_reads']
      var_reads[varidx] = variant['var_reads']

  return (A, ref_reads, var_reads)
