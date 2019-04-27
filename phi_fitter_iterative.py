import numpy as np
import scipy.stats
import common
import binom
from progressbar import progressbar
from numba import njit

# Matrices: M mutations, K clusters, S samples
#   adj: KxK, adjacency matrix -- adj[a,b]=1 iff a is parent of b (with 1 on diagonal)
#   Z: KxK, Z[a,b]=1 iff cluster a is ancestral to cluster b (with 1 on diagonal)
#   A: MxK, A[m,k]=1 iff mut m is in cluster k
#   phi: Kx1, per-cluster phis
#   mut_phi: Mx1, per-mutation phis

@njit
def calc_mut_p(A, Z, psi):
  eta = softmax(psi) # Kx1
  phi = np.dot(Z, eta) # Kx1
  mut_phi = np.dot(A, phi) # Mx1
  # TODO: should this use omega_v?
  mut_p = 0.5 * mut_phi

  # Avoid numerical issues.
  delta = 1e-30
  mut_p[binom.isclose(mut_p, 0)]   += delta
  mut_p[binom.isclose(mut_p, 0.5)] -= delta

  return mut_p

@njit
def calc_llh(var_reads, ref_reads, A, Z, psi):
  total_reads = var_reads + ref_reads
  mut_p = calc_mut_p(A, Z, psi)
  mut_probs = binom.logpmf(var_reads, total_reads, mut_p)
  return np.sum(mut_probs)

def fit_phis(adj, superclusters, supervars, method, iterations, parallel):
  svids = common.extract_vids(supervars)
  ref_reads = np.array([supervars[svid]['ref_reads'] for svid in svids])
  var_reads = np.array([supervars[svid]['var_reads'] for svid in svids])

  # Supervar `i` is in cluster `i`. Cluster 0 is empty.
  assert len(supervars) == len(adj) - 1
  A = np.insert(np.eye(len(supervars)), 0, 0, axis=1)

  return fit_all_phis(adj, A, ref_reads, var_reads, method, iterations, parallel)

@njit
def fit_phi_S(eta_S, var_reads_S, ref_reads_S, A, Z, method, iterations):
  eta_S = np.maximum(common._EPSILON, eta_S)
  psi_S = np.log(eta_S)

  if method == 'graddesc':
    psi_S = grad_desc(var_reads_S, ref_reads_S, A, Z, psi_S, iterations)
  elif method == 'rprop':
    psi_S = rprop(var_reads_S, ref_reads_S, A, Z, psi_S, iterations)
  else:
    raise Exception('Unknown psi fitter')

  eta_S = softmax(psi_S)
  return eta_S

def fit_all_phis(adj, A, ref_reads, var_reads, method, iterations, parallel):
  # TODO: I can probably rip out all the parallelism machinery, since I don't
  # use this any longer.
  Z = common.make_ancestral_from_adj(adj)
  M, K = A.shape
  _, S = ref_reads.shape

  phi_implied = 2*(var_reads / (ref_reads + var_reads))
  # Make first row 1 to account for normal root.
  phi_implied = np.insert(phi_implied, 0, 1, axis=0)
  # Since phi = Z.eta, we have eta = (Z^-1)phi.
  eta = np.dot(np.linalg.inv(Z), phi_implied).T
  assert eta.shape == (S, K)

  if parallel > 0:
    modified_samples = []
    futures = []
    import concurrent.futures
    with concurrent.futures.ProcessPoolExecutor(max_workers=parallel) as ex:
      for s in range(S):
        # If no negative elements are in eta, we have the analytic solution.
        # Otherwise, set the negative elements to zero to provide an initial
        # starting point, then run the optimizer.
        if np.any(eta[s] < 0):
          modified_samples.append(s)
          futures.append(ex.submit(fit_phi_S, eta[s], var_reads[:,s], ref_reads[:,s], A, Z, method, iterations))
      with progressbar(total=len(futures), desc='Fitting phis', unit='sample', dynamic_ncols=True) as pbar:
        for F in concurrent.futures.as_completed(futures):
          pbar.update()
    for s, F in zip(modified_samples, futures):
      eta[s] = F.result()
  else:
    for s in range(S):
      if np.any(eta[s] < 0):
        eta[s] = fit_phi_S(eta[s], var_reads[:,s], ref_reads[:,s], A, Z, method, iterations)

  phi = np.dot(Z, eta.T)
  return (phi, eta.T)

@njit
def _tile_rows(vec, repeats):
  assert vec.ndim == 1
  mat = np.empty((repeats, len(vec)))
  mat[:] = vec
  return mat

@njit
def calc_grad(var_reads, ref_reads, A, Z, psi):
  M, K = A.shape
  AZ = np.dot(A, Z) # MxK

  mut_p = calc_mut_p(A, Z, psi)
  binom_grad = (var_reads / mut_p) - (ref_reads / (1 - mut_p)) # Mx1
  assert len(binom_grad) == M

  eta = softmax(psi)
  eta_rows = _tile_rows(eta, K)
  softmax_grad = eta_rows * (np.eye(K) - eta_rows.T) # KxK

  grad_elem = np.zeros((M, K))
  for m in range(M):
    active_b = softmax_grad * _tile_rows(AZ[m], K)
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

@njit
def grad_desc(var_reads, ref_reads, A, Z, psi, iterations, convergence_threshold=1e-30):
  learn_rate = 1e-4
  last_llh = -np.inf
  last_psi = psi

  for I in range(iterations):
    grad = calc_grad(var_reads, ref_reads, A, Z, psi)
    #grad_num = calc_grad_numerical(var_reads, ref_reads, A, Z, psi)
    #print(grad - grad_num)

    step = learn_rate * grad
    if np.max(np.abs(step)) < convergence_threshold:
      break
    psi += step
    llh = calc_llh(var_reads, ref_reads, A, Z, psi)

    if llh > last_llh:
      #print(I, 'keep', '%.2e' % learn_rate, '%.2e' % llh, sep='\t')
      last_llh = llh
      last_psi = psi
      learn_rate *= 1.1
    else:
      #print(I, 'revert', '%.2e' % learn_rate, '%.2e' % llh, sep='\t')
      llh = last_llh
      psi = last_psi
      learn_rate *= 0.5

  return psi

@njit
def rprop(var_reads, ref_reads, A, Z, psi, iterations, convergence_threshold=1e-4):
  # For explanation of `rprop`, see
  # http://www.cs.toronto.edu/~tijmen/csc321/slides/lecture_slides_lec6.pdf.
  step = 1e-4 * np.ones(psi.shape)
  last_sign = np.ones(psi.shape)

  for T in range(iterations):
    grad = calc_grad(var_reads, ref_reads, A, Z, psi)
    if np.max(np.abs(grad)) < convergence_threshold:
      break
    sign = np.ones(psi.shape)
    sign[grad < 0] = -1
    sign_agree = sign == last_sign
    sign_disagree = np.logical_not(sign_agree)
    last_sign = sign

    delta = 1.2*sign_agree + 0.5*sign_disagree
    delta = np.minimum(50, delta)
    step *= delta
    psi += sign * step

  return psi

@njit
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
