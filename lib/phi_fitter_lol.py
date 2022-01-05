import numpy as np
from numba import njit, jit, prange
import common
import binom
import util

_EPS = 1e-10

# Matrices: M mutations, K clusters, S samples
#   adj: KxK, adjacency matrix -- adj[a,b]=1 iff a is parent of b (with 1 on diagonal)
#   Z: KxK, Z[a,b]=1 iff cluster a is ancestral to cluster b (with 1 on diagonal)
#   A: MxK, A[m,k]=1 iff mut m is in cluster k
#   phi: Kx1, per-cluster phis
#   mut_phi: Mx1, per-mutation phis

def fit_etas(adj, superclusters, supervars, method, iterations, parallel, eta_init='mle'):
  # Supervar `i` is in cluster `i`. Cluster 0 is empty.
  assert len(supervars) == len(adj) - 1
  A = np.insert(np.eye(len(supervars)), 0, 0, axis=1)

  svids = common.extract_vids(supervars)
  arrs = {key: np.array([supervars[svid][key] for svid in svids]) \
    for key in ('ref_reads', 'var_reads', 'omega_v')}
  return _fit_etas(
    adj,
    A,
    arrs['ref_reads'],
    arrs['var_reads'],
    arrs['omega_v'],
    method,
    iterations,
    parallel,
    eta_init,
  )

def _fit_etas(adj, A, ref_reads, var_reads, omega, method, iterations, parallel, eta_init):
  M, K = A.shape
  _, S = ref_reads.shape
  assert ref_reads.shape == var_reads.shape == omega.shape == (M, S)

  Z = util.make_ancestral_from_adj(adj)
  # Numba only supports dot products on float arrays, not int arrays.
  A = A.astype(np.float64)
  Z = Z.astype(np.float64)
  
  # we want to operate on float arrays when computing phi_implied 
  var_reads = var_reads.astype(np.float64)
  ref_reads = ref_reads.astype(np.float64)

  if isinstance(eta_init, str) and eta_init == 'mle':

    # To prevent the true_divide error, we use numpy divide with the condition to only perform division if omega and var_reads are not zero.
    tmp_vaf = np.divide(var_reads, (ref_reads + var_reads), out=np.zeros_like(var_reads), where=var_reads!=0)
    phi_implied = np.divide(tmp_vaf, omega, out=np.zeros_like(tmp_vaf), where=omega!=0)
    phi_implied = np.minimum(1, phi_implied)
    phi_implied = np.maximum(0, phi_implied)
    # Make first row 1 to account for normal root.
    phi_implied = np.insert(phi_implied, 0, 1, axis=0)
    # Since phi = Z.eta, we have eta = (Z^-1)phi.
    eta = np.dot(np.linalg.inv(Z), phi_implied)
  elif isinstance(eta_init, str) and eta_init == 'dirichlet':
    eta = np.random.dirichlet(alpha = K*[1e0], size = S).T
  else:
    assert isinstance(eta_init, np.ndarray)
    eta = np.copy(eta_init)

  assert eta.shape == (K, S)
  assert np.allclose(1, np.sum(eta, axis=0))
  eta = np.maximum(_EPS, eta)
  eta /= np.sum(eta, axis=0)
  assert np.allclose(1, np.sum(eta, axis=0))

  # This could be parallelized. Maybe using numba.prange?
  for s in prange(S):
    eta[:,s] = _fit_eta_S(
      eta[:,s],
      var_reads[:,s],
      ref_reads[:,s],
      omega[:,s],
      A,
      Z,
      method,
      iterations,
    )

  return eta

@njit
def _tile_rows(vec, repeats):
  assert vec.ndim == 1
  mat = np.empty((repeats, len(vec)))
  mat[:] = vec
  return mat

@njit
def _calc_grad(var_reads, ref_reads, omega, A, Z, psi):
  # Prevent overflow in `exp`. Cap `psi` at 300, since
  # `np.log(np.finfo(np.float64).max) =~ 709`, and we take `exp(2*psi)` below.
  # Thus, we need `psi < 709/2`.
  psi = np.minimum(300, psi)

  M, K = A.shape
  eta = util.softmax(psi)   # Kx1
  phi = np.dot(Z, eta) # Kx1

  # Prevent division by zero when `omega = 1`, in `(1 - phi_muts*omega)`.
  # Also, because of numerical weirdness, phi_muts can be very slightly over 1.
  #
  # Setting delta to 1e-20 resulted in phi being exactly 1 sometimes. That's ... concerning.
  # Having it as 1e-10 is fine, though.
  # Observe the weirdness:
  #     In [3]: 1 - 1e-20 == 1
  #     Out[3]: True
  delta = 1e-10
  phi = np.minimum(1 - delta, phi)
  phi = np.maximum(delta,     phi)
  phi_muts = np.dot(A, phi) # Mx1
  assert np.all(phi_muts > 0) and np.all(phi_muts < 1)

  dlogpx_dphi_muts = (var_reads / phi_muts) - (omega*ref_reads)/(1 - phi_muts*omega) # Mx1
  dlogpx_dphi = np.dot(A.T, dlogpx_dphi_muts) # Kx1
  assert dlogpx_dphi.shape == (K,)

  psi_tiled = _tile_rows(psi, K)
  assert psi_tiled.shape == (K, K)
  # Note: psi_sums[i,j] = psi[i] + psi[j]
  psi_sums = psi_tiled + psi_tiled.T
  logsumexp_psi = logsumexp(psi)
  # Note: deta_dpsi[i,j] = deta[i] / dpsi[j]
  deta_dpsi = -np.exp(psi_sums - 2*logsumexp_psi) # KxK
  np.fill_diagonal(deta_dpsi, (np.exp(psi + logsumexp_psi) - np.exp(2*psi)) / np.exp(2*logsumexp_psi))
  assert deta_dpsi.shape == (K, K)

  # beta[i,j] = sum of `deta_dpsi` values in column `j` that are relevant for phi
  # `i` (i.e., where the entries in rows `a` where `a` is a descendant of node
  # `i`)
  beta = np.dot(Z, deta_dpsi)        # KxK
  grad = np.dot(dlogpx_dphi.T, beta) # 1xK

  assert grad.shape == (K,)
  assert not np.any(np.isnan(grad))
  assert not np.any(np.isinf(grad))
  return grad

@njit
def _rprop(var_reads, ref_reads, omega, A, Z, psi, iterations, convergence_threshold):
  # For explanation of `rprop`, see
  # http://www.cs.toronto.edu/~tijmen/csc321/slides/lecture_slides_lec6.pdf.
  step = 1e-4 * np.ones(psi.shape)
  last_sign = np.ones(psi.shape)

  for T in range(iterations):
    grad = _calc_grad(var_reads, ref_reads, omega, A, Z, psi)
    assert not np.any(np.isnan(grad))
    #grad_num = _calc_grad_numerical(var_reads, ref_reads, omega, A, Z, psi)
    #print(grad, grad_num)
    if np.max(np.abs(grad)) < convergence_threshold:
      break
    sign = np.ones(psi.shape)
    sign[grad < 0] = -1
    sign_agree = sign == last_sign
    sign_disagree = np.logical_not(sign_agree)
    last_sign = sign

    delta = 1.2*sign_agree + 0.5*sign_disagree
    # These limits are chosen according to the slides above.
    step *= delta
    step = np.minimum(50, step)
    step = np.maximum(1e-6, step)
    psi += sign * step

  return psi

@njit
def logsumexp(V):
  B = np.max(V)
  log_sum = B + np.log(np.sum(np.exp(V - B)))
  return log_sum

@njit
def _fit_eta_S(eta_S, var_reads_S, ref_reads_S, omega_S, A, Z, method, iterations):
  psi_S = np.log(eta_S)

  if method == 'rprop':
    psi_S = _rprop(
      var_reads_S,
      ref_reads_S,
      omega_S,
      A,
      Z,
      psi_S,
      iterations,
      convergence_threshold=1e-5,
    )
  else:
    raise Exception('Unknown psi fitter')

  eta_S = util.softmax(psi_S)
  return eta_S

@njit
def _calc_llh(var_reads, ref_reads, omega, A, Z, psi):
  K = len(psi)
  assert Z.shape == (K,K)
  assert var_reads.shape == ref_reads.shape == omega.shape

  eta = util.softmax(psi)
  phi = np.dot(Z, eta) # Kx1
  var_phis = np.dot(A, phi)
  logp = binom.logpmf(var_reads, ref_reads + var_reads, var_phis * omega)
  return np.sum(logp)

@njit
def _calc_grad_numerical(var_reads, ref_reads, omega, A, Z, psi):
  _, K = A.shape
  delta = 1e-10
  grad = np.zeros(K)
  for k in range(K):
    psi_prime = np.copy(psi)
    psi_prime[k] += delta
    llh1 = _calc_llh(var_reads, ref_reads, omega, A, Z, psi)
    llh2 = _calc_llh(var_reads, ref_reads, omega, A, Z, psi_prime)
    grad[k] = (llh2 - llh1)/delta
  return grad
