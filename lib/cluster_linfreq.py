import numpy as np
from numba import njit
import util
import inputparser

@njit
def _calc_llh(V, T_prime, Z, phi_alpha0, phi_beta0, logconc):
  uniq_Z  = set(list(Z))
  C = len(uniq_Z)
  #assert uniq_Z == set(range(C))
  N = len(Z)
  cluster_sizes = np.array([np.sum(Z == c) for c in range(C)])
  assert np.sum(cluster_sizes) == N

  llh  = C * logconc
  llh += np.sum(util.logfactorial(cluster_sizes - 1))
  llh -= np.sum(np.log(np.exp(logconc) + np.arange(N)))

  for cidx in range(C):
    members = np.flatnonzero(Z == cidx)
    V_summed = np.sum(V[members], axis=0)
    T_summed = np.sum(T_prime[members], axis=0)
    R_summed = T_summed - V_summed

    binom_coef = util.log_N_choose_K(T_prime[members], V[members])
    llh_c = np.sum(binom_coef, axis=0)
    llh_c += util.lbeta(V_summed + phi_alpha0, R_summed + phi_beta0)
    llh_c -= util.lbeta(phi_alpha0, phi_beta0)

    llh += np.sum(llh_c)

  return llh

@njit
def _calc_cweight(V, T_prime, phi_alpha0, phi_beta0, vidx, members):
  V_summed = np.sum(V[members], axis=0)
  R_summed = np.sum(T_prime[members], axis=0) - V_summed
  Vi = V[vidx]
  Ti = T_prime[vidx]
  Ri = Ti - Vi

  result  = util.log_N_choose_K(Ti, Vi)
  result += util.lbeta(V_summed + Vi + phi_alpha0, R_summed + Ri + phi_beta0)
  result -= util.lbeta(V_summed + phi_alpha0, R_summed + phi_beta0)

  return np.log(len(members)) + np.sum(result)

@njit
def _calc_new_cluster_weight(V, T_prime, phi_alpha0, phi_beta0, vidx, logconc):
  Vi = V[vidx]
  Ti = T_prime[vidx]
  Ri = Ti - Vi

  result  = util.log_N_choose_K(Ti, Vi)
  result += util.lbeta(Vi + phi_alpha0, Ri + phi_beta0)
  # This is subtracted from *each element* of `result`, which is important.
  result -= util.lbeta(phi_alpha0, phi_beta0)

  return logconc + np.sum(result)

@njit
def _compute_cweights_full(V, T_prime, Z, phi_alpha0, phi_beta0, vidx, logconc, C):
  mask = np.ones(len(V), dtype=np.bool_)
  mask[vidx] = 0
  llh_excluded = _calc_llh(V[mask], T_prime[mask], Z[mask], phi_alpha0, phi_beta0, logconc)

  cweights_full = np.full(C + 1, np.nan)
  Z_prime = np.copy(Z)

  for cidx in range(C + 1):
    Z_prime[vidx] = cidx
    llh_included = _calc_llh(V, T_prime, Z_prime, phi_alpha0, phi_beta0, logconc)
    cweights_full[cidx] = llh_included - llh_excluded

  return cweights_full

@njit
def _do_gibbs_iter(V, T_prime, phi_alpha0, phi_beta0, logconc, C, Z, check_full_llh=False):
  N = len(Z)
  Z = np.copy(Z)

  for vidx in range(N):
    old_cluster = Z[vidx]
    Z[vidx] = -1
    if not np.any(Z == old_cluster):
      # If `vidx` was the only member, remove this cluster.
      # Do so by moving the last cluster to its index.
      # (These operations are valid even if `highest == old_cluster`.)
      highest = C - 1
      Z[Z == highest] = old_cluster
      C -= 1

    # `cweights`: LLHs of each cluster destination for `vidx`
    # cweights[C] = LLH of adding new cluster
    cweights = np.empty(C + 1)
    # Consider every possible destination.
    for cidx in range(C):
      members = np.flatnonzero(Z == cidx)
      cweights[cidx] = _calc_cweight(V, T_prime, phi_alpha0, phi_beta0, vidx, members)
    # Consider adding a new cluster.
    cweights[C] = _calc_new_cluster_weight(V, T_prime, phi_alpha0, phi_beta0, vidx, logconc)
    cweights -= np.log(np.exp(logconc) + N - 1)

    if check_full_llh:
      cweights_full = _compute_cweights_full(V, T_prime, Z, phi_alpha0, phi_beta0, vidx, logconc, C)
      assert np.all(util.isclose(cweights, cweights_full))

    cprobs = util.softmax(cweights)
    new_cluster = util.sample_multinom(cprobs)
    Z[vidx] = new_cluster
    if new_cluster == C:
      C += 1

  llh = _calc_llh(V, T_prime, Z, phi_alpha0, phi_beta0, logconc)
  return (C, Z, llh)

def cluster(variants, raw_clusters, logconc, iters, seed, progress_queue):
  np.random.seed(seed % 2**32)
  vids, V, T, T_prime, omega = inputparser.load_read_counts(variants)

  # M: number of variants
  # C: number of clusters
  # Z: cluster assignments
  M = len(V)
  C = 1
  Z = np.zeros(M, np.int32)

  # Beta distribution prior for phi
  phi_alpha0 = 1.
  phi_beta0 = 1.

  clusterings = []
  llhs = []

  for I in range(iters):
    if progress_queue is not None:
      progress_queue.put(I)
    C, Z, llh = _do_gibbs_iter(V, T_prime, phi_alpha0, phi_beta0, logconc, C, Z, check_full_llh=False)
    clusterings.append(Z)
    llhs.append(llh)

  return (vids, np.array(clusterings), np.array(llhs))
