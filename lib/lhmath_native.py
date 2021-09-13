import numpy as np
import scipy.special
import binom
from common import Models

def _make_lower(phi1, midx):
  ret = {
    Models.A_B: 0,
    Models.B_A: phi1,
    Models.diff_branches: 0,
  }[midx]
  return np.broadcast_to(ret, np.array(phi1).shape)

def _make_upper(phi1, midx):
  ret = {
    Models.A_B: phi1,
    Models.B_A: 1,
    Models.diff_branches: 1 - phi1,
  }[midx]
  return np.broadcast_to(ret, np.array(phi1).shape)

def binom_logpmf_scalar(X, N, P):
  return binom.logpmf(
    np.array([X]),
    np.array([N]),
    np.array([P]),
  )[0]

def integral_separate_clusters(phi1, V1, V2, sidx, midx, logsub):
  logP = binom_logpmf_scalar(
    V1.var_reads[sidx],
    V1.total_reads[sidx],
    V1.omega_v[sidx]*phi1
  )
  lower = _make_lower(phi1, midx)
  upper = _make_upper(phi1, midx)

  A = V2.var_reads[sidx] + 1
  B = V2.ref_reads[sidx] + 1
  betainc_upper = scipy.special.betainc(A, B, V2.omega_v[sidx] * upper)
  betainc_lower = scipy.special.betainc(A, B, V2.omega_v[sidx] * lower)
  if np.isclose(betainc_upper, betainc_lower):
    return 0
  logP += np.log(betainc_upper - betainc_lower)
  logP -= logsub

  return np.exp(logP)

def integral_same_cluster(phi1, V1, V2, sidx, logsub):
  binom_params = [(
    V.var_reads[sidx],
    V.total_reads[sidx],
    V.omega_v[sidx]*phi1,
  ) for V in (V1, V2)]
  X, N, P = [np.array(A) for A in zip(*binom_params)]

  B = binom.logpmf(X, N, P)
  # The `log(sqrt(2))` comes from computing the line integral. See theorem 6.3
  # ("Evaluating a Scalar Line Integral") at
  # https://openstax.org/books/calculus-volume-3/pages/6-2-line-integrals.
  # ... except I've now removed it because my external pointed out it was wrong.
  logP = B[0] + B[1] - logsub
  return np.exp(logP)
