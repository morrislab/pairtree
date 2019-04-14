import numpy as np
import math
from numba import njit

@njit
def _logfactorial(X):
  result = np.empty(X.shape)
  for idx in range(len(X)):
    result[idx] = math.lgamma(X[idx] + 1)
  return result

@njit
def _log_N_choose_K(N, K):
  return _logfactorial(N) - (_logfactorial(K) + _logfactorial(N - K))

@njit
def isclose(arr, V, tol=1e-6):
  return np.abs(arr - V) <= tol

@njit
def logpmf(X, N, P):
  assert X.shape == N.shape == P.shape
  assert X.ndim == N.ndim == P.ndim == 1

  P_isclose_0 = isclose(0, P)
  P_isclose_1 = isclose(1, P)
  X_equals_0 = X == 0
  X_equals_N = X == N
  special = [
    (np.logical_and(P_isclose_0, X_equals_0), 0),
    (np.logical_and(P_isclose_0, X_equals_N), -np.inf),
    (np.logical_and(P_isclose_1, X_equals_0), -np.inf),
    (np.logical_and(P_isclose_1, X_equals_N), 0),
  ]

  result = np.full_like(P, np.nan)
  for idx in range(len(special)):
    idxs = special[idx][0]
    val = special[idx][1]
    assert np.all(np.isnan(result[idxs]))
    result[idxs] = val

  unfilled = np.isnan(result)
  Nu = N[unfilled]
  Xu = X[unfilled]
  Pu = P[unfilled]
  result[unfilled] = _log_N_choose_K(Nu, Xu) + Xu*np.log(Pu) + (Nu - Xu)*np.log(1 - Pu)

  return result
