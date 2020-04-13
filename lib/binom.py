import numpy as np
from numba import njit
import util

@njit
def logpmf(X, N, P):
  assert X.shape == N.shape == P.shape
  # Limit to one-dimension arrays, as indexing an array using a boolean array
  # of >1D (as I do with `results[idxs] = val`) is not supported by Numba.
  assert X.ndim == N.ndim == P.ndim == 1

  P_isclose_0 = util.isclose(0, P)
  P_isclose_1 = util.isclose(1, P)
  X_equals_0 = X == 0
  X_equals_N = X == N
  special = [
    (np.logical_and(P_isclose_0, X_equals_0), 0),
    (np.logical_and(P_isclose_0, np.logical_not(X_equals_0)), -np.inf),
    (np.logical_and(P_isclose_1, X_equals_N), 0),
    (np.logical_and(P_isclose_1, np.logical_not(X_equals_N)), -np.inf),
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
  result[unfilled] = util.log_N_choose_K(Nu, Xu) + Xu*np.log(Pu) + (Nu - Xu)*np.log(1 - Pu)
  assert not np.any(np.isnan(result))

  return result
