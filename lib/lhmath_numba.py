import ctypes
from numba.extending import get_cython_function_address
import numba
import binom
import numpy as np
from common import Models
from scipy import LowLevelCallable
import util

# TODO: Can I just replace this with `util.lbeta`? Would it be faster / less bullshit?
def _make_betainc():
  addr = get_cython_function_address('scipy.special.cython_special', 'betainc')
  betainc_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
  func = betainc_type(addr)
  return func

# Signature: betacdf(A, B, X)
betacdf = _make_betainc()

# Specifically for scalars.
@numba.njit
def _binom_logpmf(X, N, P):
  val = binom.logpmf(
    np.array([X], dtype=np.int64),
    np.array([N], dtype=np.int64),
    np.array([P], dtype=np.float64),
  )
  return val[0]

@numba.njit
def _make_lower(phi1, midx):
  if midx == Models.A_B:
    return 0
  elif midx == Models.B_A:
    return phi1
  elif midx == Models.diff_branches:
    return 0
  else:
    return np.nan

@numba.njit
def _make_upper(phi1, midx):
  if midx == Models.A_B:
    return phi1
  elif midx == Models.B_A:
    return 1
  elif midx == Models.diff_branches:
    return 1 - phi1
  else:
    return np.nan

def _integral_separate_clusters(args):
  phi1, V1_var_reads, V1_ref_reads, V1_omega, V2_var_reads, V2_ref_reads, V2_omega, midx, logsub = args

  V1_total_reads = V1_var_reads + V1_ref_reads
  logP = _binom_logpmf(
    V1_var_reads,
    V1_total_reads,
    V1_omega * phi1,
  )
  upper = _make_upper(phi1, midx)
  lower = _make_lower(phi1, midx)

  A = V2_var_reads + 1
  B = V2_ref_reads + 1
  betainc_upper = betacdf(A, B, V2_omega * upper)
  betainc_lower = betacdf(A, B, V2_omega * lower)
  if util.isclose(betainc_upper, betainc_lower):
    return 0
  logP += np.log(betainc_upper - betainc_lower)
  logP -= logsub

  return np.exp(logP)

def _integral_same_cluster(args):
  phi1, V1_var_reads, V1_ref_reads, V1_omega, V2_var_reads, V2_ref_reads, V2_omega, logsub = args

  V1_total_reads = V1_var_reads + V1_ref_reads
  V2_total_reads = V2_var_reads + V2_ref_reads

  X = np.array([V1_var_reads,    V2_var_reads])
  N = np.array([V1_total_reads,  V2_total_reads])
  P = np.array([V1_omega * phi1, V2_omega * phi1])

  B = binom.logpmf(X, N, P)
  logP = B[0] + B[1] - logsub
  return np.exp(logP)

# See:
# https://stackoverflow.com/questions/51109429
# https://stackoverflow.com/questions/49683653
def _make_jitted_integrand(integrand):
  jitted = numba.njit(integrand)
  # This is the function that scipy.integrate.quad can call.
  @numba.cfunc(numba.types.float64(numba.types.intc, numba.types.CPointer(numba.types.float64)))
  def integrand_cfunc(N, XX):
    vals = numba.carray(XX, N)
    return jitted(vals)
  return LowLevelCallable(integrand_cfunc.ctypes)

integral_separate_clusters = _make_jitted_integrand(_integral_separate_clusters)
integral_same_cluster      = _make_jitted_integrand(_integral_same_cluster)
