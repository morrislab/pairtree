import numpy as np
import common
import ctypes
import numpy.ctypeslib as npct
import subprocess
import os
import sys

MIN_VARIANCE = 1e-4

def _convert_adjm_to_adjlist(adjm):
  adjm = np.copy(adjm)
  assert np.all(np.diag(adjm) == 1)
  # Make undirected.
  adjm += adjm.T
  np.fill_diagonal(adjm, 0)
  assert np.all(np.logical_or(adjm == 0, adjm == 1))

  adjl = []
  for I, row in enumerate(adjm):
    adjl.append(np.flatnonzero(row))

  return adjl

def fit_etas(adj, superclusters, supervars):
  svids = common.extract_vids(supervars)
  R = np.array([supervars[svid]['ref_reads'] for svid in svids], dtype=np.float64)
  V = np.array([supervars[svid]['var_reads'] for svid in svids], dtype=np.float64)
  T = R + V
  omega = np.array([supervars[svid]['omega_v'] for svid in svids], dtype=np.float64)
  M, S = T.shape

  # To prevent the true_divide error, we use numpy divide with the condition to only perform division if omega * T is not zero.
  phi_hat = np.divide(V, omega * T, out=np.zeros_like(V), where=(omega * T)!=0)
  phi_hat = np.maximum(0, phi_hat)
  phi_hat = np.minimum(1, phi_hat)
  phi_hat = np.insert(phi_hat, 0, 1, axis=0)

  # Make Quaid happy with `V_hat` and `T_hat`. I don't really understand why
  # we're doing this, but I believe that, when `V = 0` and `T` is relatively
  # small, this will result in `var_phi_hat` being larger than it would
  # otherwise be. Without using `V_hat` and `T_hat`, `var_phi_hat` would be
  # zero here, and so would be bumped up to the floor of 1e-8 below.
  V_hat = V + 1
  T_hat = T + 2
  var_phi_hat = V_hat*(1 - V_hat/T_hat) / (T_hat*omega)**2
  var_phi_hat = np.insert(var_phi_hat, 0, MIN_VARIANCE, axis=0)
  var_phi_hat = np.maximum(MIN_VARIANCE, var_phi_hat)
  assert var_phi_hat.shape == (M+1, S)

  eta = np.zeros((M+1, S))
  for sidx in range(S):
    eta[:,sidx] = _fit_eta_S(adj, phi_hat[:,sidx], var_phi_hat[:,sidx])

  assert not np.any(np.isnan(eta))
  assert np.allclose(0, eta[eta < 0])
  eta[eta < 0] = 0
  return eta

def _fit_eta_S_nancheck(adj, phi_hat, var_phi_hat):
  # sim_K100_S100_T50_M1000_G100_run1 revealed a weird issue where a particular
  # call to _project_ppm() will return NaNs. This is reproducible insofar as
  # running Pairtree again with the same seed ("1") will cause this failure at
  # exactly the same point for exactly the same tree sample. In this instance,
  # 101 of the 10,000 entries of the `eta` matrix will be NaN. However, calling
  # _project_ppm() immediately after with the same inputs will work fine.  Out
  # of 705 simulations, each of which sampled 3000 trees (meaning >2M trees
  # sampled), this was the only case where I saw this failure.
  #
  # To work around this, if the first call to _project_ppm() returns NaNs, make
  # two additional attempts.
  max_attempts = 3
  for attempt in range(max_attempts):
    eta = _fit_eta_S_ctypes(adj, phi_hat, var_phi_hat)
    if not np.any(np.isnan(eta)):
      return eta
    print('eta contains NaN, retrying ...', file=sys.stderr)
  raise Exception('eta still contains NaN after %s attempt(s)' % max_attempts)

def _project_ppm(adjm, phi_hat, var_phi_hat, root):
  assert phi_hat.ndim == var_phi_hat.ndim == 1
  inner_flag = 0
  compute_eta = 1
  M = len(phi_hat)
  S = 1
  eta = np.empty(M, dtype=np.double)
  assert M >= 1
  assert var_phi_hat.shape == (M,)

  # This is called `gamma_init` in the C code from B&J.
  gamma_init = var_phi_hat
  phi_hat = phi_hat / gamma_init

  adjl = _convert_adjm_to_adjlist(adjm)
  deg = np.array([len(children) for children in adjl], dtype=np.short)
  adjl_mat = np.zeros((M,M), dtype=np.short)
  for rowidx, row in enumerate(adjl):
    adjl_mat[rowidx,:len(row)] = row

  # Method signature:
  # realnumber tree_cost_projection(
  #   shortint inner_flag,
  #   shortint compute_M_flag,
  #   realnumber *M,
  #   shortint num_nodes,
  #   shortint T,
  #   realnumber *data,
  #   realnumber gamma_init[],
  #   shortint root_node,
  #   edge *tree,
  #   shortint *adjacency_mat,
  #   shortint *final_degrees,
  #   shortint *adj_list
  # );
  c_double_p = ctypes.POINTER(ctypes.c_double)
  c_short_p = ctypes.POINTER(ctypes.c_short)

  # Ensure arrays are C-contiguous.
  eta = np.require(eta, requirements='C')
  phi_hat = np.require(phi_hat, requirements='C')
  gamma_init = np.require(gamma_init, requirements='C')
  deg = np.require(deg, requirements='C')
  adjl_mat = np.require(adjl_mat, requirements='C')

  cost = _project_ppm.tree_cost_projection(
    inner_flag,
    compute_eta,
    eta,
    M,
    S,
    phi_hat,
    gamma_init,
    root,
    None,
    None,
    deg,
    adjl_mat,
  )
  return eta

def _init_project_ppm():
  real_arr_1d = npct.ndpointer(dtype=np.float64, ndim=1, flags='C')
  short_arr_1d = npct.ndpointer(dtype=ctypes.c_short, ndim=1, flags='C')
  short_arr_2d = npct.ndpointer(dtype=ctypes.c_short, ndim=2, flags='C')
  class Edge(ctypes.Structure):
    _fields_ = [('first', ctypes.c_short), ('second', ctypes.c_short)]
  c_edge_p = ctypes.POINTER(Edge)
  c_short_p = ctypes.POINTER(ctypes.c_short)

  lib_path = os.path.join(os.path.dirname(__file__), 'projectppm', 'bin', 'libprojectppm.so')
  assert lib_path is not None, 'Could not find libprojectppm'
  lib = ctypes.cdll.LoadLibrary(lib_path)
  func = lib.tree_cost_projection
  func.argtypes = [
    ctypes.c_short,
    ctypes.c_short,
    real_arr_1d,
    ctypes.c_short,
    ctypes.c_short,
    real_arr_1d,
    real_arr_1d,
    ctypes.c_short,
    c_edge_p,
    c_short_p,
    short_arr_1d,
    short_arr_2d,
  ]
  func.restype = ctypes.c_double
  _project_ppm.tree_cost_projection = func
_init_project_ppm()

def _fit_eta_S_ctypes(adj, phi_hat, var_phi_hat):
  assert phi_hat.ndim == var_phi_hat.ndim == 1
  M = len(phi_hat)
  assert M >= 1
  assert var_phi_hat.shape == (M,)
  root = 0

  eta = _project_ppm(adj, phi_hat, var_phi_hat, root)
  return eta

def _prepare_subprocess_inputs(adjm, phi, prec_sqrt):
  _arr2floatstr = lambda arr: ' '.join([f'{E:.10f}' for E in np.array(arr).flatten()])
  _arr2intstr = lambda arr: ' '.join([f'{E:d}' for E in np.array(arr).flatten()])

  assert phi.ndim == prec_sqrt.ndim == 1
  M = len(phi)
  assert prec_sqrt.shape == (M,)
  assert M >= 1

  root = 0
  adjl = _convert_adjm_to_adjlist(adjm)
  deg = [len(children) for children in adjl]
  should_calc_eta = 1

  calcphi_input = [
    _arr2intstr((M, 1)),
    _arr2floatstr(phi),
    _arr2floatstr(prec_sqrt),
    str(root),
    _arr2intstr(deg),
  ]
  calcphi_input += [_arr2intstr(row) for row in adjl]
  calcphi_input += [
    str(should_calc_eta),
    '' # There will be a trailing newline, as B&J's code requires.
  ]

  joined = '\n'.join(calcphi_input)
  return joined

def _fit_eta_S_subprocess(adj, phi_hat_S, var_phi_hat_S):
  prec_sqrt_S = 1 / np.sqrt(var_phi_hat_S)
  calcphi_input = _prepare_subprocess_inputs(adj, phi_hat_S, prec_sqrt_S)

  result = subprocess.run([os.path.join(os.path.dirname(__file__), 'projectppm', 'bin', 'projectppm')], input=calcphi_input, capture_output=True, encoding='UTF-8')
  result.check_returncode()
  lines = result.stdout.strip().split('\n')
  assert len(lines) == 2
  cost = float(lines[0])
  eta_S = np.array([float(E) for E in lines[1].split(' ')])
  assert len(eta_S) == len(phi_hat_S)

  return eta_S

_fit_eta_S = _fit_eta_S_nancheck
