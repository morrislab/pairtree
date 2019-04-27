import numpy as np
import common
import ctypes
import ctypes.util
import numpy.ctypeslib as npct
import subprocess

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

def fit_phis(adj, superclusters, supervars):
  svids = common.extract_vids(supervars)
  R = np.array([supervars[svid]['ref_reads'] for svid in svids])
  V = np.array([supervars[svid]['var_reads'] for svid in svids])
  T = R + V
  omega = np.array([supervars[svid]['omega_v'] for svid in svids])
  M, S = T.shape

  phi_hat = V / (omega * T)
  var_phi_hat = V*(1 - V/T) / (T*omega)**2

  phi_hat = np.insert(phi_hat, 0, 1, axis=0)
  var_phi_hat = np.maximum(1e-8, var_phi_hat)
  prec_sqrt = np.sqrt(1 / var_phi_hat)
  prec_sqrt = np.insert(prec_sqrt, 0, 1e3, axis=0)
  assert prec_sqrt.shape == (M+1, S)

  eta = np.zeros((M+1, S))
  for sidx in range(S):
    eta[:,sidx] = _fit_phi_S(adj, phi_hat[:,sidx], prec_sqrt[:,sidx])
  assert np.allclose(0, eta[eta < 0])
  eta[eta < 0] = 0

  Z = common.make_ancestral_from_adj(adj)
  phi = np.dot(Z, eta)
  return (phi, eta)

def _project_ppm(adjm, phi_hat, prec_sqrt, root):
  assert phi_hat.ndim == prec_sqrt.ndim == 1
  inner_flag = 1
  compute_eta = 1
  M = len(phi_hat)
  S = 1
  eta = np.empty(M, dtype=np.double)
  assert M >= 1
  assert prec_sqrt.shape == (M,)

  gamma_init = 1 / prec_sqrt**2 # This is now the variance.
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

  lib_path = ctypes.util.find_library('projectppm')
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

def _fit_phi_S_ctypes(adj, phi_hat, prec_sqrt):
  assert phi_hat.ndim == prec_sqrt.ndim == 1
  M = len(phi_hat)
  assert M >= 1
  assert prec_sqrt.shape == (M,)

  root = 0
  eta = _project_ppm(adj, phi_hat, prec_sqrt, root)
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
    ''
  ]

  joined = '\n'.join(calcphi_input)
  return joined

def _fit_phi_S_subprocess(adj, phi_hat_S, prec_sqrt_S):
  calcphi_input = _prepare_subprocess_inputs(adj, phi_hat_S, prec_sqrt_S)

  result = subprocess.run(['projectppm'], input=calcphi_input, capture_output=True, encoding='UTF-8')
  result.check_returncode()
  lines = result.stdout.strip().split('\n')
  assert len(lines) == 2
  cost = float(lines[0])
  eta_S = np.array([float(E) for E in lines[1].split(' ')])
  assert len(eta_S) == len(phi_hat_S)

  return eta_S

_fit_phi_S = _fit_phi_S_ctypes
