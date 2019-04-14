import numpy as np
import common
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

def _prepare_inputs(adjm, phi, var_phi):
  # order='F' means to flatten in column-major order, which Bei's code needs.
  _arr2floatstr = lambda arr: ' '.join([f'{E:.10f}' for E in np.array(arr).flatten(order='C')])
  _arr2intstr = lambda arr: ' '.join([f'{E:d}' for E in np.array(arr).flatten(order='C')])

  phi[:] = 0.3
  phi[-1,:] = 0.2
  var_phi[:] = 0.01
  adjm[:,2] = np.array([0, 1, 1, 0])

  M, S = phi.shape
  prec_sqrt = np.sqrt(1 / var_phi)
  phi = np.insert(phi, 0, 1, axis=0)
  prec_sqrt = np.insert(prec_sqrt, 0, 1e3, axis=0)
  prec_sqrt = np.mean(prec_sqrt, axis=1)
  prec_sqrt /= np.max(prec_sqrt)
  assert prec_sqrt.shape == (M+1,)

  root = 0
  adjl = _convert_adjm_to_adjlist(adjm)
  deg = [len(children) for children in adjl]
  should_calc_eta = '1'

  calcphi_input = [
    _arr2intstr((M+1,S)),
    _arr2floatstr(phi),
    _arr2floatstr(prec_sqrt),
    str(root),
    _arr2intstr(deg),
  ]
  calcphi_input += [_arr2intstr(row) for row in _convert_adjm_to_adjlist(adjm)]
  calcphi_input += [
    should_calc_eta,
    ''
  ]

  joined = '\n'.join(calcphi_input)
  print(joined)
  return joined

def fit_phis(adj, superclusters, supervars, iterations, parallel):
  svids = common.extract_vids(supervars)
  R = np.array([supervars[svid]['ref_reads'] for svid in svids])
  V = np.array([supervars[svid]['var_reads'] for svid in svids])
  T = R + V
  omega = np.array([supervars[svid]['omega_v'] for svid in svids])
  M, S = T.shape
  print(M, S)

  phi_hat = V / (omega * T)
  var_phi_hat = V*(1 - V/T) / (T*omega)**2
  var_phi_hat = np.maximum(1e-2, var_phi_hat)
  calcphi_input = _prepare_inputs(adj, phi_hat, var_phi_hat)

  result = subprocess.run(['calcphi'], input=calcphi_input, capture_output=True, encoding='UTF-8')
  result.check_returncode()
  lines = result.stdout.strip().split('\n')
  assert len(lines) == 2
  cost = float(lines[0])
  eta = np.array([float(E) for E in lines[1].split(' ')])
  print(eta)
  eta = eta.reshape((S, M+1)).T
  print(eta)
  import sys
  sys.exit()

  Z = common.make_ancestral_from_adj(adj)
  phi = np.dot(Z, eta)
  return (phi, eta)
