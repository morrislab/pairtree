import numpy as np
import argparse
import scipy.stats

import inputparser
import clustermaker
import phi_fitter
import common

MIN_FLOAT = np.finfo(np.float).min

def calc_binom_params(supervars):
  svids = common.extract_vids(supervars)
  V = np.array([supervars[svid]['var_reads'] for svid in svids])
  R = np.array([supervars[svid]['ref_reads'] for svid in svids])
  omega_v = np.array([supervars[svid]['omega_v'] for svid in svids])
  assert np.all(omega_v == 0.5)

  N = V + R
  return (V, N, omega_v)

def _calc_llh_phi_binom(phi, supervars):
  V, N, omega_v = calc_binom_params(supervars)
  K, S = phi.shape
  for arr in V, N, omega_v:
    assert arr.shape == (K-1, S)

  assert np.allclose(1, phi[0])
  P = omega_v * phi[1:]

  phi_llh = scipy.stats.binom.logpmf(V, N, P)
  phi_llh = np.sum(phi_llh)
  assert not np.isnan(phi_llh)
  # Prevent LLH of -inf.
  phi_llh = np.maximum(phi_llh, MIN_FLOAT)
  return phi_llh

def calc_beta_params(supervars):
  svids = common.extract_vids(supervars)
  V = np.array([supervars[svid]['var_reads'] for svid in svids])
  R = np.array([supervars[svid]['ref_reads'] for svid in svids])
  omega_v = np.array([supervars[svid]['omega_v'] for svid in svids])
  assert np.all(omega_v == 0.5)

  # Since these are supervars, we can just take 2*V and disregard omega_v, since
  # supervariants are always diploid (i.e., omega_v = 0.5).
  alpha = 2*V + 1
  # Must ensure beta is > 0.
  beta = np.maximum(1, R - V + 1)
  assert np.all(alpha > 0) and np.all(beta > 0)
  return (alpha, beta)

def _calc_llh_phi_beta(phi, supervars):
  alpha, beta =  calc_beta_params(supervars)
  K, S = phi.shape
  assert alpha.shape == beta.shape == (K-1, S)
  assert np.allclose(1, phi[0])
  phi_llh = scipy.stats.beta.logpdf(phi[1:,:], alpha, beta)
  phi_llh = np.sum(phi_llh)

  # I had NaNs creep into my LLH when my alpha and beta params were invalid
  # (i.e., when I had elements of beta that were <= 0).
  assert not np.isnan(phi_llh)
  # Prevent LLH of -inf.
  phi_llh = np.maximum(phi_llh, MIN_FLOAT)
  return phi_llh

def _adj2parents(adj):
  adj = np.copy(adj)
  np.fill_diagonal(adj, 0)
  return np.argmax(adj[:,1:], axis=0)

def _parents2adj(parents):
  M = len(parents) + 1
  adjm = np.eye(M)
  adjm[parents, range(1, M)] = 1
  return adjm

def print_init(supervars, adj):
  svids = common.extract_vids(supervars)
  R = np.array([supervars[svid]['ref_reads'] for svid in svids])
  V = np.array([supervars[svid]['var_reads'] for svid in svids])
  T = R + V
  omega = np.array([supervars[svid]['omega_v'] for svid in svids])
  M, S = T.shape

  phi_hat = V / (omega * T)
  phi_hat = np.insert(phi_hat, 0, 1, axis=0)
  print('parents', _adj2parents(adj))
  print('V')
  print(V)
  print('T')
  print(T)
  print()
  print_method('phi_hat', phi_hat, supervars)
  print()

def print_method(method, phi, supervars):
  llh_binom = _calc_llh_phi_binom(phi, supervars)
  llh_beta = _calc_llh_phi_beta(phi, supervars)
  print(f'{method} llh_binom={llh_binom:.3f} llh_beta={llh_beta:.3f}')
  print(phi)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.ssm_fn)
  params = inputparser.load_params(args.params_fn)
  clusters = params['clusters']

  supervars = clustermaker.make_cluster_supervars(clusters, variants)
  superclusters = clustermaker.make_superclusters(supervars)
  # Add empty initial cluster, which serves as tree root.
  superclusters.insert(0, [])
  M = len(superclusters)

  iterations = 1000
  parallel = 0

  parents = [[0, 0, 0], [0, 1, 2]]
  for P in parents:
    adj = _parents2adj(P)
    print_init(supervars, adj)
    for method in ('projection', 'rprop', 'graddesc'):
      phi, eta = phi_fitter._fit_phis(adj, superclusters, supervars, method, iterations, parallel)
      # Sometimes the `projection` fitter will return zeros, which result in an
      # LLH of -inf if the number of variant reads `V` is non-zero, since
      # `Binom(X=V > 0, | N=V+R, p=0) = 0`. To avoid this, set a floor of 1e-6
      # on phi values.
      phi = np.maximum(1e-6, phi)
      print_method(method, phi, supervars)
      print()

if __name__ == '__main__':
  main()
