import argparse
import sklearn.metrics
import numpy as np
from numba import njit

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import common
import inputparser
import util

def convert_clustering_to_assignment(clusters):
  mapping = {vid: cidx for cidx, cluster in enumerate(clusters) for vid in cluster}
  vids = common.sort_vids(mapping.keys())
  assign = np.array([mapping[vid] for vid in vids])
  return (vids, assign)

def extract_assignment(paramsfn):
  params = inputparser.load_params(paramsfn)
  C = len(clusters)
  vids, assign = convert_clustering_to_assignment(params['clusters'])
  return (C, vids, assign)

# The two `_calc_llh` functions are copied  from `bin/clustervars`. I should
# move them to a central module once I decide what the right LLH is.
@njit
def _calc_llh_pants(V, T_prime, Z, phi_alpha0, phi_beta0, conc):
  uniq_Z  = set(list(Z))
  C = len(uniq_Z)
  #assert uniq_Z == set(range(C))
  N = len(Z)
  cluster_sizes = np.array([np.sum(Z == c) for c in range(C)])
  assert np.sum(cluster_sizes) == N

  llh  = np.log(C)
  llh += np.sum(util.logfactorial(cluster_sizes - 1))
  llh -= np.sum(np.log(conc + np.arange(N)))

  for cidx in range(C):
    members = np.flatnonzero(Z == cidx)
    alpha_c = phi_alpha0 + np.sum(V[members],                    axis=0)
    beta_c  = phi_beta0  + np.sum(T_prime[members] - V[members], axis=0)
    for vidx in members:
      llh += np.sum(util.beta_binom_logpmf(V[vidx], T_prime[vidx], alpha_c, beta_c))
  return llh

@njit
def _calc_llh_socks(V, T_prime, Z, phi_alpha0, phi_beta0, conc):
  uniq_Z  = set(list(Z))
  C = len(uniq_Z)
  #assert uniq_Z == set(range(C))
  N = len(Z)
  cluster_sizes = np.array([np.sum(Z == c) for c in range(C)])
  assert np.sum(cluster_sizes) == N

  llh  = np.log(C)
  llh += np.sum(util.logfactorial(cluster_sizes - 1))
  llh -= np.sum(np.log(conc + np.arange(N)))
  #llh_first = llh

  for cidx in range(C):
    members = np.flatnonzero(Z == cidx)
    V_summed = np.sum(V[members], axis=0)
    T_summed = np.sum(T_prime[members], axis=0)
    R_summed = T_summed - V_summed

    llh_c  = util.log_N_choose_K(T_summed, V_summed)
    llh_c += util.lbeta(V_summed + phi_alpha0, R_summed + phi_beta0)
    llh_c -= util.lbeta(phi_alpha0, phi_beta0)

    llh += np.sum(llh_c)
  #print(llh_first, llh - llh_first, sep='\t')

  return llh

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('ssm_fn')
  parser.add_argument('params1_fn')
  parser.add_argument('params2_fn')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.ssm_fn)
  vids = common.extract_vids(variants)
  V, T, T_prime, omega = inputparser.load_read_counts(variants)
  M, S  = V.shape

  C1, vids1, assign1 = extract_assignment(args.params1_fn)
  C2, vids2, assign2 = extract_assignment(args.params2_fn)
  assert vids1 == vids2 == vids

  hparams = {
    'phi_alpha0': 1.,
    'phi_beta0': 1.,
    'conc': 1e-2,
  }
  _calc_llh = _calc_llh_socks
  llh1 = _calc_llh(V, T_prime, assign1, hparams['phi_alpha0'], hparams['phi_beta0'], hparams['conc'])
  llh2 = _calc_llh(V, T_prime, assign2, hparams['phi_alpha0'], hparams['phi_beta0'], hparams['conc'])
  nlglh1 = -llh1 / (M*S*np.log(2))
  nlglh2 = -llh2 / (M*S*np.log(2))

  homo, comp, vm = sklearn.metrics.homogeneity_completeness_v_measure(assign1, assign2)
  ami = sklearn.metrics.adjusted_mutual_info_score(assign1, assign2)
  print(C1, C2, llh1, llh2, nlglh1, nlglh2, homo, comp, vm, ami, sep=',')

if __name__ == '__main__':
  main()
