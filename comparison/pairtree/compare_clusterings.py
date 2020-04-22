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
import clustermaker

def convert_clustering_to_assignment(clusters):
  mapping = {vid: cidx for cidx, cluster in enumerate(clusters) for vid in cluster}
  vids = common.sort_vids(mapping.keys())
  assign = np.array([mapping[vid] for vid in vids])
  return (vids, assign)

def extract_assignment(paramsfn):
  params = inputparser.load_params(paramsfn)
  clusters = params['clusters']
  C = len(clusters)
  vids, assign = convert_clustering_to_assignment(clusters)
  return (C, vids, assign)

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
  llh1 = clustermaker.calc_llh(V, T_prime, assign1, hparams['phi_alpha0'], hparams['phi_beta0'], hparams['conc'])
  llh2 = clustermaker.calc_llh(V, T_prime, assign2, hparams['phi_alpha0'], hparams['phi_beta0'], hparams['conc'])
  nlglh1 = -llh1 / (M*S*np.log(2))
  nlglh2 = -llh2 / (M*S*np.log(2))

  homo, comp, vm = sklearn.metrics.homogeneity_completeness_v_measure(assign1, assign2)
  ami = sklearn.metrics.adjusted_mutual_info_score(assign1, assign2)
  print(C1, C2, llh1, llh2, nlglh1, nlglh2, homo, comp, vm, ami, sep=',')

if __name__ == '__main__':
  main()
