import argparse
import pickle
import numpy as np

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import evalutil

def write_mutphi(cluster_phi, clusters, mutphifn):
  vids, membership = evalutil.make_membership_mat(clusters)
  mphi = np.dot(membership, cluster_phi)
  mutphi = evalutil.Mutphi(vids=vids, phi=mphi)
  evalutil.save_mutphi(mutphi, mutphifn)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sim_data_fn')
  parser.add_argument('mutphi_fn')
  args = parser.parse_args()

  with open(args.sim_data_fn, 'rb') as dataf:
    simdata = pickle.load(dataf)
  clusters = [[]] + simdata['clusters']

  write_mutphi(simdata['phi'], clusters, args.mutphi_fn)

if __name__ == '__main__':
  main()
