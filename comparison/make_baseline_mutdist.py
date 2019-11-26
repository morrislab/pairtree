import argparse
import numpy as np
import pickle

import evalutil

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sim_data_fn')
  parser.add_argument('baseline_mutdist_fn')
  args = parser.parse_args()

  with open(args.sim_data_fn, 'rb') as dataf:
    simdata = pickle.load(dataf)

  clusters = [[]] + simdata['clusters']
  vids, membership = evalutil.make_membership_mat(clusters)
  mphi = np.dot(membership, simdata['phi'])
  np.savez_compressed(args.baseline_mutdist_fn, phi=mphi, vids=vids, assays=simdata['sampnames'])

if __name__ == '__main__':
  main()
