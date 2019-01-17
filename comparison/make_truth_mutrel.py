import argparse
import numpy as np
import pickle

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import mutrel
import evalutil

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sim_data_fn')
  parser.add_argument('mutrel_fn')
  args = parser.parse_args()

  with open(args.sim_data_fn, 'rb') as dataf:
    simdata = pickle.load(dataf)
  clusters = [[]] + simdata['clusters']

  mrel = mutrel.make_mutrel_tensor_from_cluster_adj(simdata['adjm'], clusters)
  mrel = evalutil.add_garbage(mrel, simdata['vids_garbage'])
  assert set(mrel.vids) == set(simdata['vids_good'] + simdata['vids_garbage'])
  evalutil.save_sorted_mutrel(mrel, args.mutrel_fn)

if __name__ == '__main__':
  main()
