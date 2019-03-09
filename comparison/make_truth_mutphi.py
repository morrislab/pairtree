import argparse
import pickle
import numpy as np

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import mutphi

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sim_data_fn')
  parser.add_argument('ssm_fn')
  parser.add_argument('mutphi_fn')
  args = parser.parse_args()

  with open(args.sim_data_fn, 'rb') as dataf:
    simdata = pickle.load(dataf)
  clusters = [[]] + simdata['clusters']

  llh = 0
  weight_trees_by = 'uniform'
  mphi = mutphi.calc_mutphi([simdata['phi']], [llh], [clusters], weight_trees_by, args.ssm_fn)
  mutphi.write_mutphi(mphi, args.mutphi_fn)

if __name__ == '__main__':
  main()
