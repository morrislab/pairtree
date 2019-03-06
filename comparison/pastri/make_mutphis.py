import argparse
import numpy as np

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import mutphi
import pastri_util

def write_mutphis(sampid, params_fn, trees_fn, weight_trees_by, ssms_fn, mutphi_fn):
  adjms, llhs, phis, clusterings = pastri_util.load_results(sampid, params_fn, trees_fn)
  if len(adjms) == 0:
    return
  mphi = mutphi.calc_mutphi(phis, llhs, clusterings, weight_trees_by, ssms_fn)
  mutphi.write_mutphi(mphi, mutphi_fn)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--weight-trees-by', choices=('llh', 'uniform'), default='uniform')
  parser.add_argument('sampid')
  parser.add_argument('pairtree_ssms_fn')
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('trees_fn')
  parser.add_argument('mutphi_fn')
  args = parser.parse_args()

  write_mutphis(args.sampid, args.pairtree_params_fn, args.trees_fn, args.weight_trees_by, args.pairtree_ssms_fn, args.mutphi_fn)

if __name__ == '__main__':
  main()
