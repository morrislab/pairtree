import argparse
import numpy as np

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import neutree
import pastri_util
import util

def convert(sampid, params_fn, trees_fn, neutree_fn):
  adjms, llhs, phis, clusterings = pastri_util.load_results(sampid, params_fn, trees_fn)
  if len(adjms) == 0:
    return
  structs = [util.convert_adjmatrix_to_parents(A) for A in adjms]
  N = len(structs)
  ntree = neutree.Neutree(
    structs = structs,
    phis = phis,
    counts = np.ones(N),
    logscores = llhs,
    clusterings = clusterings,
    garbage = [[] for idx in range(N)],
  )
  neutree.save(ntree, neutree_fn)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sampid')
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('trees_fn')
  parser.add_argument('neutree_fn')
  args = parser.parse_args()

  convert(args.sampid, args.pairtree_params_fn, args.trees_fn, args.neutree_fn)

if __name__ == '__main__':
  main()
