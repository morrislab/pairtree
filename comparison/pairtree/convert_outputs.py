import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import resultserializer
import inputparser
import neutree

def convert(results, garbage):
  N = len(results['struct'])
  clusters = [[]] + list(results['clusters'])
  ntree = neutree.Neutree(
    structs = results['struct'],
    phis = results['phi'],
    counts = results['count'],
    logscores = results['llh'],
    clusterings = [clusters for idx in range(N)],
    garbage = [garbage for idx in range(N)],
  )
  return ntree

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('pairtree_results_fn')
  parser.add_argument('params_fn')
  parser.add_argument('neutree_fn')
  args = parser.parse_args()

  results = resultserializer.load(args.pairtree_results_fn)
  params = inputparser.load_params(args.params_fn)
  ntree = convert(results, params['garbage'])
  neutree.save(ntree, args.neutree_fn)

if __name__ == '__main__':
  main()
