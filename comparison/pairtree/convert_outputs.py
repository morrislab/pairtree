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
  N = len(results.get('struct'))
  clusters = [[]] + list(results.get('clusters'))
  ntree = neutree.Neutree(
    structs = results.get('struct'),
    phis = results.get('phi'),
    counts = results.get('count'),
    logscores = results.get('llh'),
    clusterings = [clusters for idx in range(N)],
    garbage = garbage,
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

  results = resultserializer.Results(args.pairtree_results_fn)
  params = inputparser.load_params(args.params_fn)
  ntree = convert(results, params['garbage'])
  neutree.save(ntree, args.neutree_fn)

if __name__ == '__main__':
  main()
