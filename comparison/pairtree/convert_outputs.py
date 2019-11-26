import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import resultserializer

def convert(results):
  N = len(results['struct'])
  clusters = [[]] + list(results['clusters'])
  neutree = {
    'structs': results['struct'],
    'phis': results['phi'],
    'counts': results['count'],
    'logscores': results['llh'],
    'clusterings': [clusters for idx in range(N)],
    'garbage': [[] for idx in range(N)],
  }
  return neutree

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('pairtree_results_fn')
  parser.add_argument('neutree_fn')
  args = parser.parse_args()

  results = resultserializer.load(args.pairtree_results_fn)
  neutree = convert(results)
  np.savez_compressed(args.neutree_fn, **neutree)

if __name__ == '__main__':
  main()
