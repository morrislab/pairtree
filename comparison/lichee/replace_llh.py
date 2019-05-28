import argparse
import json
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import common
import inputparser

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('params_fn')
  parser.add_argument('pairtree_results_fn')
  args = parser.parse_args()

  params = inputparser.load_params(args.params_fn)
  pairtree_results = np.load(args.pairtree_results_fn)
  pairtree_results = {K: pairtree_results[K] for K in pairtree_results}

  for adjm, adjl in zip(pairtree_results['adjm'], params['structures']):
    adjl = inputparser.load_structure(adjl)
    implied_adjm = common.convert_adjlist_to_adjmatrix(adjl)
    assert np.array_equal(implied_adjm, adjm)

  pairtree_results['llh'] = -1*np.array(params['scores'])
  np.savez_compressed(args.pairtree_results_fn, **pairtree_results)

if __name__ == '__main__':
  main()
