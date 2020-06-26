import argparse
import numpy as np
import pickle

import mutstat

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import resultserializer
import util

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--tree-index', type=int, default=0)
  parser.add_argument('results_fn')
  parser.add_argument('baseline_mutdist_fn')
  args = parser.parse_args()

  results = resultserializer.Results(args.results_fn)

  clusters = [[]] + results.get('clusters')
  vids, membership = util.make_membership_mat(clusters)
  mphi = np.dot(membership, results.get('phi')[args.tree_index])

  baseline = mutstat.Mutstat(stats=mphi, vids=vids, assays=results.get('sampnames'))
  mutstat.write(baseline, args.baseline_mutdist_fn)

if __name__ == '__main__':
  main()
