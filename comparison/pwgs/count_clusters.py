import argparse
import numpy as np

import sys
import os
sys.path += [
  os.path.join(os.path.dirname(__file__), '..', '..', 'lib'),
  os.path.expanduser('~/.apps/phylowgs')
]
from pwgsresults.result_loader import ResultLoader
import util

def count_clusters(results):
  tidxs = np.array(sorted(results.tree_summary.keys()))
  llhs = np.array([results.tree_summary[tidx]['llh'] for tidx in tidxs])
  probs = util.softmax(llhs)
  clusters = np.array([len(results.tree_summary[tidx]['populations']) for tidx in tidxs]) - 1
  expected_clusters = np.sum(probs * clusters)
  return expected_clusters

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('tree_summary',
    help='JSON-formatted tree summaries')
  parser.add_argument('mutation_list',
    help='JSON-formatted list of mutations')
  parser.add_argument('mutation_assignment',
    help='JSON-formatted list of SSMs and CNVs assigned to each subclone')
  args = parser.parse_args()

  results = ResultLoader(args.tree_summary, args.mutation_list, args.mutation_assignment)
  C = count_clusters(results)
  print(C)

if __name__ == '__main__':
  main()
