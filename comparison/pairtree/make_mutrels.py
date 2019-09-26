import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import resultserializer
import evalutil
import pairtree_util
import util
import numpy as np

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--use-subset', type=int)
  parser.add_argument('--clustrel-mutrel')
  parser.add_argument('--trees-mutrel')
  parser.add_argument('pairtree_results_fn')
  args = parser.parse_args()

  results = resultserializer.load(args.pairtree_results_fn)
  clusters = [[]] + list(results['clusters'])
  garbage = list(results['garbage'])
  all_vids = set([V for C in results['clusters'] for V in C] + garbage)

  if args.use_subset is not None:
    pairtree_util.choose_subset(results, args.use_subset)

  if args.trees_mutrel is not None:
    tree_mutrel = evalutil.make_mutrel_from_trees_and_single_clustering(results['struct'], results['llh'], results['count'], clusters)
    tree_mutrel = evalutil.add_garbage(tree_mutrel, garbage)
    assert set(tree_mutrel.vids) == all_vids
    evalutil.save_sorted_mutrel(tree_mutrel, args.trees_mutrel)

  if args.clustrel_mutrel is not None:
    clustrel_mutrel = evalutil.make_mutrel_from_clustrel(results['clustrel_posterior'], clusters)
    clustrel_mutrel = evalutil.add_garbage(clustrel_mutrel, garbage)
    assert set(clustrel_mutrel.vids) == all_vids
    evalutil.save_sorted_mutrel(clustrel_mutrel, args.clustrel_mutrel)

if __name__ == '__main__':
  main()
