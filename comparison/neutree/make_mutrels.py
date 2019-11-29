import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import evalutil
import neutree

def are_same_clusterings(clusterings):
  # Tuple conversion probably isn't necessary, but it makes everything
  # hashable, so it's probably a good idea.
  _convert_to_tuples = lambda C: tuple([tuple(cluster) for cluster in C])
  first_C = _convert_to_tuples(clusterings[0])
  for C in clusterings[1:]:
    if _convert_to_tuples(C) != first_C:
      return False
  return True

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('neutree_fn')
  parser.add_argument('mutrel_fn')
  args = parser.parse_args()

  ntree = neutree.load(args.neutree_fn)
  clusterings = ntree.clusterings

  if are_same_clusterings(clusterings):
    mrel = evalutil.make_mutrel_from_trees_and_single_clustering(ntree.structs, ntree.logscores, ntree.counts, clusterings[0])
  else:
    mrel = evalutil.make_mutrel_from_trees_and_unique_clusterings(ntree.structs, ntree.logscores, clusterings)
  mrel = evalutil.add_garbage(mrel, ntree.garbage)
  evalutil.save_sorted_mutrel(mrel, args.mutrel_fn)

if __name__ == '__main__':
  main()
