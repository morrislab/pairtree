import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import evalutil
import neutree

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('neutree_fn')
  parser.add_argument('mutrel_fn')
  args = parser.parse_args()

  ntree = neutree.load(args.neutree_fn)
  best_idx = np.argmax(ntree.logscores)
  structs = [ntree.structs[best_idx]]
  clusterings = [ntree.clusterings[best_idx]]
  logscores = [0]

  mrel = evalutil.make_mutrel_from_trees_and_unique_clusterings(structs, logscores, clusterings)
  # Note we don't add garbage to this structure. We assume comparisons won't use garbage.
  evalutil.save_sorted_mutrel(mrel, args.mutrel_fn)

if __name__ == '__main__':
  main()
