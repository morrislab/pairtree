import argparse

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import evalutil
import pastri_util
import inputparser

def write_mutrels(sampid, params_fn, trees_fn, weight_trees_by, trees_mutrel_fn):
  adjms, llhs, phis, clusterings = pastri_util.load_results(sampid, params_fn, trees_fn)
  if len(adjms) == 0:
    return

  params = inputparser.load_params(params_fn)
  # On occasion, PASTRI won't assign mutations to all clusters -- it can have
  # empty clusters. This will cause some of the assertions in my parsing code
  # to fail, which is fine -- we just won't use its results in those cases.
  mrel = evalutil.calc_mutrel_from_trees(adjms, llhs, clusterings, weight_trees_by)
  mrel = evalutil.add_garbage(mrel, params['garbage'])
  evalutil.save_sorted_mutrel(mrel, trees_mutrel_fn)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--weight-trees-by', choices=('llh', 'uniform'), default='uniform')
  # This takes Pairtree rather than PWGS inputs, which seems a little weird,
  # but it's okay -- the PWGS inputs are the supervariants, ut we need to know
  # hich variants correspond to each cluster in the original Pairtree inputs.
  parser.add_argument('sampid')
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('trees_fn')
  parser.add_argument('mutrel_fn')
  args = parser.parse_args()

  write_mutrels(args.sampid, args.pairtree_params_fn, args.trees_fn, args.weight_trees_by, args.mutrel_fn)

if __name__ == '__main__':
  main()
