import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from common import Models
import resultserializer
import mutphi
import pairtree_util

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--use-subset', type=int)
  parser.add_argument('--weight-trees-by', choices=('llh', 'uniform'), default='uniform')
  parser.add_argument('--ssms', dest='ssms_fn', required=True)
  parser.add_argument('pairtree_results_fn')
  parser.add_argument('mutphi_fn')
  args = parser.parse_args()

  results = resultserializer.load(args.pairtree_results_fn)
  clusters = [[]] + list(results['clusters'])
  if args.use_subset is not None:
    pairtree_util.choose_subset(results, args.use_subset)

  clusterings = [clusters for idx in range(len(results['adjm']))]
  mphi = mutphi.calc_mutphi(results['phi'], results['llh'], clusterings, args.weight_trees_by, args.ssms_fn)
  mutphi.write_mutphi(mphi, args.mutphi_fn)

if __name__ == '__main__':
  main()
