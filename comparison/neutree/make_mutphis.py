import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import resultserializer
import inputparser
import mutphi
import mutstat
import neutree

def _impute(vid, variants):
  V = variants[vid]
  # This is 1/total_reads
  return -np.log(V['total_reads'].astype(np.float))

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--impute-garbage', action='store_true')
  parser.add_argument('neutree_fn')
  parser.add_argument('pairtree_ssm_fn')
  parser.add_argument('mutphi_fn')
  args = parser.parse_args()

  ntree = neutree.load(args.neutree_fn)
  mphi = mutphi.calc_mutphi(ntree.phis, ntree.logscores, ntree.clusterings, args.pairtree_ssm_fn, ntree.counts)
  variants = inputparser.load_ssms(args.pairtree_ssm_fn)
  if args.impute_garbage:
    mphi = mutstat.impute_garbage(mphi, ntree.garbage, lambda vid: _impute(vid, variants))
  mutstat.write(mphi, args.mutphi_fn)

if __name__ == '__main__':
  main()
