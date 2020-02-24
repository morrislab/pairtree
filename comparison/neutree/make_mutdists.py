import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import mutdist
import mutstat
import neutree

def remove_garbage_from_baseline(baseline, garbage):
  garbage = set(garbage)
  garb_idxs = []
  new_vids = []
  for idx, vid in enumerate(baseline.vids):
    if vid in garbage:
      garb_idxs.append(idx)
    else:
      new_vids.append(vid)

  new_phis = np.delete(baseline.stats, garb_idxs, axis=0)
  new_baseline = mutstat.Mutstat(stats=new_phis, vids=new_vids, assays=baseline.assays)
  return new_baseline

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--impute-garbage', action='store_true')
  parser.add_argument('neutree_fn')
  parser.add_argument('baseline_fn')
  parser.add_argument('mutdist_fn')
  args = parser.parse_args()

  ntree = neutree.load(args.neutree_fn)
  baseline = mutstat.load(args.baseline_fn)
  if args.impute_garbage:
    baseline = remove_garbage_from_baseline(baseline, ntree.garbage)

  mdist = mutdist.calc_mutdist(ntree.phis, ntree.logscores, ntree.clusterings, baseline, ntree.counts)
  if args.impute_garbage:
    mdist = mutstat.impute_garbage(mdist, ntree.garbage, lambda vid: 1)
  mutstat.write(mdist, args.mutdist_fn)

if __name__ == '__main__':
  main()
