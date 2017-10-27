from __future__ import print_function
from pwgsresults.result_loader import ResultLoader as PwgsResults
import re
import argparse
import numpy as np
import itertools

def lpdist(A, B, p=2):
  A, B = np.array(A), np.array(B)
  dist = np.sum(np.abs(A - B)**p)**(1./p)
  assert dist >= 0
  return dist

def compare(s1idx, s2idx, tidx, treesumm):
  pops = treesumm[tidx]['populations']
  s1cp = []
  s2cp = []

  for pop in pops.values():
    s1cp.append(pop['cellular_prevalence'][s1idx])
    s2cp.append(pop['cellular_prevalence'][s2idx])
  return lpdist(s1cp, s2cp, p=1)

def calc_concordance(samples, treesumm):
  #Use handbuilt tree.
  tidx = 0

  allsampidxs = [idx for idx, S in enumerate(samples) if re.search(r'^Diagnosis Xeno', S)]
  distlist = []
  for A, B in itertools.combinations(allsampidxs, 2):
    dist = compare(A, B, tidx, treesumm)
    distlist.append(( (A, B), dist  ))
  distlist.sort(key = lambda E: E[-1])
  dists = {}
  for pos, ((A, B), dist) in enumerate(distlist):
    dists[(A, B)] = dists[(B, A)] = (dist, (pos + 1.) / len(distlist))

  for sampidx, sampname in enumerate(samples):
    if not re.search(r'^Diagnosis Xeno \d+$', sampname):
      continue
    matches = [('%s %s' % (sampname, suffix)) in samples for suffix in ('CNS', 'Spleen')]
    for suffix in ('CNS', 'Spleen'):
      other = '%s %s' % (sampname, suffix)
      if other not in samples:
        print(sampname, other, 'absent', sep='\t')
        continue
      otheridx = samples.index(other)
      # Since every CP is in [0, 1], we can get normalized distance by dividing
      # by numer of samples.
      dist, rank = dists[(sampidx, otheridx)]
      print(
        sampname,
        other,
        dist,
        dist / len(samples),
        rank,
        sep='\t'
      )

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('treesumm_fn')
  parser.add_argument('mutlist_fn')
  args = parser.parse_args()

  loader = PwgsResults(args.treesumm_fn, args.mutlist_fn, None)
  samples = loader.params['samples']
  treesumm = loader.tree_summary

  calc_concordance(samples, treesumm)

main()
