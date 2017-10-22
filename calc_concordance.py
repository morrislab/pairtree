from __future__ import print_function
from pwgsresults.result_loader import ResultLoader as PwgsResults
import re
import argparse
import numpy as np

def euclid(A, B):
  A, B = np.array(A), np.array(B)
  return np.sqrt(np.sum((A - B)**2))

def compare(s1idx, s2idx, tidx, treesumm):
  pops = treesumm[tidx]['populations']
  s1cp = []
  s2cp = []

  for pop in pops.values():
    s1cp.append(pop['cellular_prevalence'][s1idx])
    s2cp.append(pop['cellular_prevalence'][s2idx])
  return euclid(s1cp, s2cp)

def calc_concordance(samples, treesumm):
  #Use handbuilt tree.
  tidx = 0

  for sampidx, sampname in enumerate(samples):
    if not re.search(r'^Diagnosis Xeno \d+$', sampname):
      continue
    matches = [('%s %s' % (sampname, suffix)) in samples for suffix in ('CNS', 'Spleen')]
    for suffix in ('CNS', 'Spleen'):
      other = '%s %s' % (sampname, suffix)
      if other not in samples:
        continue
      otheridx = samples.index(other)
      dist = compare(sampidx, otheridx, tidx, treesumm)
      # Since every CP is in [0, 1], we can get normalized distance by dividing
      # by numer of samples.
      print(sampname, other, dist, dist / len(samples), sep='\t')

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
