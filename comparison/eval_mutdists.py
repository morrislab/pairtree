import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import inputparser
import mutdist

MISSING = -1

# Much of the code here is shared with the `mutphi` module. Several of the
# functions (at least when I developed `mutphi`) are verbatim copies. I should
# consider merging these eventually.

def load_mutdists(mdist_args):
  mutdists = {}
  for mdist_arg in mdist_args:
    name, mdist_path = mdist_arg.split('=', 1)
    assert name not in mutdists
    if os.path.exists(mdist_path):
      mutdists[name] = mutdist.load_mutdist(mdist_path)
      assert not np.any(np.isinf(mutdists[name].dists)) and not np.any(np.isnan(mutdists[name].dists))
    else:
      mutdists[name] = None
  return mutdists

def score(dists, P):
  assert np.all(dists >= 0) and np.all(dists <= 1)
  total = np.sum(dists**P)**(1/P)
  total /= dists.size
  return total

def remove_garbage(mutdists, garbage):
  # Remove garbage if present. Some mutdists have it (e.g., PWGS run on all
  # variants, without benefit of handbuilt clusters) while others don't.
  garbage = set(garbage)
  revised = {}

  for name, mdist in mutdists.items():
    if mdist is None:
      revised[name] = None
      continue
    garb_idxs = [idx for idx, vid in enumerate(mdist.vids) if vid in garbage]
    new_vids = [vid for idx, vid in enumerate(mdist.vids) if idx not in set(garb_idxs)]
    new_dists = np.delete(mdist.dists, garb_idxs, axis=0)
    revised[name] = mutdist.Mutdist(dists=new_dists, vids=new_vids, assays=mdist.assays)

  return revised

def check_incomplete(mutdists, clustered):
  names = list(mutdists.keys())
  for name in names:
    if mutdists[name] is None:
      continue
    vids = set(mutdists[name].vids)
    assert vids.issubset(clustered)
    if vids != clustered:
      missing = clustered - vids
      msg = '%s lacks fraction=%s variants (%s)' % (name, len(missing) / len(clustered), missing)
      mutdists[name] = None
      raise Exception(msg)

def compare(mutdists, P):
  names = list(mutdists.keys())
  scores = {N: MISSING for N in names}
  present = [N for N in names if mutdists[N] is not None]
  if len(present) == 0:
    return (names, scores)
  first_present = present[0]
  vids = mutdists[first_present].vids
  assays = mutdists[first_present].assays

  for name in present:
    mdist = mutdists[name]
    assert np.array_equal(mdist.vids, vids)
    assert np.array_equal(mdist.assays, assays)
    scores[name] = score(mdist.dists, P)
  return (names, scores)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('-p', dest='p', type=float, required=True)
  parser.add_argument('--params', dest='paramsfn', required=True)
  parser.add_argument('mutdists', nargs='+')
  args = parser.parse_args()

  params = inputparser.load_params(args.paramsfn)
  clustered = set([vid for C in params['clusters'] for vid in C])

  mutdists = load_mutdists(args.mutdists)
  mutdists = remove_garbage(mutdists, params['garbage'])
  check_incomplete(mutdists, clustered)

  names, scores = compare(mutdists, args.p)
  print(*names, sep=',')
  print(*[scores[name] for name in names], sep=',')

if __name__ == '__main__':
  main()
