import numpy as np
import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import inputparser
import common
import mutphi

MISSING = -1

def load_mutphis(mphi_args):
  mutphis = {}
  for mphi_arg in mphi_args:
    name, mphi_path = mphi_arg.split('=', 1)
    assert name not in mutphis
    if os.path.exists(mphi_path):
      mphi = mutphi.load_mutphi(mphi_path)
      # Sometimes, an atrociously bad phi for a mutation will result in it
      # getting a mutphi score of -inf. If the corresponding tree is assigned
      # an LLH of > -inf, then it will "pollute" the mutphis for *every* treein
      # a run for that mutation/sample, making all of them -inf. That of course
      # renders the overall mutphi score for that run -inf. In such instances,
      # consider the run as having failed.
      if np.any(np.isinf(mphi.logprobs)):
        print('%.5f of mphis in %s are inf' % (
          np.sum(np.isinf(mphi.logprobs)) / mphi.logprobs.size,
          mphi_path,
        ), file=sys.stderr)
        mphi = None
      mutphis[name] = mphi
    else:
      mutphis[name] = None
  return mutphis

def score(logprobs):
  assert np.all(logprobs <= 0)
  score = -np.sum(logprobs)
  score /= logprobs.size
  # Convert to bits.
  score /= np.log(2)
  return score

def compare(mutphis):
  names = list(mutphis.keys())
  scores = {}
  present = [N for N in names if mutphis[N] is not None]
  first_present = present[0]
  vids = mutphis[first_present].vids
  assays = mutphis[first_present].assays

  for name in names:
    mphi = mutphis[name]
    if mphi is None:
      scores[name] = MISSING
      continue
    assert np.array_equal(mphi.vids, vids)
    assert np.array_equal(mphi.assays, assays)
    scores[name] = score(mphi.logprobs)

  print(*names, sep=',')
  print(*[scores[name] for name in names], sep=',')

def check_incomplete(mutphis, clustered):
  names = list(mutphis.keys())
  for name in names:
    if mutphis[name] is None:
      continue
    vids = set(mutphis[name].vids)
    assert vids.issubset(clustered)
    if vids != clustered:
      missing = clustered - vids
      msg = '%s lacks fraction=%s variants (%s)' % (name, len(missing) / len(clustered), missing)
      mutphis[name] = None
      raise Exception(msg)

def remove_garbage(mutphis, garbage):
  # Remove garbage if present. Some mutphis have it (e.g., PWGS run on all
  # variants, without benefit of handbuilt clusters) while others don't.
  garbage = set(garbage)
  revised = {}

  for name, mphi in mutphis.items():
    if mphi is None:
      revised[name] = None
      continue
    garb_idxs = [idx for idx, vid in enumerate(mphi.vids) if vid in garbage]
    new_vids = [vid for idx, vid in enumerate(mphi.vids) if idx not in set(garb_idxs)]
    new_logprobs = np.delete(mphi.logprobs, garb_idxs, axis=0)
    revised[name] = mutphi.Mutphi(logprobs=new_logprobs, vids=new_vids, assays=mphi.assays)

  return revised

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--params', dest='paramsfn', required=True)
  parser.add_argument('mutphis', nargs='+')
  args = parser.parse_args()

  params = inputparser.load_params(args.paramsfn)
  clustered = set([vid for C in params['clusters'] for vid in C])

  mutphis = load_mutphis(args.mutphis)
  mutphis = remove_garbage(mutphis, params['garbage'])
  check_incomplete(mutphis, clustered)
  compare(mutphis)

if __name__ == '__main__':
  main()
