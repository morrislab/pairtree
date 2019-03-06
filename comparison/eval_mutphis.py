import numpy as np
import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import inputparser
import common
import evalutil

def load_mutphis(mphi_args):
  mutphis = {}
  for phi_arg in phi_args:
    name, phi_path = phi_arg.split('=', 1)
    assert name not in mutphis
    if os.path.exists(phi_path):
      mutphis[name] = mutphi.load_mutphi(mphi_path)
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
  vids = mutphis[names[0]].vids
  assays = mutphis[names[0]].assays

  for name in names:
    mphi = mutphis[name]
    if mphi is None:
      scores[name] = -1
      continue
    assert np.array_equal(mphi.vids, vids)
    assert np.array_equal(mphi.assays, assays)
    scores[name] = score(mphi.logprobs)

  print(*names, sep=',')
  print(*[scores[name] for name in names], sep=',')

def check_complete(mutphis, clustered):
  for name, mphi in mutphis.items():
    if mphi is None:
      continue
    assert set(mphi.vids) == clustered

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
  check_complete(mutphis, clustered)
  compare(mutphis)

if __name__ == '__main__':
  main()
