import numpy as np
import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import inputparser
import common
import evalutil

def load_mutphis(phi_args):
  mutphis = {}
  for phi_arg in phi_args:
    name, phi_path = phi_arg.split('=', 1)
    assert name not in mutphis
    if os.path.exists(phi_path):
      mutphi = np.load(phi_path)
      mutphis[name] = evalutil.Mutphi(phi=mutphi['phi'], vids=mutphi['vids'])
    else:
      mutphis[name] = None
  return mutphis

def compare(mutphis):
  names = sorted(mutphis.keys())
  scores = {}
  truth = mutphis['truth']

  for name in names:
    mutphi = mutphis[name]
    if mutphi is None:
      scores[name] = -1
      continue
    assert np.array_equal(mutphi.vids, truth.vids)
    assert mutphi.phi.shape == truth.phi.shape
    score = np.mean(1 - (np.abs(mutphi.phi - truth.phi)))
    if name == 'truth':
      assert np.allclose(1, score)
    scores[name] = score

  names.remove('truth')
  print(*names, sep=',')
  print(*[scores[name] for name in names], sep=',')

def check_complete(mutphis, clustered):
  for name, mutphi in mutphis.items():
    if mutphi is None:
      continue
    assert set(mutphi.vids) == clustered

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
  check_complete(mutphis, clustered)
  compare(mutphis)

if __name__ == '__main__':
  main()
