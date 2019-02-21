import numpy as np
import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import mutrel
import inputparser
from common import Models

def load_mutrels(mutrel_args):
  mutrels = {}
  for mutrel_arg in mutrel_args:
    mutrel_name, mutrel_path = mutrel_arg.split('=', 1)
    assert mutrel_name not in mutrels
    if os.path.exists(mutrel_path):
      mrel = np.load(mutrel_path)
      mutrels[mutrel_name] = mutrel.Mutrel(vids=mrel['vids'], rels=mrel['rels'])
    else:
      mutrels[mutrel_name] = None
  return mutrels

def discard_garbage(mutrels, clustered, garbage):
  for name in list(mutrels.keys()):
    if mutrels[name] is None:
      continue
    vids = mutrels[name].vids
    assert set(vids) == clustered | garbage

    gidxs = [idx for idx, vid in enumerate(vids) if vid in garbage]
    for garbrels in (mutrels[name].rels[gidxs,:,Models.garbage], mutrels[name].rels[:,gidxs,Models.garbage].T):
      # Garbage mutations should have posterior that coclusters with
      # themselves, but that renders them garbage to every other mutation
      # (including other garbage)
      G, M = garbrels.shape
      # `expected` shape, given `G` garbage variants and `M` total variants: `GxM`
      expected = np.ones((G, M))
      expected[np.arange(G),gidxs] = 0
      assert np.allclose(expected, garbrels), '%s garbage relations are wrong' % name
    mutrels[name] = mutrel.remove_variants(mutrels[name], gidxs)
    assert set(mutrels[name].vids) == clustered

def compare(mutrels):
  assert 'truth' in mutrels
  M, _, num_models = mutrels['truth'].rels.shape
  assert mutrels['truth'].rels.shape == (M, M, num_models)
  correct = np.argmax(mutrels['truth'].rels, axis=2)

  names = sorted(mutrels.keys())
  scores = {}

  for name in names:
    mrel = mutrels[name]
    if mrel is None:
      scores[name] = -1
      continue

    assert mrel.rels.shape == (M, M, num_models)
    assert np.array_equal(mrel.vids, mutrels['truth'].vids)
    mutrel.check_posterior_sanity(mrel.rels)

    scores[name] = np.mean(1 - np.abs(mrel.rels - mutrels['truth'].rels))
    if name == 'truth':
      assert np.isclose(1, scores[name])
    
  names.remove('truth')
  print(*names, sep=',')
  print(*[scores[name] for name in names], sep=',')

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--discard-garbage', dest='discard_garbage', action='store_true')
  parser.add_argument('--params', dest='paramsfn', required=True)
  parser.add_argument('mutrels', nargs='+')
  args = parser.parse_args()

  params = inputparser.load_params(args.paramsfn)
  garbage = set(params['garbage'])
  clustered = set([vid for C in params['clusters'] for vid in C])

  mutrels = load_mutrels(args.mutrels)
  if args.discard_garbage:
    discard_garbage(mutrels, clustered, garbage)
  compare(mutrels)

if __name__ == '__main__':
  main()
