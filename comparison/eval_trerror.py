# Uses tree error defined in CALDER paper. This assumes a single solution tree
# and a single true tree. It just computes the fraction of all mutation pairs
# in the candidate tree that have the same relation as in the true tree.

import numpy as np
import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import mutrel
from common import Models

MISSING = -1

def load_mutrels(mutrel_args):
  mutrels = {}
  for mutrel_arg in mutrel_args:
    mutrel_name, mutrel_path = mutrel_arg.split('=', 1)
    assert mutrel_name not in mutrels, '%s is duplicate' % mutrel_name
    if os.path.exists(mutrel_path):
      mrel = np.load(mutrel_path)
      mutrels[mutrel_name] = mutrel.Mutrel(vids=mrel['vids'], rels=mrel['rels'])
    else:
      mutrels[mutrel_name] = None
  return mutrels

def _ensure_hard_relations(mrel):
  assert np.logical_or(mrel.rels == 0., mrel.rels == 1.)

def _score_treeror(candidate, truth):
  M, _, num_models = truth.rels.shape
  assert M >= 2
  for R in (candidate, truth):
    _ensure_hard_relations(R)
    assert R.rels.shape == (M,M,num_models)

  _extract_rel = lambda R: np.tril(np.argmax(R.rels, axis=2), -1)
  correct = _extract_rel(candidate)
  soln = _extract_rel(truth)
  num_pairs = 0.5*M*(M-1)
  prop_correct = sum(correct == soln) / num_pairs
  assert 0 <= prop_correct <= 1
  trerror = 1 - prop_correct
  return trerror

def compare(mutrels):
  assert 'truth' in mutrels
  M, _, num_models = mutrels['truth'].rels.shape
  assert mutrels['truth'].rels.shape == (M, M, num_models)

  names = sorted(mutrels.keys())
  scores = {}

  for name in names:
    mrel = mutrels[name]
    if mrel is None:
      scores[name] = MISSING
      continue
    assert mrel.rels.shape == (M, M, num_models)
    assert np.array_equal(mrel.vids, mutrels['truth'].vids)
    mutrel.check_posterior_sanity(mrel.rels)

    scores[name] = _score_trerror(mrel.rels, mutrels['truth'].rels)
    if name == 'truth':
      assert np.isclose(0, scores[name])
    
  print(*names, sep=',')
  print(*[scores[name] for name in names], sep=',')

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('mutrels', nargs='+')

  mutrels = load_mutrels(args.mutrels)
  compare(mutrels)

if __name__ == '__main__':
  main()
