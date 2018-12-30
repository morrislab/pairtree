import numpy as np
import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import mutrel

def load_mutrels(mutrel_args):
  mutrels = {}
  for mutrel_arg in mutrel_args:
    mutrel_name, mutrel_path = mutrel_arg.split('=', 1)
    assert mutrel_name not in mutrels
    mutrels[mutrel_name] = np.load(mutrel_path)['mutrel']
    assert np.all(0 <= mutrels[mutrel_name]) and np.all(mutrels[mutrel_name] <= 1)
  return mutrels

def compare(mutrels):
  assert 'truth' in mutrels
  M, _, num_models = mutrels['truth'].shape
  assert mutrels['truth'].shape == (M, M, num_models)
  correct = np.argmax(mutrels['truth'], axis=2)

  names = sorted(mutrels.keys())
  scores = {}

  for name in names:
    mrel = mutrels[name]
    assert mrel.shape == (M, M, num_models)
    mutrel.check_posterior_sanity(mrel)
    score = np.nan * np.zeros((M, M))
    for I in range(M):
      for J in range(M):
        # I don't know how to use `correct` in this format to index into
        # `mutrels[name]`, hence the `I` and `J` for loops.
        score[I,J] = mrel[I,J,correct[I,J]]
    assert not np.any(np.isnan(score))
    if name == 'truth':
      assert np.allclose(1, score)
    scores[name] = np.mean(score)
    
  names.remove('truth')
  print(*names, sep=',')
  print(*[scores[name] for name in names], sep=',')

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('mutrels', nargs='+')
  args = parser.parse_args()

  mutrels = load_mutrels(args.mutrels)
  compare(mutrels)

if __name__ == '__main__':
  main()
