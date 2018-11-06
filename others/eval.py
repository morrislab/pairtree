import numpy as np
import argparse

def load_mutrels(mutrel_args):
  mutrels = {}
  for mutrel_arg in mutrel_args:
    mutrel_name, mutrel_path = mutrel_arg.split('=', 1)
    assert mutrel_name not in mutrels
    mutrels[mutrel_name] = np.load(mutrel_path)['soft_mutrel']
    assert np.all(0 <= mutrels[mutrel_name]) and np.all(mutrels[mutrel_name] <= 1)
  return mutrels

def compare(mutrels):
  assert 'truth' in mutrels
  M, _, num_models = mutrels['truth'].shape
  assert mutrels['truth'].shape == (M, M, num_models)

  #correct = np.argmax(mutrels['truth'], axis=2)
  #assert correct.shape == (M, M)

  names = sorted(mutrels.keys())
  scores = {}

  for name in names:
    mutrel = mutrels[name]
    assert mutrel.shape == (M, M, num_models)
    score = np.nan * np.zeros((M, M))
    for I in range(M):
      for J in range(M):
        correct = np.argmax(mutrels['truth'][I,J])
        score[I,J] = mutrel[I,J,correct]
    assert not np.any(np.isnan(score))
    scores[name] = np.mean(score)
    
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
