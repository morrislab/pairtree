import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import inputparser
import mutstat

def score(stats, P):
  assert np.all(stats >= 0) and np.all(stats <= 1)
  total = (np.sum(np.abs(stats)**P) / stats.size)**(1/P)
  return total

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
  mutdists = mutstat.load_mutstats(args.mutdists)
  mutdists = mutstat.remove_garbage(mutdists, params['garbage'])
  mutstat.check_incomplete(mutdists, params['clusters'])

  names, scores = mutstat.score_mutstats(mutdists, _score = lambda stats: score(stats, args.p))
  print(*names, sep=',')
  print(*[scores[name] for name in names], sep=',')

if __name__ == '__main__':
  main()
