import os
import sys
for name in ('comparison', 'lib'):
  sys.path.append(os.path.join(os.path.dirname(__file__), '..', name))

import argparse
import mutphi
import resultserializer
import numpy as np

def score(logprobs):
  assert np.all(logprobs <= 0)
  score = -np.sum(logprobs)
  score /= logprobs.size
  # Convert to bits.
  score /= np.log(2)
  return score

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--only-best', action='store_true')
  parser.add_argument('pairtree_ssm_fn')
  parser.add_argument('results_fn')
  args = parser.parse_args()

  results = resultserializer.Results(args.results_fn)
  phi = results.get('phi')
  clusters = [[]] + results.get('clusters')
  llh = results.get('llh')
  counts = results.get('count')
  clusterings = [clusters for _ in range(len(llh))]

  if args.only_best:
    phi = [phi[0]]
    clusterings = [clusterings[0]]
    llh = [llh[0]]
    counts = [1]

  mphi = mutphi.calc_mutphi(phi, llh, clusterings, args.pairtree_ssm_fn, counts)
  print(score(mphi.stats))

if __name__ == '__main__':
  main()
