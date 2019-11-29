import numpy as np
import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import inputparser
import mutstat

MISSING = -1

def _check_infs(mphis):
  for name, mphi in mphis.items():
    if mphi is None:
      continue
    # Sometimes, an atrociously bad phi for a mutation will result in it
    # getting a mutphi score of -inf. If the corresponding tree is assigned
    # an LLH of > -inf, then it will "pollute" the mutphis for *every* treein
    # a run for that mutation/sample, making all of them -inf. That of course
    # renders the overall mutphi score for that run -inf. In such instances,
    # consider the run as having failed.
    if np.any(np.isinf(mphi.stats)):
      print('%.5f of logprobs in %s are inf' % (
        np.sum(np.isinf(mphi.logprobs)) / mphi.logprobs.size,
        name,
      ), file=sys.stderr)
      mphis[name] = None

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
  parser.add_argument('--params', dest='paramsfn', required=True)
  parser.add_argument('mutphis', nargs='+')
  args = parser.parse_args()

  params = inputparser.load_params(args.paramsfn)
  mutphis = mutstat.load_mutstats(args.mutphis, check_inf=False)
  # We do our own NaN check rather than relying on the one in
  # `mutstat.load_mutstats`, since we want to handle NaN logprobs that we
  # sometimes get from PASTRI.
  _check_infs(mutphis)
  mutphis = mutstat.remove_garbage(mutphis, params['garbage'])
  mutstat.check_incomplete(mutphis, params['clusters'])

  names, scores = mutstat.score_mutstats(mutphis, _score=score)
  print(*names, sep=',')
  print(*[scores[name] for name in names], sep=',')

if __name__ == '__main__':
  main()
