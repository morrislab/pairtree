import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import inputparser
import resultserializer
import mutphi
import common

def sort_mutphi(mphi):
  sorted_vids = common.sort_vids(mphi.vids)
  mapping = [mphi.vids.index(V) for V in sorted_vids]
  assert sorted_vids == [mphi.vids[idx] for idx in mapping]
  sorted_logprobs = np.array([mphi.logprobs[idx] for idx in mapping])
  return mutphi.Mutphi(
    vids = sorted_vids,
    assays = mphi.assays,
    logprobs = sorted_logprobs,
  )

def impute(ssmfn, params, mphi):
  clustered = set([V for C in params['clusters'] for V in C])
  mphi_vids = set(mphi.vids)
  missing = list(clustered - mphi_vids)
  if len(missing) == 0:
    sys.exit()

  variants = inputparser.load_ssms(ssmfn)
  missing_reads = np.array([variants[V]['total_reads'] for V in missing]).astype(np.float)
  assert np.all(missing_reads >= 1)
  # Assign uniform probability based on total read count.
  missing_logprobs = np.log(1 / missing_reads)

  combined = mutphi.Mutphi(
    vids = list(mphi.vids) + missing,
    assays = mphi.assays,
    logprobs = np.vstack((mphi.logprobs, missing_logprobs)),
  )
  return combined

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
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('mutphi_fn')
  args = parser.parse_args()

  params = inputparser.load_params(args.params_fn)
  orig_mphi = mutphi.load_mutphi(args.mutphi_fn)
  mphi = impute(args.ssm_fn, params, orig_mphi)
  mphi = sort_mutphi(mphi)
  mutphi.write_mutphi(mphi, args.mutphi_fn)

  old, new = score(orig_mphi.logprobs), score(mphi.logprobs)
  #print('score_cmp', old, new, new - old, (new - old) > 0)

if __name__ == '__main__':
  main()
