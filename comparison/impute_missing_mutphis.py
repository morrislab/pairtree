import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
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
  _, S = mphi.logprobs.shape
  
  clustered = set([V for C in params['clusters'] for V in C])
  mphi_vids = set(mphi.vids)
  missing = list(clustered - mphi_vids)
  if len(missing) == 0:
    sys.exit()

  # Incorrectly declaring a mutation as garbage is tantamount to not putting it
  # in the tree, meaning it should get a phi of zero.
  missing_phi = np.zeros((1, S))
  missing_mutphi = mutphi.calc_mutphi(
    [missing_phi],
    llhs = [0],
    clusterings = [[missing]],
    weight_trees_by = 'llh',
    ssmsfn = ssmfn,
  )

  combined = mutphi.Mutphi(
    vids = list(mphi.vids) + list(missing),
    assays = mphi.assays,
    logprobs = np.vstack((mphi.logprobs, missing_mutphi.logprobs)),
  )
  return combined

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

if __name__ == '__main__':
  main()
