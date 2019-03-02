import numpy as np
import argparse
import scipy.stats

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
      mutphis[name] = evalutil.load_mutphi(phi_path)
    else:
      mutphis[name] = None
  return mutphis

def score_phis(mutphi, variants):
  def _extract_arr(K, vids):
    return np.array([variants[V][K] for V in vids])
  vardata = {K: _extract_arr(K, mutphi.vids) for K in ('var_reads', 'total_reads', 'omega_v')}
  assert len(set([V.shape for V in vardata.values()] + [mutphi.phis[0].shape])) == 1

  log_px = scipy.stats.binom.logpmf(
    vardata['var_reads'][np.newaxis,:,:],
    vardata['total_reads'][np.newaxis,:,:],
    vardata['omega_v'][np.newaxis,:,:] * mutphi.phis,
  )
  assert np.all(log_px <= 0)

  scores = np.sum(log_px, axis=(1,2))
  scores /= mutphi.phis[0].size
  assert np.all(scores < 0)
  # Convert to bits.
  scores /= np.log(2)

  assert np.all(mutphi.weights >= 0)
  assert np.isclose(1, np.sum(mutphi.weights))
  score = np.sum(mutphi.weights * scores)
  return -score

def compare(mutphis, variants):
  names = list(mutphis.keys())
  scores = {}
  vids = mutphis[names[0]].vids
  mutphi_shape = mutphis[names[0]].phis[0].shape

  for name in names:
    mutphi = mutphis[name]
    if mutphi is None:
      scores[name] = -1
      continue
    assert np.array_equal(mutphi.vids, vids)
    assert mutphi.phis[0].shape == mutphi_shape
    scores[name] = score_phis(mutphi, variants)

  print(*names, sep=',')
  print(*[scores[name] for name in names], sep=',')

def check_complete(mutphis, clustered):
  for name, mutphi in mutphis.items():
    if mutphi is None:
      continue
    assert set(mutphi.vids) == clustered

def remove_garbage(mutphis, garbage):
  # Remove garbage if present. Some mutphis have it (e.g., PWGS run on all
  # variants, without benefit of handbuilt clusters) while others don't.
  garbage = set(garbage)
  revised = {}

  for name, mutphi in mutphis.items():
    if mutphi is None:
      revised[name] = None
      continue
    garb_idxs = [idx for idx, vid in enumerate(mutphi.vids) if vid in garbage]
    new_vids = [vid for idx, vid in enumerate(mutphi.vids) if idx not in set(garb_idxs)]
    new_phi = np.delete(mutphi.phis, garb_idxs, axis=1)
    revised[name] = evalutil.Mutphi(phis=new_phi, vids=new_vids, weights=mutphi.weights)

  return revised

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--params', dest='paramsfn', required=True)
  parser.add_argument('--ssms', dest='ssmsfn', required=True)
  parser.add_argument('mutphis', nargs='+')
  args = parser.parse_args()

  params = inputparser.load_params(args.paramsfn)
  variants = inputparser.load_ssms(args.ssmsfn)
  clustered = set([vid for C in params['clusters'] for vid in C])

  mutphis = load_mutphis(args.mutphis)
  mutphis = remove_garbage(mutphis, params['garbage'])
  check_complete(mutphis, clustered)
  compare(mutphis, variants)

if __name__ == '__main__':
  main()
