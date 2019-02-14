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
      mutphi = np.load(phi_path)
      mutphis[name] = evalutil.Mutphi(phi=mutphi['phi'], vids=mutphi['vids'])
    else:
      mutphis[name] = None
  return mutphis

def score_phis(mutphi, variants):
  def _extract_arr(K, vids):
    return np.array([variants[V][K] for V in vids])
  vardata = {K: _extract_arr(K, mutphi.vids) for K in ('var_reads', 'total_reads', 'omega_v')}
  assert len(set([V.shape for V in vardata.values()] + [mutphi.phi.shape])) == 1

  log_px = scipy.stats.binom.logpmf(
    vardata['var_reads'],
    vardata['total_reads'],
    vardata['omega_v'] * mutphi.phi,
  )
  assert np.all(log_px <= 0)

  score = np.sum(log_px)
  score /= mutphi.phi.size
  # Convert to bits.
  score /= np.log(2)
  return -score

def compare(mutphis, variants):
  names = list(mutphis.keys())
  scores = {}
  vids = mutphis[names[0]].vids
  mutphi_shape = mutphis[names[0]].phi.shape

  for name in names:
    mutphi = mutphis[name]
    if mutphi is None:
      scores[name] = -1
      continue
    assert np.array_equal(mutphi.vids, vids)
    assert mutphi.phi.shape == mutphi_shape
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
    new_phi = np.delete(mutphi.phi, garb_idxs, axis=0)
    revised[name] = evalutil.Mutphi(phi=new_phi, vids=new_vids)

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
