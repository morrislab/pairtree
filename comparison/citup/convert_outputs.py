import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import common
import inputparser

def replace_supervar_with_variants(clusters, mutass):
  expanded = {}
  for cid in mutass.keys():
    expanded[cid] = []
    for svid in mutass[cid]:
      assert svid.startswith('S')
      svidx = int(svid[1:])
      expanded[cid] += clusters[svidx]
  return expanded

def load_tree(tid, resultfn):
  adjl = pd.read_hdf(resultfn, f'trees/{tid}/adjacency_list').to_numpy()
  parents, children = adjl[:,0], adjl[:,1]
  K = np.max(children) + 1
  adjm = np.eye(K)
  adjm[parents,children] = 1
  common.ensure_valid_tree(adjm)
  return adjm

def load_phi(tid, resultfn):
  phi = pd.read_hdf(resultfn, f'trees/{tid}/clade_freq').to_numpy().T
  return phi

def load_vids(vidfn):
  with open(vidfn) as F:
    vids = [vid.strip() for vid in F.readlines()]
  return vids

def load_mutass(tid, resultfn, vidfn):
  assvec = pd.read_hdf(resultfn, f'trees/{tid}/variant_assignment').to_numpy()
  vids = load_vids(vidfn)
  assert len(assvec) == len(vids)

  mutass = defaultdict(list)
  for vid, cid in zip(vids, assvec):
    mutass[cid].append(vid)
  return dict(mutass)

def load_results(resultfn, vidfn, clusters, use_supervars):
  tree_ids = pd.read_hdf(resultfn, 'results/tree_id').to_numpy()
  llh = pd.read_hdf(resultfn, 'results/likelihood').to_numpy()
  results = []

  for tid, tree_llh in zip(tree_ids, llh):
    adjm = load_tree(tid, resultfn)
    phi = load_phi(tid, resultfn)
    mutass = load_mutass(tid, resultfn, vidfn)

    K, S = phi.shape
    assert adjm.shape == (K, K)
    if False and len(mutass) != K - 1:
      # This condition seems common, so suppress the print.
      print('Nodes', set(range(1, K)) - set(mutass.keys()), 'have no assigned variants')

    if use_supervars:
      mutass = replace_supervar_with_variants(clusters, mutass)

    results.append({
      'adjm': adjm,
      'phi': phi,
      'mutass': mutass,
      'llh': tree_llh,
    })

  return results

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--weight-trees-by', choices=('llh', 'uniform'), default='uniform')
  parser.add_argument('--mutrel', dest='mutrel_fn')
  parser.add_argument('--mutphi', dest='mutphi_fn')
  parser.add_argument('--use-supervars', action='store_true')
  parser.add_argument('citup_result_fn')
  parser.add_argument('citup_vid_fn')
  parser.add_argument('pairtree_ssm_fn')
  parser.add_argument('pairtree_params_fn')
  args = parser.parse_args()

  params = inputparser.load_params(args.pairtree_params_fn)
  clusters = params['clusters']
  results = load_results(args.citup_result_fn, args.citup_vid_fn, clusters, args.use_supervars)
  from IPython import embed
  embed()
  import sys
  sys.exit()

if __name__ == '__main__':
  main()
