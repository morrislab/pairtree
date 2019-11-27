import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import common
import inputparser
import neutree
import util

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

def load_list(listfn):
  with open(listfn) as F:
    vals = [L.strip() for L in F.readlines()]
  return vals

def load_mutass(tid, resultfn, vidfn):
  # Assignments are variants to tree nodes.
  # Returns `mutass` dict, where `mutass[i]` is nodex index of variant `i`

  # Definitions:
  # `M` variants
  # `C` clusters
  # `N` nodes
  # (we should have `N = C`)

  # `varass`: `M`-length vector, where `varass[i]` is node index of variant `i`
  # `vids`: `M`-length vector, where `vids[i]` is VID of variant `i` (CITUP doesn't use VIDs, so this is my own mapping that their code never sees)
  varass = pd.read_hdf(resultfn, f'trees/{tid}/variant_assignment').to_numpy()
  vids = load_list(vidfn)
  assert len(varass) == len(vids)

  # `mutass`: dict, where `mutass[n]` is VIDs associated with node `n`
  mutass = defaultdict(list)
  for vid, nid in zip(vids, varass):
    mutass[nid].append(vid)
  return dict(mutass)

def load_mutass_using_clusters(tid, resultfn, vidfn, clustfn):
  # Assignments are variants to tree nodes.
  # Returns `mutass` dict, where `mutass[i]` is nodex index of variant `i`

  # Definitions:
  # `M` variants
  # `C` clusters
  # `N` nodes
  # (we should have `N = C`)

  # `clustass`: `C`-length vector, where `clustass[i]` is cluster index of variant `i`
  # `vids`: `M`-length vector, where `vids[i]` is VID of variant `i` (CITUP doesn't use VIDs, so this is my own mapping that their code never sees)
  clustass = [int(C) for C in load_list(clustfn)]
  vids = load_list(vidfn)
  assert len(clustass) == len(vids)

  clusters = defaultdict(list)
  for C, V in zip(clustass, vids):
    clusters[C].append(V)
  # If subsequent code accesses a nonexistent index in `clusters`, I want it to
  # throw an exception, so convert to regular dict.
  clusters = dict(clusters)

  # Assignments are clusters to tree nodes.
  # Returns `mutass`, which is dictionary giving tree nodes to VIDs.
  # `clust2node`: `C`-length vector, where `clust2node[c]` is nodex index of cluster `c`
  # `mutass`: dict, where `mutass[n]` is VIDs associated with node `n`
  clust2node = pd.read_hdf(resultfn, f'trees/{tid}/cluster_assignment').to_numpy()
  mutass = defaultdict(list)
  for C, N in enumerate(clust2node):
    # This mapping isn't guaranteed to be one-to-one -- sometimes we can have
    # multiple clusters assigned to a single node.
    mutass[N] += clusters[C]
  return dict(mutass)

def _extract_clustering(mutass, K):
  clustering = [list(mutass[C]) if C in mutass else [] for C in range(K)]
  return clustering

def _check_results(results):
  def _extract_members(clusters):
    return set([vid for C in clusters for vid in C])

  K, S = results[0]['phi'].shape
  members = _extract_members(results[0]['clusters'])

  for R in results:
    assert R['adjm'].shape == (K, K)
    assert R['phi'].shape == (K, S)
    assert len(R['clusters']) == K
    M = _extract_members(R['clusters'])
    assert M == members

def load_results(resultfn, vidfn, citup_clusters_fn, clusters, use_supervars):
  # NB:
  # * `clusters` are the Pairtree clusters, which we need only if we're using supervars
  # * `citup_clusters_fn` lists the assignment of variants to clusters, which
  #   we need only if CITUP used fixed clusters, which it does only when we're
  #   using the QIP rather than the iterative approximation mode in CITUP.

  tree_ids = pd.read_hdf(resultfn, 'results/tree_id').to_numpy()
  llh = pd.read_hdf(resultfn, 'results/likelihood').to_numpy()
  llh *= -1
  results = []

  for tid, tree_llh in zip(tree_ids, llh):
    adjm = load_tree(tid, resultfn)
    phi = load_phi(tid, resultfn)
    if citup_clusters_fn is None:
      mutass = load_mutass(tid, resultfn, vidfn)
    else:
      mutass = load_mutass_using_clusters(tid, resultfn, vidfn, citup_clusters_fn)

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
      'clusters': _extract_clustering(mutass, K),
      'llh': tree_llh,
    })

  assert len(results) > 0
  _check_results(results)
  return results

def write_neutree(results, neutree_fn):
  adjms = [R['adjm'] for R in results]
  structs = [util.convert_adjmatrix_to_parents(A) for A in adjms]
  llhs = [R['llh'] for R in results]
  clusterings = [R['clusters'] for R in results]
  phis = [R['phi'] for R in results]
  N = len(structs)
  counts = np.ones(N)
  ntree = neutree.Neutree(
    structs = structs,
    phis = phis,
    counts = counts,
    logscores = llhs,
    clusterings = clusterings,
    garbage = [[] for idx in range(N)],
  )
  neutree.save(ntree, neutree_fn)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--use-supervars', action='store_true')
  parser.add_argument('--citup-clusters')
  parser.add_argument('citup_result_fn')
  parser.add_argument('citup_vid_fn')
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('neutree_fn')
  args = parser.parse_args()

  params = inputparser.load_params(args.pairtree_params_fn)
  clusters = params['clusters']
  results = load_results(args.citup_result_fn, args.citup_vid_fn, args.citup_clusters, clusters, args.use_supervars)
  write_neutree(results, args.neutree_fn)

if __name__ == '__main__':
  main()
