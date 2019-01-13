import argparse
import sys
import os
import numpy as np

sys.path += [
  os.path.join(os.path.dirname(__file__), '..', '..'),
  os.path.expanduser('~/.apps/phylowgs')
]
from pwgsresults.result_loader import ResultLoader
import common
import mutrel
import clustermaker
import inputparser

def softmax(V):
  V = np.copy(V) - np.max(V)
  return np.exp(V) / np.sum(np.exp(V))

def replace_supervar_with_variants(supervars, clusters, mutass):
  expanded = {}
  for nodeid in mutass.keys():
    assert len(mutass[nodeid]['cnvs']) == 0
    cidxs = [int(cid[1:]) for cid in mutass[nodeid]['ssms']]
    expanded[nodeid] = [vid for cidx in cidxs for vid in clusters[cidx]]
  return expanded

def calc_mutrel(results, base_clusters, supervars, tree_weights):
  tidxs = sorted(results.tree_summary.keys())
  soft_mutrel = None

  llhs = np.array([results.tree_summary[tidx]['llh'] for tidx in tidxs])
  if tree_weights == 'llh':
    weights = softmax(llhs)
  elif tree_weights == 'uniform':
    weights = np.ones(len(llhs)) / len(llhs)
  else:
    raise Exception('Unknown tree_weights=%s' % tree_weights)

  for tidx, mutass in results.load_all_mut_assignments():
    weight = weights[tidx]
    full_mutass = replace_supervar_with_variants(supervars, base_clusters, mutass)
    adjm = common.convert_adjlist_to_adjmatrix(results.tree_summary[tidx]['structure'])
    pwgs_clusters = [full_mutass[cidx] if cidx in full_mutass else [] for cidx in range(len(adjm))]
    mrel = mutrel.make_mutrel_tensor_from_cluster_adj(adjm, pwgs_clusters)

    if soft_mutrel is None:
      soft_mutrel = np.zeros(mrel.rels.shape)
    soft_mutrel += weight * mrel.rels

  return soft_mutrel

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--weight-trees-by', choices=('llh', 'uniform'), default='uniform')
  parser.add_argument('pairtree_ssm_fn')
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('tree_summary',
    help='JSON-formatted tree summaries')
  parser.add_argument('mutation_list',
    help='JSON-formatted list of mutations')
  parser.add_argument('mutation_assignment',
    help='JSON-formatted list of SSMs and CNVs assigned to each subclone')
  parser.add_argument('mutrel_fn')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.pairtree_ssm_fn)
  params = inputparser.load_params(args.pairtree_params_fn)
  supervars = clustermaker.make_cluster_supervars(params['clusters'], variants)

  base_clusters = params['clusters']
  results = ResultLoader(args.tree_summary, args.mutation_list, args.mutation_assignment)
  mrel = calc_mutrel(results, base_clusters, supervars, args.weight_trees_by)

  np.savez_compressed(args.mutrel_fn, mutrel=mrel)

if __name__ == '__main__':
  main()
