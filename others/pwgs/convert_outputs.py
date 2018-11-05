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
import inputparser

def softmax(V):
  V = np.copy(V) - np.max(V)
  return np.exp(V) / np.sum(np.exp(V))

def replace_supervar_with_variants(supervars, clusters, mutass):
  expanded = {}
  for cidx in mutass.keys():
    assert len(mutass[cidx]['cnvs']) == 0
    cidxs = [int(cid[1:]) for cid in mutass[cidx]['ssms']]
    expanded[cidx] = [ssmid for cidx in cidxs for ssmid in clusters[supervars['C%s' % cidx]['cluster']]]
  return expanded

def calc_soft_mutrel(results, base_clusters, supervars):
  tidxs = sorted(results.tree_summary.keys())
  llhs = np.array([results.tree_summary[tidx]['llh'] for tidx in tidxs])
  weights = softmax(llhs)
  soft_mutrel = None

  for tidx, weight in zip(tidxs, weights):
    mutass = results.load_mut_assignments(tidx)
    full_mutass = replace_supervar_with_variants(supervars, base_clusters, mutass)
    adjm = common.convert_adjlist_to_adjmatrix(results.tree_summary[tidx]['structure'])
    pwgs_clusters = [full_mutass[cidx] if cidx in full_mutass else [] for cidx in range(len(adjm))]
    mutrel = common.make_mutrel_tensor_from_cluster_adj(adjm, pwgs_clusters)

    if soft_mutrel is None:
      soft_mutrel = np.zeros(mutrel.shape)
    soft_mutrel += weight * mutrel

  return soft_mutrel

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('pairtree_ssm_fn')
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('tree_summary',
    help='JSON-formatted tree summaries')
  parser.add_argument('mutation_list',
    help='JSON-formatted list of mutations')
  parser.add_argument('mutation_assignment',
    help='JSON-formatted list of SSMs and CNVs assigned to each subclone')
  parser.add_argument('soft_mutrel_fn')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.pairtree_ssm_fn)
  params = inputparser.load_params(args.pairtree_params_fn)
  supervars = common.make_cluster_supervars(params['clusters'], variants)

  base_clusters = params['clusters']
  results = ResultLoader(args.tree_summary, args.mutation_list, args.mutation_assignment)
  soft_mutrel = calc_soft_mutrel(results, base_clusters, supervars)

  np.savez_compressed(args.soft_mutrel_fn, soft_mutrel=soft_mutrel)

if __name__ == '__main__':
  main()
