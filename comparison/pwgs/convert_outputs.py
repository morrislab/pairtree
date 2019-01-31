import argparse
import sys
import os
import numpy as np

sys.path += [
  os.path.join(os.path.dirname(__file__), '..'),
  os.path.join(os.path.dirname(__file__), '..', '..'),
  os.path.expanduser('~/.apps/phylowgs')
]
from pwgsresults.result_loader import ResultLoader
import common
import mutrel
import inputparser
import evalutil

def replace_supervar_with_variants(clusters, mutass):
  expanded = {}
  for nodeid in mutass.keys():
    assert len(mutass[nodeid]['cnvs']) == 0
    cidxs = [int(cid[1:]) for cid in mutass[nodeid]['ssms']]
    expanded[nodeid] = [vid for cidx in cidxs for vid in clusters[cidx]]
  return expanded

def extract_phi(pops):
  phi = [pops[C]['cellular_prevalence'] for C in sorted(pops.keys())]
  return np.array(phi)

def calc_mutrel_and_mutphi(results, base_clusters, tree_weights):
  tidxs = sorted(results.tree_summary.keys())
  llhs = np.array([results.tree_summary[tidx]['llh'] for tidx in tidxs])
  adjms = []
  phis = []
  clusterings = []

  for idx, (tidx, mutass) in enumerate(results.load_all_mut_assignments()):
    # Ensure trees are provided in sorted order.
    assert idx == tidx
    full_mutass = replace_supervar_with_variants(base_clusters, mutass)
    adjm = common.convert_adjlist_to_adjmatrix(results.tree_summary[tidx]['structure'])
    phi = extract_phi(results.tree_summary[tidx]['populations'])
    pwgs_clusters = [full_mutass[cidx] if cidx in full_mutass else [] for cidx in range(len(adjm))]
    adjms.append(adjm)
    phis.append(phi)
    clusterings.append(pwgs_clusters)

  return (
    evalutil.calc_mutrel_from_trees(adjms, llhs, clusterings, tree_weights),
    evalutil.calc_mutphi(phis, llhs, clusterings, tree_weights),
  )

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--weight-trees-by', choices=('llh', 'uniform'), default='uniform')
  # This takes Pairtree rather than PWGS inputs, which seems a little weird,
  # but it's okay -- the PWGS inputs are the supervariants, ut we need to know
  # hich variants correspond to each cluster in the original Pairtree inputs.
  parser.add_argument('--trees-mutrel')
  parser.add_argument('--phi', dest='mutphifn')
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('tree_summary',
    help='JSON-formatted tree summaries')
  parser.add_argument('mutation_list',
    help='JSON-formatted list of mutations')
  parser.add_argument('mutation_assignment',
    help='JSON-formatted list of SSMs and CNVs assigned to each subclone')
  args = parser.parse_args()

  results = ResultLoader(args.tree_summary, args.mutation_list, args.mutation_assignment)
  params = inputparser.load_params(args.pairtree_params_fn)
  base_clusters = params['clusters']
  mrel, mphi = calc_mutrel_and_mutphi(results, base_clusters, args.weight_trees_by)
  mrel = evalutil.add_garbage(mrel, params['garbage'])

  if args.trees_mutrel is not None:
    evalutil.save_sorted_mutrel(mrel, args.trees_mutrel)
  if args.mutphifn is not None:
    evalutil.save_mutphi(mphi, args.mutphifn)

if __name__ == '__main__':
  main()
