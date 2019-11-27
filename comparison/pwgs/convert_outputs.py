import argparse
import sys
import os
import numpy as np

sys.path += [
  os.path.join(os.path.dirname(__file__), '..'),
  os.path.join(os.path.dirname(__file__), '..', '..', 'lib'),
  os.path.expanduser('~/.apps/phylowgs')
]
from pwgsresults.result_loader import ResultLoader
import common
import mutrel
import inputparser
import evalutil
import mutphi
import util
import neutree

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

def convert_results(results, base_clusters, garbage, use_supervars):
  tidxs = sorted(results.tree_summary.keys())
  llhs = np.array([results.tree_summary[tidx]['llh'] for tidx in tidxs])
  adjms = []
  phis = []
  clusterings = []

  for idx, (tidx, mutass) in enumerate(results.load_all_mut_assignments()):
    # Ensure trees are provided in sorted order.
    assert idx == tidx
    if use_supervars:
      full_mutass = replace_supervar_with_variants(base_clusters, mutass)
    else:
      full_mutass = {cid: mutass[cid]['ssms'] for cid in mutass.keys()}
    adjm = common.convert_adjlist_to_adjmatrix(results.tree_summary[tidx]['structure'])
    phi = extract_phi(results.tree_summary[tidx]['populations'])
    pwgs_clusters = [full_mutass[cidx] if cidx in full_mutass else [] for cidx in range(len(adjm))]
    adjms.append(adjm)
    phis.append(phi)
    clusterings.append(pwgs_clusters)

  structs = [util.convert_adjmatrix_to_parents(A) for A in adjms]
  counts = np.ones(len(structs))
  ntree = neutree.Neutree(
    structs = structs,
    phis = phis,
    counts = counts,
    logscores = llhs,
    clusterings = clusterings,
    garbage = [[] for idx in range(len(structs))],
  )
  return ntree

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--use-supervars', action='store_true')
  # This takes Pairtree rather than PWGS inputs, which seems a little weird,
  # but it's okay -- the PWGS inputs are the supervariants, but we need to know
  # which variants correspond to each cluster in the original Pairtree inputs.
  parser.add_argument('tree_summary',
    help='JSON-formatted tree summaries')
  parser.add_argument('mutation_list',
    help='JSON-formatted list of mutations')
  parser.add_argument('mutation_assignment',
    help='JSON-formatted list of SSMs and CNVs assigned to each subclone')
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('neutree_fn')
  args = parser.parse_args()

  results = ResultLoader(args.tree_summary, args.mutation_list, args.mutation_assignment)
  if args.use_supervars:
    params = inputparser.load_params(args.pairtree_params_fn)
    base_clusters = params['clusters']
    garbage = params['garbage']
  else:
    base_clusters = None
    garbage = []

  ntree = convert_results(results, base_clusters, garbage, args.use_supervars)
  neutree.save(ntree, args.neutree_fn)

if __name__ == '__main__':
  main()
