import argparse
from io import StringIO
import numpy as np

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import inputparser
import evalutil

def load_matrices(matfn, dtype=np.float):
  mats = {}
  matname = None

  with open(matfn) as outf:
    for line in outf:
      line = line.strip()
      if line == '':
        continue
      elif line.startswith('>'):
        matname = line[1:].strip()
        assert matname not in mats
        mats[matname] = ''
      elif line.startswith('(') and line.endswith(')'):
        # This is just the shape of the matrix. Ignore.
        # Note that "tree labelling" files don't have this header. We don't
        # require it as a consequence.
        pass
      else:
        mats[matname] += line + '\n'

  for matname in mats:
    mats[matname] = np.loadtxt(StringIO(mats[matname]), dtype=dtype)
  return mats

def load_prelim_trees(trees_fn):
  # These are initial trees, representing trees with distinct LHs. Structure is
  # an ancestral matrix, not an adjacency matrix. Multiple permutations of the
  # nodes may be possible (all with the same LH, as their frequencies stay the
  # same), but those aren't listed here. We load the trees here to figure out
  # which sampled frequencies have a non-zero LH.
  mats = load_matrices(trees_fn, dtype=np.int)
  trees = {}

  for K in mats:
    fields = K.split(':')
    assert len(fields) == 3

    tree_id = int(fields[1])
    tree = {
      'rank': int(fields[0]),
      'llh': float(fields[2]),
      'structure': mats[K],
    }
    if tree['llh'] > -np.inf:
      assert tree_id not in trees
      trees[tree_id] = tree

  return trees

def load_final_trees(trees_fn):
  # These are final trees, representing all possible permutations of nodes with
  # the sampled frequencies that are consistent with any tree structure. Tree
  # structure is provided as a Kx2 matrix of edges.
  mats = load_matrices(trees_fn, dtype=np.int)
  adjms = []

  for edges in mats.values():
    parents, children = edges[:,0], edges[:,1]
    max_child = np.max(children)
    assert 0 in parents and 0 not in children

    adjm = np.eye(max_child + 1)
    adjm[parents,children] = 1
    adjms.append(adjm)

  return adjms

def load_clusters(clustfn, K):
  # K: number of clusters
  # Clusters aren't guaranteed to be contiguous -- even if you have `K` nodes
  # in the tree, some of those can have zero mutations assigned. Hence, I force
  # a list at every possible cluster index.
  clusters = {cidx: [] for cidx in range(K)}

  with open(clustfn) as F:
    for line in F:
      tokens = [int(T) for T in line.strip().split()]
      parent, children = tokens[0], tokens[1:]
      clusters[parent] += children

  clusters = [clusters[C] for C in sorted(clusters.keys())]
  return clusters

def convert_pastri_clusters(pastri_clusters, variant_clusters):
  converted = []
  for pcluster in pastri_clusters:
    converted.append([])
    for vcidx in pcluster:
      converted[-1] += variant_clusters[vcidx]
  return converted

def convert_results(sampid, prelim_trees, variant_clusters, outdir):
  _makefn = lambda tidx, ext: os.path.join(outdir, '%s.%s.%s' % (sampid, tidx, ext))
  llhs = []
  adjms = []
  clusterings = []

  for tree in prelim_trees.values():
    tidx = tree['rank']

    freqfn = _makefn(tidx, 'F')
    # Sometimes, if PASTRI's get_F_and_C.py crashed, this file won't exist,
    # even though the corresponding tree in the .trees file will have non-zero
    # LH.
    if not os.path.exists(freqfn):
      continue
    freqs = load_matrices(freqfn, np.float)
    assert len(freqs) == 1
    F = next(iter(freqs.values()))

    K = len(F)
    clustfn = _makefn(tidx, 'C')
    C = load_clusters(clustfn, K)
    clusters = convert_pastri_clusters(C, variant_clusters)

    treefn = _makefn(tidx, 'labeled_trees')
    tree_adjms = load_final_trees(treefn)
    for A in tree_adjms:
      adjms.append(A)
      # There are multiple possible permutations of tree nodes within this
      # structure. All have the same LLH and clustering.
      llhs.append(tree['llh'])
      clusterings.append(clusters)

  assert len(adjms) == len(llhs) == len(clusterings)
  return (adjms, llhs, clusterings)

def write_mutrel(sampid, params_fn, trees_fn, tree_weights, trees_mutrel_fn):
  params = inputparser.load_params(params_fn)
  outdir = os.path.dirname(trees_fn)
  prelim_trees = load_prelim_trees(trees_fn)

  adjms, llhs, clusterings = convert_results(sampid, prelim_trees, params['clusters'], outdir)
  if len(adjms) == 0:
    return
  mrel = evalutil.calc_mutrel_from_trees(adjms, llhs, clusterings, tree_weights)
  evalutil.save_sorted_mutrel(mrel, trees_mutrel_fn)

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
  parser.add_argument('sampid')
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('trees_fn')
  args = parser.parse_args()

  if args.trees_mutrel is not None:
    write_mutrel(args.sampid, args.pairtree_params_fn, args.trees_fn, args.weight_trees_by, args.trees_mutrel)

if __name__ == '__main__':
  main()
