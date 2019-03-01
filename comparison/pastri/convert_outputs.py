import argparse
from io import StringIO
import numpy as np

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import inputparser
import evalutil

class UnexpectedShapeError(Exception):
  pass

def load_matrices(matfn, dtype=np.float, check_shapes=True):
  mats = {}
  matname = None
  shapes = {}

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
        # Ensure shapeensions not specified multiple times.
        assert matname not in shapes
        shapes[matname] = tuple([int(T.strip()) for T in line[1:-1].split(',')])
      else:
        mats[matname] += line + '\n'

  for matname in mats:
    mats[matname] = np.loadtxt(StringIO(mats[matname]), dtype=dtype)

    # If the matrix is `Rx1`, the parsing code will actually have loaded it as
    # an `R`-length vector, rather than a `Rx1` matrix. Correct this.
    if mats[matname].ndim == 1:
      mats[matname] = mats[matname][:,np.newaxis]

    if check_shapes:
      assert matname in shapes
      if shapes[matname] != mats[matname].shape:
        raise UnexpectedShapeError('expected=%s actual=%s' % (shapes[matname], mats[matname].shape))
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
  # structure is provided as a Kx2 matrix of edges. Matrix shape is not
  # provided.
  mats = load_matrices(trees_fn, dtype=np.int, check_shapes=False)
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
  phis = []
  clusterings = []

  for tree in prelim_trees.values():
    tidx = tree['rank']

    freqfn = _makefn(tidx, 'F')
    # Sometimes, if PASTRI's get_F_and_C.py crashed, this file won't exist,
    # even though the corresponding tree in the .trees file will have non-zero
    # LH.
    if not os.path.exists(freqfn):
      continue
    try:
      freqs = load_matrices(freqfn, np.float)
    except UnexpectedShapeError:
      # Sometimes the stated shape of the matrix is `Kx1`, but PASTRI gives us
      # `kx1`, where `k < K`. Ignore these cases.
      continue
    assert len(freqs) == 1
    F = next(iter(freqs.values()))

    K = len(F)
    clustfn = _makefn(tidx, 'C')
    C = load_clusters(clustfn, K)
    clusters = convert_pastri_clusters(C, variant_clusters)

    treefn = _makefn(tidx, 'labeled_trees')
    tree_adjms = load_final_trees(treefn)
    assert len(tree_adjms) > 0
    for A in tree_adjms:
      adjms.append(A)
      # There are multiple possible permutations of tree nodes within this
      # structure. All have the same LLH, clustering, and phis.
      llhs.append(tree['llh'])
      phis.append(F)
      clusterings.append(clusters)

  assert len(adjms) == len(llhs) == len(clusterings)
  return (adjms, llhs, phis, clusterings)

def write_results(sampid, params_fn, trees_fn, tree_weights, trees_mutrel_fn, mutphi_fn):
  params = inputparser.load_params(params_fn)
  outdir = os.path.dirname(trees_fn)
  prelim_trees = load_prelim_trees(trees_fn)

  adjms, llhs, phis, clusterings = convert_results(sampid, prelim_trees, params['clusters'], outdir)
  if len(adjms) == 0:
    return

  if trees_mutrel_fn is not None:
    # On occasion, PASTRI won't assign mutations to all clusters -- it can have
    # empty clusters. This will cause some of the assertions in my parsing code
    # to fail, which is fine -- we just won't use its results in those cases.
    mrel = evalutil.calc_mutrel_from_trees(adjms, llhs, clusterings, tree_weights)
    mrel = evalutil.add_garbage(mrel, params['garbage'])
    evalutil.save_sorted_mutrel(mrel, trees_mutrel_fn)

  if mutphi_fn is not None:
    # PASTRI has some phis > 1 (e.g., 1.00079953). Ugh.
    phis = np.minimum(1, phis)
    mphi = evalutil.calc_mutphi(phis, llhs, clusterings, tree_weights)
    evalutil.save_mutphi(mphi, mutphi_fn)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--weight-trees-by', choices=('llh', 'uniform'), default='uniform')
  # This takes Pairtree rather than PWGS inputs, which seems a little weird,
  # but it's okay -- the PWGS inputs are the supervariants, ut we need to know
  # hich variants correspond to each cluster in the original Pairtree inputs.
  parser.add_argument('--trees-mutrel', dest='mutrelfn')
  parser.add_argument('--phi', dest='mutphifn')
  parser.add_argument('sampid')
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('trees_fn')
  args = parser.parse_args()

  if args.mutrelfn is not None:
    write_results(args.sampid, args.pairtree_params_fn, args.trees_fn, args.weight_trees_by, args.mutrelfn, args.mutphifn)

if __name__ == '__main__':
  main()
