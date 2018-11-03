import argparse
from io import StringIO
import numpy as np

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
        pass
      else:
        mats[matname] += line + '\n'

  for matname in mats:
    mats[matname] = np.loadtxt(StringIO(mats[matname]), dtype=dtype)
  return mats

def find_parents_from_ancestral(anc):
  # Basic idea: determine how many ancestors each node `i` has by examining the
  # nonzero entries in column `i`. If I do this in order of number of
  # ancestors, the parent of `i` is guaranteed to have already been placed in
  # the tree. When looking at all possible parents of `i` from amongst `i`'s
  # ancestors, choose the parent that has the greatest depth in the tree, as
  # this must be the direct parent.
  anc = np.copy(anc)
  np.fill_diagonal(anc, 0)
  K = len(anc)
  num_anc = np.sum(anc, axis=0)

  parents = np.nan * np.ones(K)
  depth = {}

  for I in np.flatnonzero(num_anc == 0):
    parents[I] = -1
    depth[I] = 0
  assert np.sum(parents == -1) == 1

  for A in range(1, K):
    if not np.any(np.isnan(parents)):
      break
    for I in np.flatnonzero(num_anc == A):
      possible_parents = np.flatnonzero(anc[:,I])
      assert len(possible_parents) > 0
      pp_depth = [(P, depth[P]) for P in possible_parents]
      # Sort by maximum depth.
      pp_depth.sort(key = lambda pair: -pair[1])
      # Ensure all possible parents are at unique depths.
      assert len(set([D for P, D in pp_depth])) == len(possible_parents)
      parent = pp_depth[0][0]
      parents[I] = parent
      depth[I] = A

  assert not np.any(np.isnan(parents))
  return parents
    
def load_trees(trees_fn):
  mats = load_matrices(trees_fn, dtype=np.int)
  trees = {}
  for K in mats:
    fields = K.split(':')
    assert len(fields) == 3
    tree_id = int(fields[1])
    assert tree_id not in trees
    trees[tree_id] = {
      'rank': int(fields[0]),
      'llh': float(fields[2]),
      'structure': find_parents_from_ancestral(mats[K]),
    }
  return trees

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('trees_fn')
  args = parser.parse_args()
  trees = load_trees(args.trees_fn)
  print(trees)

if __name__ == '__main__':
  main()
