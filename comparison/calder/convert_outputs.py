import argparse
import numpy as np
import re

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import neutree
import inputparser
import common

def _load_mat(F):
  name = None
  while True:
    line = next(F).strip()
    if len(line) > 0:
      name = line
      mat = []
      break

  col_labels = next(F).strip().split(',')
  row_labels = []
  while True:
    try:
      line = next(F).strip()
    except StopIteration:
      break
    if len(line) == 0:
      break
    fields = line.split(',')
    row_labels.append(fields[0])
    row_vals = np.array([float(V) for V in fields[1:]])
    mat.append(row_vals)

  assert name is not None and len(col_labels) > 0 and len(row_labels) > 0
  mat = np.array(mat)
  # Return transpose.
  return (name, col_labels, row_labels, mat.T)

def _load_mats(matfn):
  mats = {}
  row_labels = {}
  col_labels = {}

  with open(matfn) as F:
    while True:
      try:
        name, row_lab, col_lab, mat  = _load_mat(F)
      except StopIteration:
        break
      mats[name] = mat
      row_labels[name] = row_lab
      col_labels[name] = col_lab
  return (mats, row_labels, col_labels)

def _load_struct(svids, structfn):
  # NB: Calder seems to always assume monoprimary tree, so it doesn't represent
  # normal root in tree structure (implying you can't have polyprimary).
  patterns = {
    'label':  r'^(v\d+) \[label="([^"]+)"\];$',
    'edge':   r'^([a-z0-9]+) -> ([a-z0-9]+);$',
    'eof':    r'^}$',
    'header': r'^digraph sim_tree1 {$',
  }
  node_map = {}
  edges = []

  with open(structfn) as F:
    while True:
      line = next(F).strip()
      for line_type, pat in patterns.items():
        R = re.search(pat, line)
        if R:
          break
      else:
        raise Exception('No match for line: %s' % line)

      if line_type == 'header':
        continue
      if line_type == 'label':
        node_map[R.group(1)] = R.group(2)
      elif line_type == 'edge':
        edges.append((R.group(1), R.group(2)))
      elif line_type == 'eof':
        break
      else:
        raise Exception('Unknown line_type %s' % line_type)

  parents = {node_map[B]: node_map[A] for (A,B) in edges}
  struct = np.array([svids.index(parents[child]) for child in parents.keys()])
  assert len(struct) == len(svids) - 1
  return struct

def convert(params_fn, calder_mats_fn, calder_trees_fn, neutree_fn):
  params = inputparser.load_params(params_fn)

  mats, row_labels, col_labels = _load_mats(calder_mats_fn)
  assert row_labels['Fhat'][0] == 'samples'
  svids = row_labels['Fhat'][1:]
  assert svids == common.sort_vids(svids)

  struct = _load_struct(svids, calder_trees_fn)
  ntree = neutree.Neutree(
    structs = [struct],
    phis = [mats['Fhat']],
    counts = np.array([1]),
    logscores = np.array([0.]),
    clusterings = [params['clusters']],
    garbage = params['garbage'],
  )
  neutree.save(ntree, neutree_fn)

def _is_full_tree(stdout_fn):
  # No way to tell without checking CALDER STDOUT whether it built a tree
  # incorporating all mutations.
  with open(stdout_fn) as F:
    S = F.read()
  M = re.search(r'Maximal tree contains (\d+) out of (\d+) mutations/clusters\.', S)
  if not M:
    return False
  A, B = M.groups()
  if A != B:
    assert int(A) < int(B)
    return False
  return True

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('calder_mats_fn')
  parser.add_argument('calder_trees_fn')
  parser.add_argument('calder_stdout_fn')
  parser.add_argument('neutree_fn')
  args = parser.parse_args()

  #if not _is_full_tree(args.calder_stdout_fn):
  #  raise Exception('not full tree')
  convert(args.pairtree_params_fn, args.calder_mats_fn, args.calder_trees_fn, args.neutree_fn)

if __name__ == '__main__':
  main()
