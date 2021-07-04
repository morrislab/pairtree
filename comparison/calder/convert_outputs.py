import argparse
import numpy as np

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import neutree
import inputparser

def _load_mat(F):
  name = None
  while True:
    line = next(F).strip()
    if len(line) > 0:
      name = line
      mat = []
      break

  col_labels = next(F).strip().split(',')
  while True:
    fields = next(F).strip().split(',')
    row_label = fields[0]
    row_vals = np.array([float(V) for V in fields[1:]])
    mat.append(row_vals)
  mat = np.array(mat)

  return (name, mat)

def _load_mats(matfn):
  mats = {}
  with open(matfn) as F:
    while True:
      name, mat = _load_mat(F)
      mats[name] = mat
  return mats

def _load_struct(structfn):
  pass


def convert(params_fn, trees_fn, neutree_fn):
  N = len(structs)
  params = inputparser.load_params(params_fn)
  ntree = neutree.Neutree(
    structs = structs,
    phis = phis,
    counts = np.ones(N),
    logscores = llhs,
    clusterings = clusterings,
    garbage = params['garbage'],
  )
  neutree.save(ntree, neutree_fn)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('calder_mats_fn')
  parser.add_argument('calder_trees_fn')
  parser.add_argument('neutree_fn')
  args = parser.parse_args()

  convert(args.sampid, args.pairtree_params_fn, args.trees_fn, args.neutree_fn)

if __name__ == '__main__':
  main()
