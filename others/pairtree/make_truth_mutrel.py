import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import common
import inputparser

def make_mutrel(structure, clusters):
  adjm = common.convert_adjlist_to_adjmatrix(structure)
  mutrel = common.make_mutrel_tensor_from_cluster_adj(adjm, clusters)
  return mutrel

def load_structure(params):
  return {int(P): C for P, C in params['structure'].items()}

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('pairtree_params_fn')
  parser.add_argument('mutrel_fn')
  args = parser.parse_args()

  params = inputparser.load_params(args.pairtree_params_fn)
  structure = load_structure(params)
  mutrel = make_mutrel(structure, params['clusters'])
  np.savez_compressed(args.mutrel_fn, soft_mutrel=mutrel)

if __name__ == '__main__':
  main()
