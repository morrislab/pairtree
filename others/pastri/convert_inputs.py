import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import argparse

import common
import inputparser

def extract_matrix(variants, key):
  return np.array([variants[K][key] for K in sorted(variants.keys(), key = lambda vid: int(vid[1:]))])

def write_matrices(*matrices, outfn):
  with open(outfn, 'w') as F:
    for name, matrix in matrices:
      print('> %s' % name, file=F)
      print(matrix.shape, file=F)
      for row in matrix:
        print(*row, sep='\t', file=F)
      print('', file=F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('pastri_allele_counts_fn')
  parser.add_argument('pastri_proposal_fn')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.ssm_fn)
  params = inputparser.load_params(args.params_fn)
  sampnames = params['samples']
  clusters = params['clusters']
  supervars = common.make_cluster_supervars(clusters, variants)

  matrices = {
    'var_reads': extract_matrix(variants, 'var_reads'),
    'total_reads': extract_matrix(variants, 'total_reads'),
    'alpha': extract_matrix(supervars, 'var_reads'),
    'beta': extract_matrix(supervars, 'total_reads'),
  }

  C_max = 10
  matrices['alpha'] = matrices['alpha'][:C_max,]
  matrices['beta']  = matrices['beta'][:C_max,]

  write_matrices(('A', matrices['var_reads']), ('D', matrices['total_reads']), outfn = args.pastri_allele_counts_fn)
  write_matrices(('A', matrices['alpha']),     ('B', matrices['beta']),        outfn = args.pastri_proposal_fn)

if __name__ == '__main__':
  main()
