import numpy as np
import argparse

import common
import handbuilt
import json

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

def load_sampnames(paramsfn):
  with open(paramsfn) as F:
    params = json.load(F)
  return params['samples']

def filter_samples(mat, desired_samps, all_samps):
  _, cols = mat.shape
  assert len(all_samps) == cols
  assert set(desired_samps).issubset(set(all_samps))
  sidxs = [all_samps.index(S) for S in desired_samps]
  return mat[:,sidxs]

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--tree-type', dest='tree_type', choices = ('handbuilt.xeno', 'handbuilt.patient'), default='handbuilt.patient')
  parser.add_argument('sampid')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('handbuilt_json_fn')
  parser.add_argument('pastri_inputs_fn')
  parser.add_argument('pastri_proposal_fn')
  args = parser.parse_args()

  hbjson = handbuilt._load_handbuilt(args.handbuilt_json_fn, args.tree_type)
  variants = common.parse_ssms(args.sampid, args.ssm_fn)
  clusters = handbuilt._load_clusters(hbjson, variants)
  supervars, svid2svidx, svidx2svid = common.make_cluster_supervars(clusters, variants)
  sampnames = load_sampnames(args.params_fn)

  matrices = {
    'var_reads': extract_matrix(variants, 'var_reads'),
    'total_reads': extract_matrix(variants, 'total_reads'),
    'alpha': extract_matrix(supervars, 'var_reads'),
    'beta': extract_matrix(supervars, 'total_reads'),
  }

  if args.tree_type == 'handbuilt.patient':
    desired_samps = ('D', 'R1')
    for name, mat in matrices.items():
      matrices[name] = filter_samples(mat, desired_samps, sampnames)

  #N, S = var_reads.shape
  #C = len(supervars)
  #alpha = 3*np.ones((C, S))
  #beta  = 3*np.ones((C, S))

  C_max = 10
  matrices['alpha'] = matrices['alpha'][:C_max,]
  matrices['beta']  = matrices['beta'][:C_max,]

  write_matrices(('A', matrices['var_reads']), ('D', matrices['total_reads']), outfn = args.pastri_inputs_fn)
  write_matrices(('A', matrices['alpha']),     ('B', matrices['beta']),        outfn = args.pastri_proposal_fn)

  from IPython import embed
  embed()

if __name__ == '__main__':
  main()
