import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import common
import inputparser

def write_snvs(variants, sampnames, garbage, snv_fn, normal_vaf=0.0):
  sampnames = ['Normal'] + sampnames

  with open(snv_fn, 'w') as F:
    print('#chr', 'position', 'description', *sampnames, sep='\t', file=F)
    vids = common.sort_vids(variants.keys())

    idx = 1
    for vid in vids:
      if vid in garbage:
        continue
      vaf = (variants[vid]['var_reads'] / variants[vid]['total_reads']).tolist()
      vaf = [normal_vaf] + vaf
      print('1', idx, vid, *vaf, sep='\t', file=F)
      idx += 1

def extract_mat(variants, key):
  mat = np.array([V[key] for V in variants])
  return mat

def write_clusters(variants, clusters, cluster_fn, normal_vaf=0.0):
  rows = []
  for cluster in clusters:
    cvars = [variants[V] for V in cluster]
    var_reads = np.sum(extract_mat(cvars, 'var_reads'), axis=0)
    total_reads = np.sum(extract_mat(cvars, 'total_reads'), axis=0)
    cvaf = (var_reads / total_reads).tolist()
    cvaf = [normal_vaf] + cvaf

    sampmask = '0' + (len(cvaf) - 1)*'1'
    snv_idxs = [str(int(V[1:]) + 1) for V in common.sort_vids(cluster)]

    rows.append([sampmask] + cvaf + [','.join(snv_idxs)])

  with open(cluster_fn, 'w') as F:
    for row in rows:
      print(*row, sep='\t', file=F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--uniform-proposal', action='store_true')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('lichee_snv_fn')
  parser.add_argument('lichee_cluster_fn')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.ssm_fn)
  params = inputparser.load_params(args.params_fn)
  sampnames = params['samples']
  clusters = params['clusters']
  garbage = set(params['garbage'])

  write_snvs(variants, sampnames, garbage, args.lichee_snv_fn)
  write_clusters(variants, clusters, args.lichee_cluster_fn)

if __name__ == '__main__':
  main()
