import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import common
import inputparser
import clustermaker

def load_nongarb_vids(variants, garbage):
  vids = set(variants.keys())
  nongarb_vids = common.sort_vids(vids - set(garbage))
  return nongarb_vids

def write_snvs(variants, garbage, snv_fn, vid_fn):
  vids = load_nongarb_vids(variants, garbage)

  with open(snv_fn, 'w') as F:
    for vid in vids:
      ccf = (variants[vid]['var_reads'] / variants[vid]['total_reads'] / variants[vid]['omega_v']).tolist()
      print(*ccf, sep='\t', file=F)

  with open(vid_fn, 'w') as F:
    print(*vids, sep='\n', file=F)

def write_clusters(variants, garbage, clusters, cluster_fn):
  vids = load_nongarb_vids(variants, garbage)
  cluster_map = {vid: cidx for cidx, C in enumerate(clusters) for vid in C}
  assignments = [cluster_map[vid] for vid in vids]
  with open(cluster_fn, 'w') as F:
    print(*assignments, sep='\n', file=F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--use-supervars', action='store_true')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('citup_snv_fn')
  parser.add_argument('citup_vid_fn')
  parser.add_argument('citup_clusters_fn')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.ssm_fn)
  params = inputparser.load_params(args.params_fn)
  clusters = params['clusters']

  if args.use_supervars:
    supervars = clustermaker.make_cluster_supervars(clusters, variants)
    superclusters = clustermaker.make_superclusters(supervars)
    garbage = set()
    write_snvs(supervars, garbage, args.citup_snv_fn, args.citup_vid_fn)
    write_clusters(supervars, garbage, superclusters, args.citup_clusters_fn)
  else:
    garbage = set(params['garbage'])
    write_snvs(variants, garbage, args.citup_snv_fn, args.citup_vid_fn)
    write_clusters(variants, garbage, clusters, args.citup_clusters_fn)

if __name__ == '__main__':
  main()
