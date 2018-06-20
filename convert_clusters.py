import argparse
import csv
import json
import sys
from collections import defaultdict
import common

def load_handbuilt(handbuiltfn):
  with open(handbuiltfn) as F:
    J = json.load(F)
  return J

def load_steph_clusters(stephfn):
  clusters = defaultdict(set)
  with open(stephfn) as F:
    reader = csv.DictReader(F)
    for row in reader:
      C = int(row['Cluster'])
      chrom = row['Chrom']
      assert chrom.startswith('chr')
      chrom = chrom[3:]
      pos = int(row['Pos'])
      clusters[C].add((chrom, pos))
  return clusters

def remap_clusters(hb_clusters, hb_garbage, steph_clusters, variants):
  hb_clusters_and_garb = hb_clusters + [hb_garbage]
  hb_ssmids  = set(['s%s' % S for C in hb_clusters_and_garb for S in C])
  hb_vars    = set([(variants[S]['chrom'], variants[S]['pos']) for S in hb_ssmids])
  steph_vars = set([V for C in steph_clusters.values() for V in C])

  if hb_vars != steph_vars:
    print('hb_vars - steph_vars', hb_vars - steph_vars)
    print('steph_vars - hb_vars', steph_vars - hb_vars)
    print('len(steph_vars)=%s len(hb_vars)=%s len(steph_vars - hb_vars)=%s len(hb_vars - steph_vars)=%s len(intersection)=%s' % (
      len(steph_vars),
      len(hb_vars),
      len(steph_vars - hb_vars),
      len(hb_vars - steph_vars),
      len(steph_vars & hb_vars),
    ))
    write_varlist(steph_vars - hb_vars, 'SJETV010.steph_minus_jeff.csv')
    write_varlist(hb_vars - steph_vars, 'SJETV010.jeff_minus_steph.csv')
    write_varlist(hb_vars & steph_vars, 'SJETV010.common.csv')
    sys.exit(1)

def write_varlist(S, outfn):
  S = sorted(S)
  with open(outfn, 'w') as F:
    print('Chrom,Pos', file=F)
    print(*['%s,%s' % V for V in S], file=F, sep='\n')

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sampid')
  parser.add_argument('ssmfn')
  parser.add_argument('tree_type')
  parser.add_argument('steph_clustersfn')
  parser.add_argument('handbuiltfn')
  parser.add_argument('outfn')
  args = parser.parse_args()

  hb = load_handbuilt(args.handbuiltfn)
  hb_clusters = hb[args.tree_type]['clusters']
  hb_garbage = hb[args.tree_type]['garbage']
  hb_garage = hb[args.tree_type]['garbage']
  steph_clusters = load_steph_clusters(args.steph_clustersfn)
  variants = common.parse_ssms(args.sampid, args.ssmfn)

  remap_clusters(hb_clusters, hb_garbage, steph_clusters, variants)

main()
