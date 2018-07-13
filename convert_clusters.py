import argparse
import csv
import json
import sys
from collections import defaultdict
import common
import vaf_plotter

def load_handbuilt(handbuiltfn):
  with open(handbuiltfn) as F:
    J = json.load(F)
  return J

def load_steph_clusters(stephfn):
  clusters = defaultdict(set)
  with open(stephfn) as F:
    reader = csv.DictReader(F)
    for row in reader:
      C = row['Cluster']
      if C == '.CNA':
        continue
      elif C == 'Garbage' or C == 'garbage1':
        continue
      C = int(C)

      chrom = row['Chrom']
      if chrom.startswith('chr'):
        chrom = chrom[3:]
      pos = int(row['Pos'])
      clusters[C].add((chrom, pos))
  return clusters

def check_vars(hb_clusters, hb_garbage, steph_clusters, variants):
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

    with open('intersection.txt', 'w') as F:
      intersection = hb_vars & steph_vars
      for V in intersection:
        print('%s_%s' % V, file=F)

    sys.exit(1)

def map_variants_to_clusters(steph_clusters, variants):
  varidmap = {(V['chrom'], V['pos']): int(V['id'][1:]) for V in variants.values()}

  existing = set(varidmap.keys())
  steph_vars = set([S for C in steph_clusters for S in steph_clusters[C]])
  print(
    len(existing),
    len(steph_vars),
    len(existing & steph_vars),
    len(steph_vars - existing),
    len(existing - steph_vars),
  )

  remapped = {C: sorted([varidmap[S] for S in steph_clusters[C] if S in varidmap]) for C in steph_clusters.keys()}
  # Add expliict empty cluster for tree root.
  remapped[0] = []
  assert set(remapped.keys()) == set(range(len(steph_clusters) + 1))
  return remapped

def write_varlist(S, outfn):
  S = sorted(S)
  with open(outfn, 'w') as F:
    print('Chrom,Pos', file=F)
    print(*['%s,%s' % V for V in S], file=F, sep='\n')

def convert_hb_to_steph(hb_clusters, hb_garbage, variants):
  def _print_row(ssmid, cidx):
    V = variants['s%s' % ssmid]
    vals = [
      find_gene_name(V['chrom'], V['pos']),
      V['chrom'],
      str(V['pos']),
      str(cidx),
    ]
    print(','.join(vals))

  outf = sys.stdout
  print('Gene,Chrom,Pos,Cluster', file=outf)
  for cidx, C in enumerate(hb_clusters):
    for ssmid in C:
      _print_row(ssmid, cidx)
  for ssmid in hb_garbage:
    _print_row(ssmid, 'garbage')

def find_gene_name(chrom, pos):
  if not hasattr(find_gene_name, 'spreadsheet_rows'):
    raise Exception('find_gene_name not initialized')
  return vaf_plotter.find_gene_name(chrom, pos, find_gene_name.spreadsheet_rows)

def write_modified(remapped_clusters, hb, tree_type, outfn):
  hb[tree_type]['clusters'] = remapped_clusters
  hb[tree_type]['garbage'] = []

  N = len(remapped_clusters) - 1
  linear = {idx: [idx + 1] for idx in range(N)}
  hb[tree_type]['structure'] = linear

  with open(outfn, 'w') as F:
    json.dump(hb, F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--direction', choices=('steph_to_hb', 'hb_to_steph'), required=True)
  parser.add_argument('--sampid', required=True)
  parser.add_argument('--ssms', dest='ssmfn', required=True)
  parser.add_argument('--tree-type', required=True)
  parser.add_argument('--steph-clusters', dest='steph_clusters_fn')
  parser.add_argument('--handbuilt', dest='handbuilt_fn', required=True)
  parser.add_argument('--spreadsheet', dest='spreadsheet_fn', required=True)
  parser.add_argument('--output', dest='out_fn')
  args = parser.parse_args()

  find_gene_name.spreadsheet_rows = vaf_plotter.load_spreadsheet(args.spreadsheet_fn)

  hb = load_handbuilt(args.handbuilt_fn)
  hb_clusters = hb[args.tree_type]['clusters']
  hb_garbage = hb[args.tree_type]['garbage']
  hb_garage = hb[args.tree_type]['garbage']
  variants = common.parse_ssms(args.sampid, args.ssmfn)

  if args.direction == 'steph_to_hb':
    steph_clusters = load_steph_clusters(args.steph_clusters_fn)
    #check_vars(hb_clusters, hb_garbage, steph_clusters, variants)
    remapped = map_variants_to_clusters(steph_clusters, variants)
    write_modified(remapped, hb, args.tree_type, args.out_fn)
  else:
    convert_hb_to_steph(hb_clusters, hb_garbage, variants)


main()