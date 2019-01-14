import argparse
import json
import csv
from collections import defaultdict

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import inputparser

def convert_clusters(scresultsfn, varid_map):
  clusters = defaultdict(list)
  garbage = []

  with open(scresultsfn) as F:
    R = csv.DictReader(F, delimiter='\t')
    for row in R:
      chrom, pos, cluster = row['chr'], int(row['st']), row['cluster']
      varid = varid_map['%s_%s' % (chrom, pos)]
      if cluster == 'NA':
        garbage.append(varid)
      else:
        cluster = int(cluster)
        clusters[cluster].append(varid)

  cids = sorted(clusters.keys())
  assert set(cids) == set(range(1, len(cids) + 1))
  clusters = [clusters[cid] for cid in cids]
  return (clusters, garbage)

def build_variant_to_varid_map(variants):
  varid_map = {'%s_%s' % (V['chrom'], V['pos']): int(V['id'][1:]) for V in variants.values()}
  # Ensure no duplicate entries exist.
  assert len(varid_map) == len(variants)
  return varid_map

def add_missing_sex_variants_to_garbage(variants, clusters, garbage):
  # I run SciClone without sex variants, since I don't know how to specify
  # total numbers of locus according to their inputs -- maybe I need to make a
  # quasi-CNA covering all of X and Y in males, but I looked into this and
  # couldn't figure it out. As such, just mark all sex variants as garbage.
  existing = set([V for C in clusters for V in C] + list(garbage))
  vids = sorted([int(V[1:]) for V in variants.keys()])
  for vid, var in variants.items():
    vid = int(vid[1:])
    if vid in existing:
      continue
    assert var['chrom'] in ('X', 'Y')
    garbage.append(vid)

def write_results(clusters, garbage, params_fn_orig, params_fn_modified):
  params = inputparser.load_params(params_fn_orig)
  for K in ('clusters', 'garbage'):
    if K in params:
      del params[K]
  params['clusters'] = clusters
  params['garbage'] = garbage

  with open(params_fn_modified, 'w') as F:
    json.dump(params, F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('ssm_fn')
  parser.add_argument('scresults_fn')
  parser.add_argument('params_fn_orig')
  parser.add_argument('params_fn_modified')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.ssm_fn)
  varid_map = build_variant_to_varid_map(variants)
  clusters, garbage = convert_clusters(args.scresults_fn, varid_map)
  add_missing_sex_variants_to_garbage(variants, clusters, garbage)
  write_results(clusters, garbage, args.params_fn_orig, args.params_fn_modified)

main()
