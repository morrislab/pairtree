import csv
import numpy as np
import argparse
import json
from collections import OrderedDict

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import inputparser

# Convert from PWGS format to Pairtree format, discarding garbage variants in
# the process.

def load_phylowgs(pwgs_fn):
  variants = OrderedDict()

  with open(pwgs_fn) as F:
    reader = csv.DictReader(F, delimiter='\t')
    for row in reader:
      row['name'] = row['gene']
      row['ref_reads'] = np.array([float(V) for V in row['a'].split(',')], dtype=np.int)
      row['total_reads'] = np.array([float(V) for V in row['d'].split(',')], dtype=np.int)
      row['omega_v'] = 1 - float(row['mu_v'])

      assert np.all(row['total_reads'] >= row['ref_reads'])
      row['var_reads'] = row['total_reads'] - row['ref_reads']
      assert 0 <= row['omega_v'] <= 1
      variants[row['id']] = row

  return variants

def remove_garbage(variants, garbage_ids, clusters):
  garbage_variants = {}
  N = len(variants)
  for varid in sorted(variants.keys(), key = lambda S: int(S[1:])):
    if int(varid[1:]) in garbage_ids:
      garbage_variants[varid] = variants[varid]
      del variants[varid]
  assert len(variants) == N - len(garbage_ids)
  return garbage_variants

def make_varids_contiguous(variants, garbage_ids, clusters):
  mapping = {}
  for idx, old_varid in enumerate(sorted(variants.keys(), key = lambda S: int(S[1:]))):
    old_varidx = int(old_varid[1:])
    new_varidx = idx
    mapping[old_varidx] = new_varidx

  new_variants = {'s%s' % mapping[int(V[1:])]: variants[V] for V in variants.keys()}
  for V in new_variants.keys():
    new_variants[V]['id'] = V

  new_clusters = [sorted([mapping[V] for V in C]) for C in clusters]

  assert set([int(V[1:]) for V in new_variants.keys()]) == \
    set([V for C in new_clusters for V in C]) == \
    set([int(V['id'][1:]) for V in new_variants.values()])
  assert len(new_clusters[0]) == 0 and not np.any(np.array([len(C) for C in new_clusters[1:]]) == 0)

  return (new_variants, new_clusters)

def load_handbuilt(hbfn, tree_type):
  with open(hbfn) as F:
    return json.load(F)[tree_type]

def write_pairtree_params(sampnames, clusters, structure, pairtree_params_fn):
  params = {
    'samples': sampnames,
    'clusters': clusters,
    'structure': structure,
  }
  with open(pairtree_params_fn, 'w') as F:
    json.dump(params, F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--handbuilt', dest='handbuilt_fn', required=True)
  parser.add_argument('pwgs_ssm_fn')
  parser.add_argument('pwgs_params_fn')
  parser.add_argument('pairtree_ssm_fn')
  parser.add_argument('pairtree_params_fn')
  args = parser.parse_args()

  tree_type = 'handbuilt.xeno'
  hb = load_handbuilt(args.handbuilt_fn, tree_type)
  clusters = hb['clusters']
  garbage = hb['garbage']
  structure = hb['structure']

  pwgs_params = inputparser.load_params(args.pwgs_params_fn)
  variants = load_phylowgs(args.pwgs_ssm_fn)
  remove_garbage(variants, garbage, clusters)
  variants, clusters = make_varids_contiguous(variants, garbage, clusters)

  inputparser.write_ssms(variants, args.pairtree_ssm_fn)
  write_pairtree_params(pwgs_params['samples'], clusters, structure, args.pairtree_params_fn)

if __name__ == '__main__':
  main()
