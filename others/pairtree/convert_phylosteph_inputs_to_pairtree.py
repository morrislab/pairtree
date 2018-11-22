import csv
import numpy as np
import argparse
import json
from collections import OrderedDict

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import inputparser
import common

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
      S = len(row['total_reads'])
      row['omega_v'] = (1 - float(row['mu_v'])) * np.ones(S)

      assert np.all(row['total_reads'] >= row['ref_reads'])
      row['var_reads'] = row['total_reads'] - row['ref_reads']
      assert np.all(0 <= row['omega_v']) and np.all(row['omega_v'] <= 1)
      variants[row['id']] = row

  return variants

def remove_garbage(variants, garbage):
  garbage_variants = {}
  N = len(variants)
  for varid in common.extract_vids(variants):
    if varid in garbage:
      garbage_variants[varid] = variants[varid]
      del variants[varid]
  assert len(variants) == N - len(garbage_ids)
  return garbage_variants

def make_varids_contiguous(variants, garbage, clusters):
  mapping = {}
  for new_idx, old_varid in enumerate(common.extract_vids(variants)):
    mapping[old_varid] = 's%s' % new_varidx

  new_variants = {mapping[V]: variants[V] for V in variants.keys()}
  for V in new_variants.keys():
    new_variants[V]['id'] = V

  new_clusters = [common.sort_vids([mapping[V] for V in C]) for C in clusters]

  assert set(new_variants.keys()) == \
    set([V for C in new_clusters for V in C]) == \
    set([V['id'] for V in new_variants.values()])
  assert not np.any(np.array([len(C) for C in new_clusters]) == 0)

  return (new_variants, new_clusters)

def load_handbuilt(hbfn, tree_type):
  with open(hbfn) as F:
    return json.load(F)[tree_type]

def convert_clusters(clusters):
  # Remove empty first cluster.
  assert len(clusters[0]) == 0
  return clusters[1:]

def write_pairtree_params(sampnames, garbage, clusters, structure, pairtree_params_fn):
  clusters = [['s%s' % V for V in C] for C in clusters]
  garbage = ['s%s' % V for V in garbage]
  params = {
    'samples': sampnames,
    'clusters': clusters,
    'structure': structure,
    'garbage': garbage,
  }
  with open(pairtree_params_fn, 'w') as F:
    json.dump(params, F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--discard-garbage', dest='discard_garbage', action='store_true')
  parser.add_argument('--handbuilt', dest='handbuilt_fn', required=True)
  parser.add_argument('pwgs_ssm_fn')
  parser.add_argument('pwgs_params_fn')
  parser.add_argument('pairtree_ssm_fn')
  parser.add_argument('pairtree_params_fn')
  args = parser.parse_args()

  tree_type = 'handbuilt.xeno'
  hb = load_handbuilt(args.handbuilt_fn, tree_type)
  clusters = convert_clusters(hb['clusters'])
  garbage = hb['garbage']
  # Since we remove the empty first cluster, the indexing on `structure` is now
  # a little weird -- cluster `i` is now represented by `i + 1` in `structure`.
  # That's okay.
  # TODO: I guess this means we need to change how we evaluate results in method comparison?
  structure = hb['structure']

  pwgs_params = inputparser.load_params(args.pwgs_params_fn)
  variants = load_phylowgs(args.pwgs_ssm_fn)
  if args.discard_garbage:
    remove_garbage(variants, garbage)
    variants, clusters = make_varids_contiguous(variants, garbage, clusters)
    garbage = []

  inputparser.write_ssms(variants, args.pairtree_ssm_fn)
  write_pairtree_params(pwgs_params['samples'], garbage, clusters, structure, args.pairtree_params_fn)

if __name__ == '__main__':
  main()
