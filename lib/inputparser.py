import csv
import json
import numpy as np
import numpy.ma as ma
import common

def _extract_nums(S, dtype):
  return np.array(S.split(','), dtype=dtype)

def _extract_mat(variants, key):
  vids = common.extract_vids(variants)
  arr = [variants[vid][key] for vid in vids]
  return np.array(arr)

def load_read_counts(variants):
  vids = common.extract_vids(variants)
  V = _extract_mat(variants, 'var_reads')
  T = _extract_mat(variants, 'total_reads')
  omega = _extract_mat(variants, 'omega_v')
  T_prime = np.maximum(V, omega*T)
  return (vids, V, T, T_prime, omega)

def load_ssms(ssmfn, max_ssms=None):
  variants = {}

  with open(ssmfn) as F:
    reader = csv.DictReader(F, delimiter='\t')
    for row in reader:
      if max_ssms is not None and len(variants) >= max_ssms:
        break
      variant = {
        'id': row['id'],
        'name': row['name'],
        'var_reads': _extract_nums(row['var_reads'], np.int),
        'total_reads': _extract_nums(row['total_reads'], np.int),
        'omega_v': _extract_nums(row['var_read_prob'], np.float),
      }

      assert np.all(variant['total_reads'] >= variant['var_reads'])
      assert np.all(0 <= variant['omega_v']) and np.all(variant['omega_v'] <= 1)

      variant['ref_reads'] = variant['total_reads'] - variant['var_reads']
      T = ma.masked_equal(variant['total_reads'], 0)
      variant['vaf'] = np.array(variant['var_reads'] / T)

      assert variant['id'] not in variants
      variants[variant['id']] = variant

  S = len(next(iter(variants.values()))['var_reads'])
  for vid, V in variants.items():
    for K in ('var_reads', 'total_reads', 'omega_v'):
      assert len(V[K]) == S, '%s for %s is of length %s, but %s expected' % (K, vid, len(V[K]), S)

  vids = set([int(vid[1:]) for vid in variants.keys()])
  return variants

def remove_garbage(variants, garbage):
  assert set(garbage).issubset(set(variants.keys()))
  munged = dict(variants)
  for vid in garbage:
    del munged[vid]
  return munged

def load_params(paramsfn):
  if paramsfn is None:
    return {}
  with open(paramsfn) as P:
    return json.load(P)

def load_ssms_and_params(ssmfn, paramsfn, remove_garb=True):
  variants = load_ssms(ssmfn)
  params = load_params(paramsfn)
  if 'garbage' not in params:
    params['garbage'] = []
  if remove_garb:
    variants = remove_garbage(variants, params['garbage'])
  return (variants, params)

def write_ssms(variants, ssm_fn):
  keys = ('id', 'name', 'var_reads', 'total_reads', 'var_read_prob')

  with open(ssm_fn, 'w') as outf:
    print(*keys, sep='\t', file=outf)
    for V in variants.values():
      V = dict(V) # Don't modify original variant.
      for K in ('var_reads', 'total_reads', 'omega_v'):
        V[K] = ','.join([str(R) for R in V[K]])
      V['var_read_prob'] = V['omega_v']
      print(*[V[K] for K in keys], sep='\t', file=outf)
