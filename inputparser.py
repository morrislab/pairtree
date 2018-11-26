import csv
import json
import numpy as np

def _extract_nums(S, dtype):
  return np.array(S.split(','), dtype=dtype)

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
      variant['omega_v'] = np.maximum(variant['omega_v'], 1e-5)
      variant['omega_v'] = np.minimum(variant['omega_v'], 1 - 1e-5)

      variant['ref_reads'] = variant['total_reads'] - variant['var_reads']
      variant['vaf'] = variant['var_reads'] / variant['total_reads']
      variant['chrom'], variant['pos'] = variant['name'].split('_')
      variant['pos'] = int(variant['pos'])

      assert variant['id'] not in variants
      variants[variant['id']] = variant

  S = len(next(iter(variants.values()))['var_reads'])
  for vid, V in variants.items():
    for K in ('var_reads', 'total_reads', 'omega_v'):
      assert len(V[K]) == S, '%s for %s is of length %s, but %s expected' % (K, vid, len(V[K]), S)

  vids = set([int(vid[1:]) for vid in variants.keys()])
  assert vids == set(range(len(variants)))
  return variants

def load_params(paramsfn):
  if paramsfn is None:
    return {}
  with open(paramsfn) as P:
    return json.load(P)

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
