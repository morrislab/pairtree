import csv
import json
import numpy as np

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
        'var_reads': np.array([int(V) for V in row['var_reads'].split(',')], dtype=np.int),
        'total_reads': np.array([int(V) for V in row['total_reads'].split(',')], dtype=np.int),
        'omega_v': float(row['var_read_prob']),
      }

      assert np.all(variant['total_reads'] >= variant['var_reads'])
      assert 0 <= variant['omega_v'] <= 1

      variant['ref_reads'] = variant['total_reads'] - variant['var_reads']
      variant['vaf'] = variant['var_reads'] / variant['total_reads']
      variant['chrom'], variant['pos'] = variant['name'].split('_')
      variant['pos'] = int(variant['pos'])

      assert variant['id'] not in variants
      variants[variant['id']] = variant

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
      for K in ('var_reads', 'total_reads'):
        V[K] = ','.join([str(R) for R in V[K]])
      V['var_read_prob'] = V['omega_v']
      print(*[V[K] for K in keys], sep='\t', file=outf)
