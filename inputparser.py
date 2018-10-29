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
        'var_reads': np.array([float(V) for V in row['var_reads'].split(',')], dtype=np.int),
        'total_reads': np.array([float(V) for V in row['total_reads'].split(',')], dtype=np.int),
        'omega_v': float(row['var_read_prob']),
      }
      assert np.all(variant['total_reads'] >= variant['var_reads'])
      variant['ref_reads'] = variant['total_reads'] - variant['var_reads']
      variant['vaf'] = variant['var_reads'] / variant['total_reads']
      variant['chrom'], variant['pos'] = variant['name'].split('_')
      variant['pos'] = int(variant['pos'])
      variants[row['id']] = variant

  return variants

def load_params(paramsfn):
  with open(paramsfn) as P:
    return json.load(P)
