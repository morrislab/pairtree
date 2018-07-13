import vaf_plotter
import common
import argparse
import csv
import numpy as np
import json

def load_sampnames(paramsfn):
  with open(paramsfn) as P:
    params = json.load(P)
  sampnames = params['samples']
  return sampnames

def parse_ssms(sampid, ssmfn):
  vaf = []
  ssm_ids = []
  var_names = []
  variants = {}

  with open(ssmfn) as F:
    reader = csv.DictReader(F, delimiter='\t')
    for row in reader:
      if False and len(variants) >= 3:
        break
      variant = {
        'id': row['id'],
        'name': row['gene'],
        'gene': 'LOL',
        'ref_reads': np.array([float(V) for V in row['a'].split(',')]),
        'total_reads': np.array([float(V) for V in row['d'].split(',')]),
        'mu_v': float(row['mu_v']),
      }

      zero_total = variant['total_reads'] == 0
      # Fix NaNs
      variant['ref_reads'][zero_total] = variant['total_reads'][zero_total] = 1

      variant['var_reads'] = variant['total_reads'] - variant['ref_reads']
      variant['vaf'] = variant['var_reads'] / variant['total_reads']
      if '.' in variant['name']:
        name = variant['name'].split('.', 1)[1]
      else:
        name = variant['name']
      variant['chrom'], variant['pos'] = name.split('_', 2)
      variant['pos'] = int(variant['pos'])
      variants[row['id']] = variant

  return variants

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sampid')
  parser.add_argument('ssm_fn')
  parser.add_argument('paramsfn')
  parser.add_argument('out_fn')
  args = parser.parse_args()

  variants = parse_ssms(args.sampid, args.ssm_fn)
  sampnames = load_sampnames(args.paramsfn)
  with open(args.out_fn, 'w') as outf:
    vaf_plotter.plot_unclustered_vafs(args.sampid, variants, None, sampnames, None, outf, patient_samples_only=False)

if __name__ == '__main__':
  main()
