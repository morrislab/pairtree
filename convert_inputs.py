import csv
import numpy as np
import argparse
import inputparser
from collections import OrderedDict

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

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('pwgsfn')
  parser.add_argument('pairtreefn')
  args = parser.parse_args()

  variants = load_phylowgs(args.pwgsfn)
  inputparser.write_ssms(variants, args.pairtreefn)

if __name__ == '__main__':
  main()
