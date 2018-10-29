import csv
import numpy as np
import argparse

def load_phylowgs(pwgs_fn):
  variants = []

  with open(pwgs_fn) as F:
    reader = csv.DictReader(F, delimiter='\t')
    for row in reader:
      row['name'] = row['gene']
      row['ref_reads'] = np.array([float(V) for V in row['a'].split(',')], dtype=np.int)
      row['total_reads'] = np.array([float(V) for V in row['d'].split(',')], dtype=np.int)
      row['var_read_prob'] = 1 - float(row['mu_v'])

      assert np.all(row['total_reads'] >= row['ref_reads'])
      row['var_reads'] = row['total_reads'] - row['ref_reads']
      assert 0 <= row['var_read_prob'] <= 1
      variants.append(row)

  return variants

def write_pairtree(variants, pairtree_fn):
  keys = ('id', 'name', 'var_reads', 'total_reads', 'var_read_prob')
  with open(pairtree_fn, 'w') as outf:
    print(*keys, sep='\t', file=outf)
    for V in variants:
      for K in ('var_reads', 'total_reads'):
        V[K] = ','.join([str(R) for R in V[K]])
      print(*[V[K] for K in keys], sep='\t', file=outf)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('pwgsfn')
  parser.add_argument('pairtreefn')
  args = parser.parse_args()

  variants = load_phylowgs(args.pwgsfn)
  write_pairtree(variants, args.pairtreefn)

if __name__ == '__main__':
  main()
