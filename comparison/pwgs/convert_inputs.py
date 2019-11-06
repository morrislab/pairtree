import sys
import os
import argparse
import json
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import inputparser
import clustermaker

def write_ssms(variants, outfn):
  _stringify = lambda A: ','.join([str(V) for V in A])
  mu_r = 0.999
  cols = ('id', 'gene', 'a', 'd', 'mu_r', 'mu_v')

  with open(outfn, 'w') as outf:
    print(*cols, sep='\t', file=outf)
    for V in variants.values():
      assert len(set(V['omega_v'])) == 1

      variant = {
        'id': 's%s' % int(V['id'][1:]),
        'gene': V['name'],
        'a': _stringify(V['ref_reads']),
        'd': _stringify(V['total_reads']),
        'mu_r': mu_r,
        'mu_v': np.mean(1 - V['omega_v']),
      }
      print(*[variant[K] for K in cols], sep='\t', file=outf)

def write_params(sampnames, outfn):
  with open(outfn, 'w') as outf:
    json.dump({'samples': sampnames}, outf)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--use-supervars', dest='use_supervars', action='store_true')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('pwgs_ssm_fn')
  parser.add_argument('pwgs_params_fn')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.ssm_fn)
  params = inputparser.load_params(args.params_fn)

  if args.use_supervars:
    variants = clustermaker.make_cluster_supervars(params['clusters'], variants)
  write_ssms(variants, args.pwgs_ssm_fn)
  write_params(params['samples'], args.pwgs_params_fn)

if __name__ == '__main__':
  main()
