import argparse
import json
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import inputparser

def _process(ssmfn, jsonfn, order):
  params = inputparser.load_params(jsonfn)
  ssms = inputparser.load_ssms(ssmfn)

  order = [int(idx) for idx in order.split(',')]
  N = len(params['samples'])
  assert set(range(N)) == set(order)
  assert len(list(ssms.values())[0]['var_reads']) == N

  params['samples'] = [params['samples'][idx] for idx in order]
  for vid in ssms.keys():
    for K in ('var_reads', 'ref_reads', 'total_reads', 'vaf', 'omega_v'):
      ssms[vid][K] = ssms[vid][K][order]

  with open(jsonfn, 'w') as F:
    json.dump(params, F)
  inputparser.write_ssms(ssms, ssmfn)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('ssm_fn')
  parser.add_argument('json_fn')
  parser.add_argument('order')
  args = parser.parse_args()

  _process(args.ssm_fn, args.json_fn, args.order)

if __name__ == '__main__':
  main()
