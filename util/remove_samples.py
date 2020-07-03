import argparse
import json
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import inputparser

def _process(ssmfn, jsonfn, to_remove):
  params = inputparser.load_params(jsonfn)
  ssms = inputparser.load_ssms(ssmfn)

  to_remove = set([int(idx) for idx in to_remove.split(',')])
  N = len(params['samples'])
  all_samps = set(range(N))
  assert to_remove.issubset(all_samps)
  to_keep = sorted(all_samps - to_remove)
  assert len(to_keep) > 0

  params['samples'] = [params['samples'][idx] for idx in to_keep]
  for vid in ssms.keys():
    for K in ('var_reads', 'ref_reads', 'total_reads', 'vaf', 'omega_v'):
      ssms[vid][K] = ssms[vid][K][to_keep]

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
  parser.add_argument('to_remove')
  args = parser.parse_args()

  _process(args.ssm_fn, args.json_fn, args.to_remove)

if __name__ == '__main__':
  main()
