import argparse
import random
import json

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import inputparser

def _is_good(cand, last):
  if last is None:
    return 'D'
  elif last == 'D':
    return cand == 'R1'
  elif last == 'R1':
    return cand.startswith('Diagnosis Xeno')
  elif last.startswith('Diagnosis Xeno'):
    return cand.startswith('Relapse Xeno')
  elif last.startswith('Relapse Xeno'):
    return cand.startswith('Diagnosis Xeno')
  else:
    raise Exception('Unknown last choice: %s' % last)

def _select_subset(sampnames, C, must_include):
  # Special-case the instance where the user wants every sample.
  if len(sampnames) == C:
    return list(sampnames)

  assert set(must_include).issubset(set(sampnames))
  # Duplicate.
  subset = list(must_include)
  candidates = set(sampnames) - set(subset)
  last = subset[-1] if len(subset) > 0 else None

  while len(subset) < C:
    cands = list(candidates)
    random.shuffle(cands)
    for cand in cands:
      if _is_good(cand, last):
        subset.append(cand)
        candidates.remove(cand)
        last = cand
        break
    else:
      raise Exception('Could not find any candidate after %s' % last)

  return subset

def _select_samp_subsets(sampnames, counts, all_must_include=None):
  subsets = []
  if all_must_include is None:
    all_must_include = []

  for C in sorted(counts):
    assert 0 < C <= len(sampnames)
    must_include = subsets[-1] if len(subsets) > 0 else all_must_include
    subsets.append(_select_subset(sampnames, C, must_include))
  return subsets

def _filter_ssms(ssms, samp_idxs):
  new_ssms = {}
  for sidx, ssm in ssms.items():
    # Duplicate so as to not modify original.
    new_ssms[sidx] = dict(ssm)
    for K in ('var_reads', 'ref_reads', 'total_reads', 'omega_v', 'vaf'):
      new_ssms[sidx][K] = ssms[sidx][K][samp_idxs]
  return new_ssms

def _find_idxs(sampnames, subset):
  idxs = [sampnames.index(mem) for mem in subset]
  return idxs

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--counts', required=True)
  parser.add_argument('in_ssm_fn')
  parser.add_argument('in_params_fn')
  parser.add_argument('out_base')
  args = parser.parse_args()

  random.seed(1337)

  counts = [int(C) for C in args.counts.split(',')]
  assert len(counts) == len(set(counts))
  ssms = inputparser.load_ssms(args.in_ssm_fn)
  params = inputparser.load_params(args.in_params_fn)
  sampnames = params['samples']

  # Always include diagnosis sample, on assumption we're working with
  # SJbALL022609 from Steph for the paper congraph figure.
  subsets = _select_samp_subsets(sampnames, counts, all_must_include=['D'])
  for subset in subsets:
    idxs = _find_idxs(sampnames, subset)
    new_ssms = _filter_ssms(ssms, idxs)
    new_params = dict(params)
    new_params['samples'] = subset

    out_base = '%s_S%s' % (args.out_base, len(subset))
    inputparser.write_ssms(new_ssms, out_base + '.ssm')
    with open(out_base + '.params.json', 'w') as F:
      json.dump(new_params, F)

if __name__ == '__main__':
  main()
