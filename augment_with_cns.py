import argparse
import json
import csv
import sys
from collections import OrderedDict

def parse_ssms(ssm_fn):
  ssms = OrderedDict()
  with open(ssm_fn) as ssmf:
    reader = csv.DictReader(ssmf, delimiter='\t')
    for ssm in reader:
      for key in ('a', 'd'):
        ssm[key] = [int(V) for V in ssm[key].split(',')]
      assert ssm['gene'] not in ssms
      ssms[ssm['gene']] = ssm
  return (ssms, reader.fieldnames)

def parse_samples(paramsfn):
  with open(paramsfn) as F:
    params = json.load(F)
  return params['samples']

def add_samples(dataset, orig_ssm_fn, orig_params_fn, augment_ssm_fn, augment_params_fn, should_add):
  orig_ssms, orig_fields = parse_ssms(orig_ssm_fn)
  augmented_ssms, augmented_fields = parse_ssms(augment_ssm_fn)
  orig_samples, augmented_samples = parse_samples(orig_params_fn), parse_samples(augment_params_fn)
  assert orig_fields == augmented_fields

  augment_mask = []
  for aug_sampname in augmented_samples:
    augment_mask.append(should_add(dataset, aug_sampname))
  num_aug_samples = len([V for V in augment_mask if V])

  added_samples = [samp for samp, should_include in zip(augmented_samples, augment_mask) if should_include]
  print(dataset, 'adding', added_samples, file=sys.stderr)
  orig_samples += added_samples

  for pos, orig_ssm in orig_ssms.items():
    if pos not in augmented_ssms:
      print(dataset, 'missing', orig_ssm['id'], pos, 'in augmented', file=sys.stderr)
      for key in ('a', 'd'):
        orig_ssm[key] += num_aug_samples * [1]
      continue
    augssm = augmented_ssms[pos]
    for key in ('a', 'd'):
      added = [val for val, should_include in zip(augssm[key], augment_mask) if should_include]
      assert len(added) == len(added_samples)
      orig_ssm[key] += added

  return (orig_fields, orig_ssms, orig_samples)

def write_munged(fields, ssms, samples, ssmfn, sampfn):
  params = {'samples': samples}
  with open(sampfn, 'w') as S:
    json.dump(params, S)

  with open(ssmfn, 'w') as F:
    print(*fields, sep='\t', file=F)
    for ssm in ssms.values():
      for K in ('a', 'd'):
        ssm[K] = ','.join([str(val) for val in ssm[K]])
      vals = [ssm[K] for K in fields]
      print(*vals, sep='\t', file=F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('dataset')
  parser.add_argument('orig_ssm_fn')
  parser.add_argument('orig_params_fn')
  parser.add_argument('augment_ssm_fn')
  parser.add_argument('augment_params_fn')
  args = parser.parse_args()

  should_add = lambda dataset, sampname: 'Diagnosis' in sampname and 'CNS' in sampname

  fields, ssms, samples = add_samples(args.dataset, args.orig_ssm_fn, args.orig_params_fn, args.augment_ssm_fn, args.augment_params_fn, should_add)
  write_munged(fields, ssms, samples, args.orig_ssm_fn, args.orig_params_fn)

main()
