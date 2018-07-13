import argparse
import json

def determine_retained(dataset, params, should_remove):
  nsamples = len(params['samples'])
  to_keep = nsamples * [True]

  for sidx, samp in enumerate(params['samples']):
    for sr in should_remove:
      if sr(dataset, samp):
        # If *any* filters fail, remove the sample.
        to_keep[sidx] = False

  return to_keep

def remove(dataset, ssmfn, paramsfn, should_remove):
  with open(paramsfn) as F:
    params = json.load(F)

  to_keep = determine_retained(dataset, params, should_remove)
  retained = [S for (S, K) in zip(params['samples'], to_keep) if K]
  retained_hidden = [S for S in params['hidden_samples'] if S in retained]
  print('%s: removing %s' % (dataset, set(params['samples']) - set(retained)))
  params['samples'] = retained
  params['hidden_samples'] = retained_hidden

  with open(paramsfn, 'w') as F:
    json.dump(params, F)

  with open(ssmfn) as F:
    records = [L.strip().split('\t') for L in F.readlines()]
  header, records = records[0], records[1:]

  with open(ssmfn, 'w') as F:
    print(*header, sep='\t', file=F)
    for rec in records:
      for kidx in (2, 3):
        vals = rec[kidx].split(',')
        assert len(vals) == len(to_keep)
        rec[kidx] = ','.join([val for (keep, val) in zip(to_keep, vals) if keep])
      print(*rec, sep='\t', file=F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('dataset')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  args = parser.parse_args()

  should_remove = [
    lambda dataset, sampname: ('CNS' in sampname),
    lambda dataset, sampname: (dataset == 'SJETV047' and sampname == 'R2'),
    lambda dataset, sampname: (dataset == 'SJBALL022611' and sampname == 'Relapse Xeno 11'),
    lambda dataset, sampname: (dataset == 'SJETV010steph' and sampname == 'R1'),
  ]
  remove(args.dataset, args.ssm_fn, args.params_fn, should_remove)

main()
