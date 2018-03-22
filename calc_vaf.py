import common
import argparse
import json
import numpy as np

def load_sampnames(paramsfn):
  with open(paramsfn) as P:
    params = json.load(P)
  sampnames = params['samples']
  return sampnames

def extract_vaf(variants):
  vaf = np.array([V['vaf'] for V in variants.values()])
  return vaf

def remove_samples(vaf, sampnames):
  M, S = vaf.shape
  kept_indices = []
  kept_samples = []
  for idx, sampname in enumerate(sampnames):
    if 'xeno' in sampname.lower():
      continue
    kept_indices.append(idx)
    kept_samples.append(sampname)
  return (vaf[:,kept_indices], kept_samples)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sampid')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  args = parser.parse_args()

  variants = common.parse_ssms(args.sampid, args.ssm_fn)
  sampnames = load_sampnames(args.params_fn)

  vaf = extract_vaf(variants)
  vaf, sampnames = remove_samples(vaf, sampnames)
  from IPython import embed
  embed()

main()
