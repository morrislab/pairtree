import argparse
import common
import json
import os

def load_sampnames(paramsfn):
  with open(paramsfn) as P:
    params = json.load(P)
  sampnames = params['samples']
  return sampnames

def convert(variants, sampnames, outdir):
  vids = variants.keys()
  for sidx, S in enumerate(sampnames):
    with open(os.path.join(outdir, '%s.dat' % S), 'w') as F:
      print(*['chr', 'pos', 'ref_reads', 'var_reads', 'vaf'], sep='\t', file=F)
      for vid in vids:
        vals = (
          variants[vid]['chrom'],
          variants[vid]['pos'],
          int(variants[vid]['ref_reads'][sidx]),
          int(variants[vid]['var_reads'][sidx]),
          variants[vid]['vaf'][sidx],
        )
        print(*vals, sep='\t', file=F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sampid')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('out_dir')
  args = parser.parse_args()

  variants = common.parse_ssms(args.sampid, args.ssm_fn)
  sampnames = load_sampnames(args.params_fn)
  convert(variants, sampnames, args.out_dir)

main()
