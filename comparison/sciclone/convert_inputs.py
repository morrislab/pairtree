import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import argparse
import inputparser

def convert(variants, sampnames, outdir):
  vids = variants.keys()
  for sidx, S in enumerate(sampnames):
    with open(os.path.join(outdir, '%s.dat' % S), 'w') as F:
      print(*['chr', 'pos', 'ref_reads', 'var_reads', 'vaf'], sep='\t', file=F)
      for vid in vids:
        vals = (
          variants[vid]['chrom'],
          variants[vid]['pos'],
          variants[vid]['ref_reads'][sidx],
          variants[vid]['var_reads'][sidx],
          variants[vid]['vaf'][sidx],
        )
        print(*vals, sep='\t', file=F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('out_dir')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.ssm_fn)
  params = inputparser.load_params(args.params_fn)
  sampnames = params['samples']

  convert(variants, sampnames, args.out_dir)

main()
