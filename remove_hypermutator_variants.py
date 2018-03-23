# Written to remove hypermutator variants from SJETV010, where ~1/3 of variants
# have a non-zero VAF in only R2. We don't care about these, so we don't want
# to cluster them or work with them.
#     cd ~/work/steph/data/inputs/steph.xenos.nocns
#     python3 ~/work/steph/protocols/pairwise/remove_hypermutator_variants.py SJETV010 SJETV010.{sampled.ssm,params.json} SJETV010nohypermut.ssm
#     for ext in params.json sampled.cnv; do cp -a SJETV010.$ext SJETV010nohypermut.$ext; done

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

def remove_hypermutator_vars(variants, sampnames, nonzero_samps):
  varids = list(variants.keys())
  nonzero_sidxs = [sampnames.index(S) for S in nonzero_samps]
  for varid in varids:
    if np.all(variants[varid]['var_reads'][nonzero_sidxs] == 0):
      del variants[varid]

def write_modified(variants, outfn):
  with open(outfn, 'w') as F:
    print('id', 'gene', 'a', 'd', 'mu_r', 'mu_v', sep='\t', file=F)
    new_varid = 0
    varids = sorted(variants.keys(), key = lambda K: int(K[1:]))
    for varid in varids:
      V = variants[varid]
      print(*[
        's%s' % new_varid,
        '%s_%s' % (V['chrom'], V['pos']),
        ','.join([str(int(R)) for R in V['ref_reads']]),
        ','.join([str(int(R)) for R in V['total_reads']]),
        V['mu_r'],
        V['mu_v']
      ], sep='\t', file=F)
      new_varid += 1

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sampid')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('out_ssm_fn')
  args = parser.parse_args()

  variants = common.parse_ssms(args.sampid, args.ssm_fn)
  sampnames = load_sampnames(args.params_fn)

  nonzero_samps = ('D', 'R1')
  remove_hypermutator_vars(variants, sampnames, nonzero_samps)
  write_modified(variants, args.out_ssm_fn)

main()
