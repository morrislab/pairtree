import sys
import os
import numpy as np
import argparse

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import clustermaker
import inputparser

def extract_matrix(variants, key):
  vids = sorted(variants.keys(), key = lambda vid: int(vid[1:]))
  arr = np.array([variants[K][key] for K in vids])
  return (vids, arr)

def _write_inputs(vids, sampnames, var_reads, ref_reads, resultfn):
  K, S = var_reads.shape
  assert ref_reads.shape == (K,S)
  assert len(sampnames) == S
  assert len(vids) == K

  sep = '\t'
  with open(resultfn, 'w') as F:
    header = [V for V in vids for _ in range(2)]
    print(*header, sep=sep, file=F)

    for sidx in range(S):
      outrow = [sampnames[sidx]]
      for kidx in range(K):
        for mat in (ref_reads, var_reads):
          outrow.append(str(mat[kidx,sidx]))
      print(*outrow, sep=sep, file=F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('calder_input_fn')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.ssm_fn)
  params = inputparser.load_params(args.params_fn)
  clusters = params['clusters']
  supervars = clustermaker.make_cluster_supervars(clusters, variants)

  vids1, var_reads = extract_matrix(supervars, 'var_reads')
  vids2, ref_reads = extract_matrix(supervars, 'ref_reads')
  assert vids1 == vids2
  vids = vids1

  _write_inputs(vids, params['samples'], var_reads, ref_reads, args.calder_input_fn)

if __name__ == '__main__':
  main()
