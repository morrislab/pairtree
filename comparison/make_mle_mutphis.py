import argparse
import numpy as np

import mutphi
import mutstat

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import inputparser
import common

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('ssm_fn')
  parser.add_argument('mutphi_fn')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.ssm_fn)
  vids = common.extract_vids(variants)
  var_reads = np.array([variants[V]['var_reads'] for V in vids])
  total_reads = np.array([variants[V]['total_reads'] for V in vids])
  omega_v = np.array([variants[V]['omega_v'] for V in vids])

  mle_phi = (1/omega_v)*(var_reads / total_reads)
  assert np.all(0 <= mle_phi)
  mle_phi = np.minimum(1, mle_phi)

  clusters = [[V] for V in vids]
  llhs = [0]
  counts = [1]
  mphi = mutphi.calc_mutphi([mle_phi], llhs, [clusters], args.ssm_fn, counts)
  mutstat.write(mphi, args.mutphi_fn)

if __name__ == '__main__':
  main()
