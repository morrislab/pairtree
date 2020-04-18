# Use instances where the maximum likelihood estimate of phi are implausibly
# large (i.e., well over 1) to indicate that uncalled copy-number losses of the
# normal allele occurred, which highlight instances where we should increase
# `omega.`
import argparse
import numpy as np
import scipy.stats

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import inputparser

def _fix_omegas(ssms, print_bad=False):
  # Note that SSMs are modified in place.
  percentile = 1e-10
  vaf_threshold = 0.5
  # To flag only extreme cases, set this above 1.0 -- e.g., to 1.5.
  phi_mle_threshold = 1.0
  fixed_omega = 1.0

  alpha0 = 0.5
  beta0 = 0.5
  bad = 0
  total = 0

  for vid, V in ssms.items():
    vaf_alpha = alpha0 + V['var_reads']
    vaf_beta = beta0 + V['ref_reads']
    phi_mle = V['var_reads'] / (V['omega_v'] * V['total_reads'])

    bad_omega = np.logical_and.reduce((
      # Only flag variants that haven't already had `omega_v` adjusted.
      np.isclose(0.5, V['omega_v']),
      # Is the true VAF extremely unlikely to be less than 0.5?
      scipy.stats.beta.cdf(vaf_threshold, vaf_alpha, vaf_beta) < percentile,
      # Is the phi MLE likely to be too high?
      phi_mle > phi_mle_threshold,
    ))

    if print_bad and np.any(bad_omega):
      print(np.vstack((
        V['var_reads'][bad_omega],
        V['total_reads'][bad_omega],
        phi_mle[bad_omega],
      )))
      print('')

    V['omega_v'][bad_omega] = fixed_omega
    bad += np.sum(bad_omega)
    total += len(bad_omega)

  fixed_prop = bad / total
  return fixed_prop

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('in_ssm_fn')
  parser.add_argument('out_ssm_fn')
  args = parser.parse_args()

  np.set_printoptions(linewidth=400, precision=3, threshold=sys.maxsize, suppress=True)
  np.seterr(divide='raise', invalid='raise', over='raise')

  ssms = inputparser.load_ssms(args.in_ssm_fn)
  fixed_prop = _fix_omegas(ssms, print_bad=False)
  print('fixed_omegas=%s' % fixed_prop)
  inputparser.write_ssms(ssms, args.out_ssm_fn)

if __name__ == '__main__':
  main()
