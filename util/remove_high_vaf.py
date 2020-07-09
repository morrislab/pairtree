# Find variants whose maximum likelihood estimate of phi are implausibly
# large (i.e., well over 1), indicating something is wrong with such variants.
# Render these as garbage.

import argparse
import numpy as np
import scipy.stats
import json

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import common
import inputparser

def _remove_bad(ssms, phi_hat_threshold, quantile, print_bad=False):
  # Note that SSMs are modified in-place.
  alpha0 = 1
  beta0 = 1
  # To flag only extreme cases, set this above 1.0 -- e.g., to 1.5.
  phi_mle_threshold = 1.0

  bad = []
  bad_count = 0
  total_count = 0

  for vid, V in ssms.items():
    phi_alpha = alpha0 + V['var_reads']
    phi_beta = beta0 + np.maximum(0, V['omega_v']*V['total_reads'] - V['var_reads'])
    # Binary decision: is the mass in the `phi ~ Beta(...)` CDF below
    # the `quantile` threshold?
    probs = scipy.stats.beta.cdf(phi_hat_threshold, phi_alpha, phi_beta)

    phi_mle = V['var_reads'] / (V['omega_v'] * V['total_reads'])
    bad_samps = np.logical_and.reduce((
      probs < quantile,
      phi_mle > phi_mle_threshold,
    ))

    if np.any(bad_samps):
      bad.append(vid)
    if print_bad and np.any(bad_samps):
      print(vid, np.vstack((
        V['var_reads'],
        V['total_reads'],
        phi_mle,
        probs,
        bad_samps,
      )), sep='\n')
      print('')

    bad_count += np.sum(bad_samps)
    total_count += len(bad_samps)

  bad_samp_prop = bad_count / total_count
  return (bad, bad_samp_prop)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--phi-hat-threshold', type=float, default = 1 - 1e-2,
    help='Blah')
  parser.add_argument('--quantile', type=float, default = 0.5,
    help='Blah')
  parser.add_argument('--print-bad-data', action='store_true')
  parser.add_argument('in_ssm_fn')
  parser.add_argument('in_params_fn')
  parser.add_argument('out_params_fn')
  args = parser.parse_args()

  np.set_printoptions(linewidth=400, precision=3, threshold=sys.maxsize, suppress=True)
  np.seterr(divide='raise', invalid='raise', over='raise')

  ssms = inputparser.load_ssms(args.in_ssm_fn)
  params = inputparser.load_params(args.in_params_fn)
  ssms = inputparser.remove_garbage(ssms, params['garbage'])

  bad_vids, bad_samp_prop = _remove_bad(ssms, args.phi_hat_threshold, args.quantile, args.print_bad_data)
  bad_ssm_prop = len(bad_vids) / len(ssms)
  if len(bad_vids) > 0:
    params['garbage'] = common.sort_vids(params['garbage'] + bad_vids)
    with open(args.out_params_fn, 'w') as F:
      json.dump(params, F)

  stats = {
    'bad_ssms': common.sort_vids(bad_vids),
    'bad_samp_prop': '%.3f' % bad_samp_prop,
    'bad_ssm_prop': '%.3f' % bad_ssm_prop,
  }
  for K, V in stats.items():
    print('%s=%s' % (K, V))

if __name__ == '__main__':
  main()
