# Find variants whose maximum likelihood estimate of phi are implausibly
# large (i.e., well over 1), indicating something is wrong with such variants.
# Render these as garbage.

import argparse
import numpy as np
import scipy.special
import scipy.stats
import json

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import common
import inputparser
import util

def _calc_logbf(ssm, omega1=1):
  V = ssm['var_reads']
  R = ssm['ref_reads']
  S = len(V)

  omega = [np.broadcast_to(omega1, S), ssm['omega_v']]
  logp = []

  for om in omega:
    integ = scipy.special.betainc(V+1, R+1, om)
    # Prevent division by zero.
    loginteg = np.log(np.maximum(1e-30, integ))
    logp.append(-np.log(om) + loginteg)

  logbf = logp[0] - logp[1]
  # Convert from nats to bits.
  logbf /= np.log(2)

  return logbf

def _remove_bad(ssms, logbf_thresh, print_bad=False):
  bad = []
  bad_count = 0
  total_count = 0

  alpha0 = 1
  beta0 = 1
  phi_hat_threshold = 1 - 1e-2

  for vid, V in ssms.items():
    logbfs = _calc_logbf(V)
    is_samp_bad = logbfs > logbf_thresh

    bad_count   += np.sum(is_samp_bad)
    total_count += len(is_samp_bad)

    # Following is just used for diagnostics and is not part of decision for bad vs. good.
    phi_mle = V['var_reads'] / (V['omega_v'] * V['total_reads'])
    phi_alpha = alpha0 + V['var_reads']
    phi_beta = beta0 + np.maximum(0, V['omega_v']*V['total_reads'] - V['var_reads'])
    # How much mass in our phi_hat posterior is in the region [phi_hat_threshold, 1]?
    probs = 1. - scipy.stats.beta.cdf(phi_hat_threshold, phi_alpha, phi_beta)

    if np.any(is_samp_bad):
      bad.append(vid)
    if print_bad and np.any(is_samp_bad):
      print(vid, np.vstack((
        V['var_reads'],
        V['total_reads'],
        phi_mle,
        probs,
        is_samp_bad,
        logbfs,
      )),
      '',
      sep='\n'
    )

  bad_samp_prop = bad_count / total_count
  assert 0 <= bad_samp_prop <= 1
  return (bad, bad_samp_prop)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--logbf-threshold', type=float, default=10.,
    help='Blah')
  parser.add_argument('--print-bad', action='store_true')
  parser.add_argument('--keep-existing-garbage', action='store_true')
  parser.add_argument('in_ssm_fn')
  parser.add_argument('in_params_fn')
  parser.add_argument('out_params_fn')
  args = parser.parse_args()

  np.set_printoptions(linewidth=400, precision=3, threshold=sys.maxsize, suppress=True)
  np.seterr(divide='raise', invalid='raise', over='raise')

  ssms = inputparser.load_ssms(args.in_ssm_fn)
  params = inputparser.load_params(args.in_params_fn)
  if args.keep_existing_garbage:
    params['garbage'] = []
  else:
    ssms = inputparser.remove_garbage(ssms, params['garbage'])

  bad_vids, bad_samp_prop = _remove_bad(ssms, args.logbf_threshold, args.print_bad)
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
