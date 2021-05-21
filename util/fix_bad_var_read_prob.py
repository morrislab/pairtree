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

def _remove_bad(ssms, logbf_thresh, var_read_prob_alt, verbose=False):
  bad = []
  bad_count = 0
  total_count = 0

  alpha0 = 1
  beta0 = 1
  phi_hat_threshold = 1 - 1e-2

  for vid, V in ssms.items():
    logbfs = _calc_logbf(V, omega1=var_read_prob_alt)
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
    if verbose and np.any(is_samp_bad):
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
    description='Find variants with likely incorrect var_read_prob by comparing model with provided var_read_prob to haploid (LOH) model using Bayes factors',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--logbf-threshold', type=float, default=10.,
    help='Logarithm of Bayes factor threshold at which the haploid model is accepted as more likely model than the model using the provided var_read_prob')
  parser.add_argument('--verbose', action='store_true',
    help='Print debugging messages')
  parser.add_argument('--ignore-existing-garbage', action='store_true',
    help='Ignore any existing garbage variants listed in in_params_fn and test all variants. If not specified, any existing garbage variants will be kept as garbage and not tested again.')
  parser.add_argument('--action', choices=('add_to_garbage', 'modify_var_read_prob'), default='add_to_garbage')
  parser.add_argument('--var-read-prob-alt', type=float, default=1.)
  parser.add_argument('in_ssm_fn',
    help='Input SSM file with mutations')
  parser.add_argument('in_params_fn',
    help='Input params file listing sample names and any existing garbage mutations')
  parser.add_argument('out_ssm_fn',
    help='Output SSM file with modified list of garbage mutations')
  parser.add_argument('out_params_fn',
    help='Output params file with modified list of garbage mutations')
  args = parser.parse_args()

  np.set_printoptions(linewidth=400, precision=3, threshold=sys.maxsize, suppress=True)
  np.seterr(divide='raise', invalid='raise', over='raise')

  if args.ignore_existing_garbage:
    variants, params = inputparser.load_ssms_and_params(args.in_ssm_fn, args.in_params_fn, remove_garb=False)
    params['garbage'] = []
  else:
    variants, params = inputparser.load_ssms_and_params(args.in_ssm_fn, args.in_params_fn)

  bad_vids, bad_samp_prop = _remove_bad(variants, args.logbf_threshold, args.var_read_prob_alt, args.verbose)
  bad_ssm_prop = len(bad_vids) / len(variants)

  if args.action == 'add_to_garbage':
    params['garbage'] = common.sort_vids(set(bad_vids) | set(params['garbage']))
  elif args.action == 'modify_var_read_prob':
    for vid in bad_vids:
      variants[vid]['omega_v'][:] = args.var_read_prob_alt
  else:
    raise Exception('Unknown action: %s' % args.action)

  inputparser.write_ssms(variants, args.out_ssm_fn)
  with open(args.out_params_fn, 'w') as F:
    json.dump(params, F)

  stats = {
    'num_bad_ssms': len(bad_vids),
    'bad_ssms': common.sort_vids(bad_vids),
    'bad_samp_prop': '%.3f' % bad_samp_prop,
    'bad_ssm_prop': '%.3f' % bad_ssm_prop,
  }
  for K, V in stats.items():
    print('%s=%s' % (K, V))

if __name__ == '__main__':
  main()
