import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import pairwise
import inputparser
import cluster_pairwise
import cluster_linfreq

def _make_init_clusters(variants):
  vids = common.extract_vids(variants)
  clusters = [[vid] for vid in vids]
  return clusters

def _make_coclust_logprior(prior, S, scale_prior_with_samps=True):
  assert 0 < prior <= 1
  logprior = {'garbage': -np.inf, 'cocluster': np.log(prior)}
  if scale_prior_with_samps:
    logprior['cocluster'] *= S
  return logprior

def _normalize_logconc(logconc, S):
  # Convert from base 10 to base e.
  logconc /= np.log10(np.e)
  logconc *= S
  logconc = np.maximum(-600, logconc)
  logconc = np.minimum(600,  logconc)
  return logconc

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--concentration', dest='logconc', type=float, default=-2,
    help='log10(alpha) for Chinese restaurant process. The larger this is, the stronger the preference for more clusters.')
  parser.add_argument('--parallel', dest='parallel', type=int, default=1,
    help='Number of tasks to run in parallel')
  parser.add_argument('--prior', type=float, default=0.25,
    help='Pairwise coclustering prior probability. Used only for --model=pairwise or --model=both.')
  parser.add_argument('--model', choices=('pairwise', 'linfreq'), required=True,
    help='Clustering model to use')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  args = parser.parse_args()

  variants = inputparser.load_ssms(args.ssm_fn)
  params = inputparser.load_params(args.params_fn)
  clusters = params['clusters']
  garbage = params.get('garbage', [])
  variants = inputparser.remove_garbage(variants, garbage)

  M = len(variants)
  S = len(list(variants.values())[0]['var_reads'])
  logconc = _normalize_logconc(args.logconc, S)

  if args.model == 'pairwise':
    vids, Z = cluster_pairwise._convert_clustering_to_assignment(clusters)
    logprior = _make_coclust_logprior(args.prior, S)
    mutrel_posterior, mutrel_evidence = pairwise.calc_posterior(variants, logprior, 'mutation', args.parallel)
    assert vids == mutrel_posterior.vids
    log_clust_probs, log_notclust_probs = cluster_pairwise._make_coclust_probs(mutrel_posterior)
    llh = cluster_pairwise._calc_llh(Z, log_clust_probs, log_notclust_probs, logconc)
  elif args.model == 'linfreq':
    vids1, V, T, T_prime, omega = inputparser.load_read_counts(variants)
    vids2, Z = cluster_pairwise._convert_clustering_to_assignment(clusters)
    assert vids1 == vids2

    # Beta distribution prior for phi
    phi_alpha0 = 1.
    phi_beta0 = 1.
    llh = cluster_linfreq._calc_llh(V, T_prime, Z, phi_alpha0, phi_beta0, logconc)
  else:
    raise Exception('Unknown model')

  nlglh = -llh / (M*S*np.log(2))
  print(llh, nlglh)

if __name__ == '__main__':
  main()
