from collections import namedtuple
import evalutil
import numpy as np
import scipy.stats

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import inputparser
import common

Mutphi = namedtuple('Mutphi', ('vids', 'assays', 'logprobs'))

def _calc_logprob(mphi, vids, variants):
  def _extract_arr(K):
    return np.array([variants[V][K] for V in vids])
  vardata = {K: _extract_arr(K) for K in ('var_reads', 'total_reads', 'omega_v')}
  assert len(set([V.shape for V in vardata.values()] + [mphi.shape])) == 1

  log_px = scipy.stats.binom.logpmf(
    vardata['var_reads'],
    vardata['total_reads'],
    vardata['omega_v'] * mphi,
  )
  assert np.all(log_px <= 0)
  return log_px

def calc_mutphi(cluster_phis, llhs, clusterings, weight_trees_by, ssmsfn):
  variants = inputparser.load_ssms(ssmsfn)
  logweights = evalutil._make_logweights(llhs, weight_trees_by)
  weights = evalutil._softmax(logweights)
  assert len(cluster_phis) == len(llhs) == len(clusterings)

  vids = None
  # TODO: make assays meaningful, rather than just always setting it to None.
  assays = None
  logprobs = None

  for (cluster_phi, clustering, weight) in zip(cluster_phis, clusterings, weights):
    cluster_phi = evalutil._fix_rounding_errors(cluster_phi)
    assert np.all(0 <= cluster_phi) and np.all(cluster_phi <= 1)
    V, membership = evalutil.make_membership_mat(clustering)
    mphi = np.dot(membership, cluster_phi)

    if vids is None:
      vids = V
    assert V == vids
    if logprobs is None:
      logprobs = np.zeros(mphi.shape)

    # TODO: should I be doing something like logsumexp?
    logprobs += weight * _calc_logprob(mphi, vids, variants)

  return Mutphi(vids=vids, assays=assays, logprobs=logprobs)

def write_mutphi(mphi, mutphifn):
  # calc_mutphi should have created `mutphi` with sorted vids, but double-check
  # this is true.
  assert list(mphi.vids) == common.sort_vids(mphi.vids)
  np.savez_compressed(mutphifn, logprobs=mphi.logprobs, vids=mphi.vids, assays=mphi.assays)

def load_mutphi(mutphifn):
  results = np.load(mutphifn)
  return Mutphi(vids=results['vids'], assays=results['assays'], logprobs=results['logprobs'])
