import numpy as np
from numba import njit

import pairwise
import common
from common import Models, debug
from mutrel import Mutrel
import util

def _check_clusters(variants, clusters, garbage):
  for C in clusters:
      assert len(C) > 0

  vids = common.extract_vids(variants)
  clustered = [child for C in clusters for child in C]
  garbage = set(garbage)
  clustered = set(clustered)
  assert len(clustered & garbage) == 0
  assert set(vids) == (clustered | garbage)

def use_pre_existing(variants, logprior, parallel, clusters, garbage):
  supervars = make_cluster_supervars(clusters, variants)
  clust_posterior, clust_evidence = pairwise.calc_posterior(supervars, logprior, rel_type='supervariant', parallel=parallel)
  _check_clusters(variants, clusters, garbage)
  return (supervars, clust_posterior, clust_evidence, clusters, garbage)

# This code is currently unused. Perhaps I can implement a garbage-detection
# algorithm in the future using it.
def _discard_garbage(clusters, mutrel_posterior, mutrel_evidence):
  garbage = []

  while True:
    M = len(clusters)
    assert len(mutrel_posterior.vids) == M
    assert mutrel_posterior.rels.shape == (M, M, NUM_MODELS)

    garbage_pairs = np.argmax(mutrel_posterior.rels, axis=2) == Models.garbage
    if not np.any(garbage_pairs):
      break
    num_garbage = np.sum(garbage_pairs, axis=1)
    most_garbage = np.argmax(num_garbage)
    debug('most_garbage=', mutrel_posterior.vids[most_garbage], 'cluster=', clusters[most_garbage], 'garb_post=', max(mutrel_posterior.rels[most_garbage,:,Models.garbage]))

    garbage += clusters[most_garbage]
    del clusters[most_garbage]
    mutrel_posterior = pairwise.remove_variants(mutrel_posterior, (most_garbage,))
    mutrel_evidence = pairwise.remove_variants(mutrel_evidence, (most_garbage,))

  return (clusters, garbage, mutrel_posterior, mutrel_evidence)

def _make_supervar(name, variants):
  assert len(variants) > 0
  N = np.array([var['total_reads'] for var in variants])
  M, S = N.shape
  V = np.array([var['var_reads'] for var in variants])
  omega_v = np.array([var['omega_v'] for var in variants])

  # In creating a supervariant, we must rescale either `V` or `N` (the variant
  # or total read counts, respectively), as we're fixing the supervariant's
  # `omega_v`. We want the ratio `V / N*omega_v` to remain constant, as this
  # represents the "fraction of informative reads that are variant reads". That
  # is, `N * omega_v` indicates the number of reads that *could* be variant,
  # given what we know of the locus' copy-number state.
  #
  # As we will fix `omega_v` at 0.5 for the rescaled variant, and we want `V /
  # N*omega_v` to remain constant, we can either change `V` or `N` to maintain
  # the relationship `V / N*omega_v = V_hat / N_hat*omega_v`. We can solve this
  # equation two ways, depending on whether we modify `V` or `N`:
  #
  # 1. V_hat = V / 2*omega_v
  # 2. N_hat = 2*N*omega_v
  #
  # The advantage to approach 2 is that, when `omega_v` is small, it preserves
  # the binomial variance. We have `var = np(1 - p) =~ np = E[X]` when `p` is
  # small, meaning the variance scales with the mean when `p` is small. To
  # maintain the variance, we should try to maintain `np`, and so since we're
  # making `p` bigger, we should make `N` smaller, rather than changing `V`.
  #
  # Also, this means we can avoid any weird hackery when `omega_v = 0`, since
  # we don't have to add edge case checks to avoid division by zero.
  #
  # Another argument: since we don't have a proper sequencing noise model, a
  # small number of variant reads can be assumed to be noise regardless of what
  # `omega_v` is. If `omega_v` is small and we observe a couple variant reads,
  # we can assume those are solely noise. So, we shouldn't rescale `V` to be
  # really large, which is what we formerly did under solution 1.
  N_hat = 2*N*omega_v
  V_hat = np.minimum(V, N_hat)
  omega_v_hat = 0.5 * np.ones(S)

  svar = {
    'id': name,
    'name': name,
    'chrom': None,
    'pos': None,
    'omega_v': omega_v_hat,
    'var_reads': np.round(np.sum(V_hat, axis=0)).astype(np.int),
    'total_reads': np.round(np.sum(N_hat, axis=0)).astype(np.int),
  }
  svar['ref_reads'] = svar['total_reads'] - svar['var_reads']
  svar['vaf'] = svar['var_reads'] / svar['total_reads']

  return svar

def make_cluster_supervars(clusters, variants):
  supervars = {}

  for cidx, cluster in enumerate(clusters):
    assert len(cluster) > 0
    cvars = [variants[vid] for vid in cluster]
    S_name = 'S%s' % len(supervars)
    supervars[S_name] = _make_supervar(S_name, cvars)

  return supervars

def make_superclusters(supervars):
  N = len(supervars)
  svids = common.extract_vids(supervars)
  superclusters = [[S] for S in svids]
  return superclusters
