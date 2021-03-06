import numpy as np
import util
import common
from common import Models
from numba import njit

@njit
def _calc_llh(Z, log_clust_probs, log_notclust_probs, logconc):
  log_clust_probs = np.triu(log_clust_probs)
  log_notclust_probs = np.triu(log_notclust_probs)

  uniq_Z  = set(list(Z))
  C = len(uniq_Z)
  #assert uniq_Z == set(range(C))
  N = len(Z)
  cluster_sizes = np.array([np.sum(Z == c) for c in range(C)])
  assert np.sum(cluster_sizes) == N

  llh  = C * logconc
  llh += np.sum(util.logfactorial(cluster_sizes - 1))
  llh -= np.sum(np.log(np.exp(logconc) + np.arange(N)))

  for cidx in range(C):
    members = np.flatnonzero(Z == cidx)
    nonmembers = np.flatnonzero(Z != cidx)
    p_clust_c = log_clust_probs[members][:,members]
    p_notclust_c = log_notclust_probs[members][:,nonmembers]
    assert p_clust_c.size + p_notclust_c.size == N*len(members)
    llh += np.sum(p_clust_c) + np.sum(p_notclust_c)

  return llh

@njit
def _compute_cweights_full(log_clust_probs, log_notclust_probs, Z, vidx, logconc, C):
  mask = np.ones(len(Z), dtype=np.bool_)
  mask[vidx] = 0
  log_clust_probs_excluded = log_clust_probs[mask][:,mask]
  log_notclust_probs_excluded = log_notclust_probs[mask][:,mask]
  Z_excluded = Z[mask]
  llh_excluded = _calc_llh(Z_excluded, log_clust_probs_excluded, log_notclust_probs_excluded, logconc)

  cweights_full = np.full(C + 1, np.nan)
  Z_prime = np.copy(Z)

  for cidx in range(C + 1):
    Z_prime[vidx] = cidx
    llh_included = _calc_llh(Z_prime, log_clust_probs, log_notclust_probs, logconc)
    cweights_full[cidx] = llh_included - llh_excluded

  return cweights_full

@njit
def _do_gibbs_iter(C, Z, log_clust_probs, log_notclust_probs, logconc, check_full_llh=False):
  N = len(Z)
  Z = np.copy(Z)

  for vidx in range(N):
    old_cluster = Z[vidx]
    Z[vidx] = -1
    if not np.any(Z == old_cluster):
      # If `vidx` was the only member, remove this cluster.
      # Do so by moving the last cluster to its index.
      # (These operations are valid even if `highest == old_cluster`.)
      highest = C - 1
      Z[Z == highest] = old_cluster
      C -= 1

    # `cweights`: LLHs of each cluster destination for `vidx`
    # cweights[C] = LLH of adding new cluster
    cweights = np.empty(C + 1)
    # Consider every possible destination.
    for cidx in range(C):
      members = Z == cidx
      nonmembers = Z != cidx
      cweights[cidx] = np.log(np.sum(members)) + np.sum(log_clust_probs[vidx][members]) + np.sum(log_notclust_probs[vidx][nonmembers])
    # Consider adding a new cluster.
    cweights[C] = logconc + np.sum(log_notclust_probs[vidx])
    cweights -= np.log(np.exp(logconc) + N - 1)

    if check_full_llh:
      cweights_full = _compute_cweights_full(log_clust_probs, log_notclust_probs, Z, vidx, logconc, C)
      assert np.all(util.isclose(cweights, cweights_full))

    cprobs = util.softmax(cweights)
    new_cluster = util.sample_multinom(cprobs)
    Z[vidx] = new_cluster
    if new_cluster == C:
      C += 1

  return (C, Z)

def _make_coclust_probs(mutrel):
  delta = 1e-20
  clust_probs    = np.maximum(delta,     mutrel.rels[:,:,Models.cocluster])
  notclust_probs = np.maximum(delta, 1 - mutrel.rels[:,:,Models.cocluster])
  log_clust_probs    = np.log(clust_probs)
  log_notclust_probs = np.log(notclust_probs)

  # These values should not play into the joint probabilities, so just set them
  # to 0.
  np.fill_diagonal(log_clust_probs, 0)
  np.fill_diagonal(log_notclust_probs, 0)
  assert np.allclose(log_clust_probs, log_clust_probs.T)
  assert np.allclose(log_notclust_probs, log_notclust_probs.T)

  return (log_clust_probs, log_notclust_probs)

def _convert_svid_assign_to_raw_assign(svid_assigns, vids, vid_clusters):
  svid_assigns = np.array(svid_assigns)

  A = len(svid_assigns)
  M = len(vids)
  C = len(vid_clusters)
  assert svid_assigns.shape == (A, C)
  assign = np.empty((A, M), dtype=np.int32)

  vidmap = {vids[idx]: idx for idx in range(M)}
  clustmap = [[vidmap[vid] for vid in vid_clusters[idx]] for idx in range(C)]

  for aidx in range(A):
    for cidx in range(C):
      assign[aidx, clustmap[cidx]] = svid_assigns[aidx,cidx]

  return assign

def _convert_clustering_to_assignment(clusters):
  mapping = {vid: cidx for cidx, cluster in enumerate(clusters) for vid in cluster}
  vids = common.sort_vids(mapping.keys())
  assign = np.array([mapping[vid] for vid in vids], dtype=np.int32)
  return (vids, assign)

def cluster(variants, raw_clusters, supervars, superclusters, clustrel_posterior, logconc, iters, seed, progress_queue):
  np.random.seed(seed % 2**32)
  vids = common.extract_vids(variants)
  assert set(vids) == set([vid for clust in raw_clusters for vid in clust])

  C = len(superclusters)
  assert C == len(raw_clusters)
  svids, Z = _convert_clustering_to_assignment(superclusters) 
  assert svids == clustrel_posterior.vids
  log_clust_probs, log_notclust_probs = _make_coclust_probs(clustrel_posterior)

  clusterings = []
  llhs = []

  for I in range(iters):
    if progress_queue is not None:
      progress_queue.put(I)
    C, Z = _do_gibbs_iter(
      C,
      Z,
      log_clust_probs,
      log_notclust_probs,
      logconc,
      check_full_llh = False,
    )
    llh = _calc_llh(Z, log_clust_probs, log_notclust_probs, logconc)
    clusterings.append(Z)
    llhs.append(llh)

  raw_clusterings = _convert_svid_assign_to_raw_assign(clusterings, vids, raw_clusters)
  return (vids, raw_clusterings, np.array(llhs))
