import numpy as np
import util
from progressbar import progressbar
import common
from common import Models
import clustermaker

def _calc_llh(Z, log_clust_probs, log_notclust_probs, conc):
  uniq_Z  = set(list(Z))
  C = len(uniq_Z)
  #assert uniq_Z == set(range(C))
  N = len(Z)
  cluster_sizes = np.array([np.sum(Z == c) for c in range(C)])
  assert np.sum(cluster_sizes) == N

  llh  = C * np.log(conc)
  llh += np.sum(util.logfactorial(cluster_sizes - 1))
  llh -= np.sum(np.log(conc + np.arange(N)))

  for cidx in range(C):
    members = np.flatnonzero(Z == cidx)
    nonmembers = np.flatnonzero(Z != cidx)
    p_clust_c = log_clust_probs[np.ix_(members, members)]
    p_notclust_c = log_notclust_probs[np.ix_(members, nonmembers)]
    assert p_clust_c.size + p_notclust_c.size == N*len(members)
    llh += np.sum(np.triu(p_clust_c)) + np.sum(np.triu(p_notclust_c))

  return llh

def _do_gibbs_iter(C, Z, log_clust_probs, log_notclust_probs, conc):
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
      cweights[cidx] = np.log(np.sum(members)) + np.sum(log_clust_probs[vidx,members]) + np.sum(log_notclust_probs[vidx,nonmembers])
    # Consider adding a new cluster.
    cweights[C] = np.log(conc) + np.sum(log_notclust_probs[vidx])
    cweights -= np.log(conc + N - 1)

    #print(cweights)
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

def cluster(variants, raw_clusters, conc, iters, clust_prior, parallel):
  logprior = {'garbage': -np.inf}
  if clust_prior is not None:
    logprior['cocluster'] = np.log(clust_prior)

  vids = common.extract_vids(variants)
  assert set(vids) == set([vid for clust in raw_clusters for vid in clust])

  #import resultserializer
  #results = resultserializer.Results('/tmp/results.npz')
  #if results.has_mutrel('clustrel_posterior'):
  #  supervars = clustermaker.make_cluster_supervars(raw_clusters, variants)
  #  clustrel_posterior = results.get_mutrel('clustrel_posterior')
  #else:
  supervars, clustrel_posterior, _, _, _ = clustermaker.use_pre_existing(
    variants,
    logprior,
    parallel,
    raw_clusters,
    [],
  )
  #  results.add_mutrel('clustrel_posterior', clustrel_posterior)
  #  results.save()
  superclusters = clustermaker.make_superclusters(supervars)

  C = len(superclusters)
  assert C == len(raw_clusters)
  svids, Z = _convert_clustering_to_assignment(superclusters) 
  assert svids == clustrel_posterior.vids
  log_clust_probs, log_notclust_probs = _make_coclust_probs(clustrel_posterior)

  clusterings = []
  llhs = []

  with progressbar(total=iters, desc='Clustering variants', unit='iteration', dynamic_ncols=True) as pbar:
    for I in range(iters):
      pbar.update()
      C, Z = _do_gibbs_iter(C, Z, log_clust_probs, log_notclust_probs, conc)
      llh = _calc_llh(Z, log_clust_probs, log_notclust_probs, conc)
      clusterings.append(Z)
      llhs.append(llh)

  raw_clusterings = _convert_svid_assign_to_raw_assign(clusterings, vids, raw_clusters)
  return (vids, raw_clusterings, np.array(llhs))
