import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import util
import mutrel
import common
from common import Models

np.seterr(divide='raise', invalid='raise')

def fix_rounding_errors(mat):
  # Floating point error means that some entires can slightly exceed 1, even if
  # constituents of average obey the constraint `0 <= entry <= 1`. Ensure this
  # doesn't happen by much.
  #
  # This issue typically arises when performing evaluations, using an "average"
  # matrix created by summing many constituent weighted matrices.
  assert np.allclose(1, mat[mat > 1])
  mat = np.minimum(1, mat)
  assert np.all(0 <= mat)
  return mat

def add_garbage(posterior, garb_svids):
  if len(garb_svids) == 0:
    return posterior
  assert len(set(posterior.vids) & set(garb_svids)) == 0
  new_vids = posterior.vids + garb_svids
  new_posterior = mutrel.init_mutrel(new_vids)
  G = len(garb_svids)
  M = len(new_posterior.vids)

  # Rather than carefully slicing and dicing the array to set it, just use a
  # series of carefully ordered overwrite operations to put it in the correct
  # state.
  new_posterior.rels[:] = 0
  new_posterior.rels[:,:,Models.garbage] = 1
  diag = range(M)
  new_posterior.rels[diag,diag,:] = 0
  new_posterior.rels[diag,diag,Models.cocluster] = 1
  new_posterior.rels[:-G,:-G,:] = posterior.rels

  mutrel.check_posterior_sanity(new_posterior.rels)
  return new_posterior

def save_sorted_mutrel(mrel, mrelfn):
  mrel = mutrel.sort_mutrel_by_vids(mrel)
  mutrel.check_posterior_sanity(mrel.rels)
  np.savez_compressed(mrelfn, rels=mrel.rels, vids=mrel.vids)

def _distinguish_unique_trees(adjms, logweights, clusterings):
  # Uncomment these two lines to short-circuit the "distinguish unique trees"
  # behaviour.
  #M = len(adjms)
  #return (adjms, clusterings, logweights, np.ones(M))
  used_logweights = set()
  uniq_adjms = []
  uniq_clusterings = []
  uniq_logweights = []
  counts = []

  for adjm, logweight, clustering in zip(adjms, logweights, clusterings):
    clustering = tuple((set(C) for C in clustering))
    for idx in range(len(uniq_adjms)):
      # Note that having the same clustering and adjm as a previous sample
      # doesn't imply that the LLH will be close -- PWGS may have assigned
      # different phis between two such samples, leading to a different LLH.
      # This wil not, however, be a problem for Pairtree, since it will always
      # have the same inter-sample LLH in these cases.
      #
      # We use `used_logweights` as a cheap O(1) lookup for whether we've
      # already seen a sampled `(adjm, logweight, clustering)` triplet,
      # avoiding a quadratic scan through the lists when unnecessary. This
      # relies on the logweights of equivalent samples being identical rather
      # than just close, but that's okay -- we expect that two truly identical
      # samples should have the same LLH, as LLH calculation should be
      # deterministic.
      if logweight in used_logweights and \
      uniq_logweights[idx] == logweight and \
      uniq_clusterings[idx] == clustering and \
      np.array_equal(uniq_adjms[idx], adjm):
        counts[idx] += 1
        break
    else:
      uniq_adjms.append(adjm)
      uniq_clusterings.append(clustering)
      uniq_logweights.append(logweight)
      used_logweights.add(logweight)
      counts.append(1)

  return (
    uniq_adjms,
    uniq_clusterings,
    np.array(uniq_logweights),
    np.array(counts),
  )

def make_logweights(llhs, weight_type):
  if weight_type == 'llh':
    return np.copy(llhs)
  elif weight_type == 'uniform':
    return np.zeros(len(llhs))
  else:
    raise Exception('Unknown weight_type=%s' % weight_type)

def make_mutrel_from_trees_and_unique_clusterings(adjms, llhs, clusterings, weight_type):
  '''
  Relative to `make_mutrel_from_trees_and_single_clustering`, this function is
  slower and more memory intensive, but also more flexible. It differs in three
  respects:

  1. It doesn't assume that the user has already computed counts for all unique
  samples -- i.e., it allows duplicate samples and performs the counting
  itself.

  2. It allows unique clusterings for every sample.

  3. It works on parent vectors rather than adjacency matrices.
  '''
  assert len(adjms) == len(llhs) == len(clusterings)
  logweights = make_logweights(llhs, weight_type)
  uniq_adjms, uniq_clusterings, uniq_logweights, counts = _distinguish_unique_trees(adjms, logweights, clusterings)

  # Oftentimes, we will have many samples of the same adjacency matrix paired
  # with the same clustering. This will produce the same mutrel. As computing
  # the mutrel from adjm + clustering is expensive, we want to avoid repeating
  # this unnecessarily. Instead, we just modify the associated weight of the
  # the pairing to reflect this.
  #
  # Weights are stored in log space. To get linear-space weights, we take the
  # softmax of the logweights. Observe that if we have `C` copies of the
  # logweight `W`, we obtain equivalent post-softmax linear-space weights under
  # either of the following two methods:
  #
  # 1. (naive) Represent the associated samples `C` separate times in the softmax
  # 2. (smart) Set `W' = W + log(C)`, as `exp(W') = Cexp(W)`
  weights = util.softmax(uniq_logweights + np.log(counts))

  vids = None
  for adjm, clustering, weight in zip(uniq_adjms, uniq_clusterings, weights):
    mrel = make_mutrel_from_cluster_adj(adjm, clustering)
    if vids is None:
      vids = mrel.vids
      soft_mutrel = np.zeros(mrel.rels.shape)
    else:
      assert mrel.vids == vids
    soft_mutrel += weight * mrel.rels

  soft_mutrel = fix_rounding_errors(soft_mutrel)
  return mutrel.Mutrel(
    vids = vids,
    rels = soft_mutrel,
  )

def make_mutrel_from_trees_and_single_clustering(structs, llhs, counts, clustering):
  weights = util.softmax(llhs + np.log(counts))
  vids = None

  for struct, weight in zip(structs, weights):
    adjm = util.convert_parents_to_adjmatrix(struct)
    crel = make_clustrel_from_cluster_adj(adjm)

    if vids is None:
      vids = crel.vids
      soft_clustrel = np.zeros(crel.rels.shape)
    else:
      assert crel.vids == vids
    soft_clustrel += weight * crel.rels

  soft_clustrel = fix_rounding_errors(soft_clustrel)
  clustrel = mutrel.Mutrel(rels=soft_clustrel, vids=vids)
  mrel = make_mutrel_from_clustrel(clustrel, clustering)
  return mrel

def make_membership_mat(clusters):
  vids = common.sort_vids([vid for C in clusters for vid in C])
  vidmap = {vid: vidx for vidx, vid in enumerate(vids)}
  N = len(vids)
  K = len(clusters)

  # membership[i,j] = 1 iff mutation `i` is in cluster `j`
  membership = np.zeros((N, K))
  for cidx, C in enumerate(clusters):
    members = [vidmap[vid] for vid in C]
    membership[members,cidx] = 1
  return (vids, membership)

def make_mutrel_from_cluster_adj(cluster_adj, clusters):
  '''
  * `M` = # of mutations
  * `K` = # of clusters

  Arguments:
  `cluster_adj`: a `KxK` adjacency matrix, where `cluster_adj[a,b] = 1` iff
  `a = b` or `b` is a child of `a`
  `clusters`: a `K`-length list of lists, forming a partition over the variant IDs `[s0, s1, ..., s(M-1)]`

  Returns:
  an `MxMx5` binary mutation relation tensor
  '''
  clustrel = make_clustrel_from_cluster_adj(cluster_adj)
  mutrel = make_mutrel_from_clustrel(clustrel, clusters)
  return mutrel

def make_clustrel_from_cluster_adj(cluster_adj):
  '''
  * `K` = # of clusters (including empty first cluster)

  Arguments:
  `cluster_adj`: a `KxK` adjacency matrix, where `cluster_adj[a,b] = 1` iff
  `a = b` or `b` is a child of `a`

  Returns:
  a `KxKx5` binary mutation relation tensor
  '''
  K = len(cluster_adj)
  assert cluster_adj.shape == (K, K)
  cluster_anc = common.make_ancestral_from_adj(cluster_adj)
  # In determining A_B relations, don't want to set mutations (i,j), where i
  # and j are in same cluster, to 1.
  assert np.all(1 == cluster_anc[0])
  np.fill_diagonal(cluster_anc, 0)

  clustrel = np.zeros((K, K, len(Models._all)))
  clustrel[:,:,Models.cocluster] = np.eye(K)
  clustrel[:,:,Models.A_B] = cluster_anc
  clustrel[:,:,Models.B_A] = clustrel[:,:,Models.A_B].T

  existing = (Models.cocluster, Models.A_B, Models.B_A)
  already_filled = np.sum(clustrel[:,:,existing], axis=2)
  clustrel[already_filled == 0, Models.diff_branches] = 1

  assert np.array_equal(np.ones((K,K)), np.sum(clustrel, axis=2))
  vids = ['S%s' % idx for idx in range(K - 1)]
  clustrel = mutrel.Mutrel(vids=vids, rels=clustrel[1:,1:])
  mutrel.check_posterior_sanity(clustrel.rels)
  return clustrel

def make_mutrel_from_clustrel(clustrel, clusters):
  mutrel.check_posterior_sanity(clustrel.rels)
  assert len(clusters[0]) == 0
  vids, membership = make_membership_mat(clusters[1:])
  # K: number of non-empty clusters
  M, K = membership.shape

  num_models = len(Models._all)
  mrel = np.zeros((M, M, num_models))
  assert clustrel.rels.shape == (K, K, num_models)

  for modelidx in range(num_models):
    mut_vs_cluster = np.dot(membership, clustrel.rels[:,:,modelidx]) # MxK
    mrel[:,:,modelidx] = np.dot(mut_vs_cluster, membership.T)
  mutrel.check_posterior_sanity(mrel)

  return mutrel.Mutrel(
    vids = vids,
    rels = mrel,
  )
