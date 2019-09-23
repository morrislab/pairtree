import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import util
import mutrel
import common
from common import Models

np.seterr(divide='raise', invalid='raise')

def _fix_rounding_errors(mat):
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

def _make_logweights(llhs, tree_weights):
  if tree_weights == 'llh':
    return np.copy(llhs)
  elif tree_weights == 'uniform':
    return np.zeros(len(llhs))
  else:
    raise Exception('Unknown tree_weights=%s' % tree_weights)

def calc_mutrel_from_trees(adjms, llhs, clusterings, tree_weights):
  logweights = _make_logweights(llhs, tree_weights)
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
    mrel = make_mutrel_tensor_from_cluster_adj(adjm, clustering)
    if vids is None:
      vids = mrel.vids
      soft_mutrel = np.zeros(mrel.rels.shape)
    else:
      assert mrel.vids == vids
    soft_mutrel += weight * mrel.rels

  soft_mutrel = _fix_rounding_errors(soft_mutrel)
  return mutrel.Mutrel(
    vids = vids,
    rels = soft_mutrel,
  )

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

def make_mutrel_tensor_from_cluster_adj(cluster_adj, clusters):
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
  K = len(clusters)
  M = sum([len(clus) for clus in clusters])
  assert cluster_adj.shape == (K, K)

  vids = common.sort_vids([vid for cluster in clusters for vid in cluster])
  vidmap = {vid: vidx for vidx, vid in enumerate(vids)}
  clusters = [[vidmap[vid] for vid in cluster] for cluster in clusters]

  cluster_anc = common.make_ancestral_from_adj(cluster_adj)
  # In determining A_B relations, don't want to set mutaitons (i,j), where i
  # and j are in same cluster, to 1.
  np.fill_diagonal(cluster_anc, 0)
  mrel = np.zeros((M, M, len(Models._all)))

  for k in range(K):
    self_muts = np.array(clusters[k])
    desc_clusters = np.flatnonzero(cluster_anc[k])
    desc_muts = np.array([vidx for cidx in desc_clusters for vidx in clusters[cidx]])

    if len(self_muts) > 0:
      mrel[self_muts[:,None,None], self_muts[None,:,None], Models.cocluster] = 1
    if len(self_muts) > 0 and len(desc_muts) > 0:
      mrel[self_muts[:,None,None], desc_muts[None,:,None], Models.A_B] = 1

  mrel[:,:,Models.B_A] = mrel[:,:,Models.A_B].T
  existing = (Models.cocluster, Models.A_B, Models.B_A)
  already_filled = np.sum(mrel[:,:,existing], axis=2)
  mrel[already_filled == 0,Models.diff_branches] = 1
  assert np.array_equal(np.ones((M,M)), np.sum(mrel, axis=2))

  return Mutrel(vids=vids, rels=mrel)
