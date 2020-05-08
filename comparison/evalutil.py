import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import util
import mutrel
import common
from common import Models, NUM_MODELS

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

def make_mutrel_from_trees_and_unique_clusterings(structs, llhs, clusterings):
  '''
  Relative to `make_mutrel_from_trees_and_single_clustering`, this function is
  slower and more memory intensive, but also more flexible. It differs in two
  respects:

  1. It doesn't assume that the user has already computed counts for all unique
  samples -- i.e., it allows duplicate samples.

  2. It allows unique clusterings for every sample.
  '''
  assert len(structs) == len(llhs) == len(clusterings)
  weights = util.softmax(llhs)
  vids = None

  for struct, clustering, weight in zip(structs, clusterings, weights):
    adjm = util.convert_parents_to_adjmatrix(struct)
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
  # Oftentimes, we will have many samples of the same adjacency matrix paired
  # with the same clustering. This will produce the same mutrel. As computing
  # the mutrel from adjm + clustering is expensive, we want to avoid repeating
  # this unnecessarily. Instead, we just modify the associated weight of the
  # the pairing to reflect this.
  #
  # Observe that if we have `C` copies of the LLH `W`, we obtain
  # equivalent post-softmax linear-space weights under either of the following
  # two methods:
  #
  # 1. (naive) Represent the associated samples `C` separate times in the softmax
  # 2. (smart) Set `W' = W + log(C)`, as `exp(W') = Cexp(W)`
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
  mrel = make_mutrel_from_clustrel(clustrel, clusters)
  return mrel

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
  cluster_anc = util.make_ancestral_from_adj(cluster_adj)
  # In determining A_B relations, don't want to set mutations (i,j), where i
  # and j are in same cluster, to 1.
  assert np.all(1 == cluster_anc[0])
  np.fill_diagonal(cluster_anc, 0)

  clustrel = np.zeros((K, K, NUM_MODELS))
  clustrel[:,:,Models.cocluster] = np.eye(K)
  clustrel[:,:,Models.A_B] = cluster_anc
  clustrel[:,:,Models.B_A] = clustrel[:,:,Models.A_B].T

  existing = (Models.cocluster, Models.A_B, Models.B_A)
  already_filled = np.sum(clustrel[:,:,existing], axis=2)
  clustrel[already_filled == 0, Models.diff_branches] = 1

  assert np.array_equal(np.ones((K,K)), np.sum(clustrel, axis=2))
  vids = ['S%s' % (idx + 1) for idx in range(K)]
  clustrel = mutrel.Mutrel(vids=vids, rels=clustrel)
  mutrel.check_posterior_sanity(clustrel.rels)
  return clustrel

def make_mutrel_from_clustrel(clustrel, clusters, check_sanity=True):
  mutrel.check_posterior_sanity(clustrel.rels)
  K = len(clusters)
  assert clustrel.rels.shape == (K, K, NUM_MODELS)

  vids, membership = make_membership_mat(clusters)
  # K: number of non-empty clusters
  M = len(membership)
  assert len(vids) == M
  assert membership.shape == (M, K)

  mrel = np.zeros((M, M, NUM_MODELS))

  for modelidx in range(NUM_MODELS):
    mut_vs_cluster = np.dot(membership, clustrel.rels[:,:,modelidx]) # MxK
    mrel[:,:,modelidx] = np.dot(mut_vs_cluster, membership.T)
  # Disable check to improve performance. Since this is called for each tree
  # (for methods that don't have a fixed clustering), it can be prohibitively
  # slow -- it was consuming >50% of the total runtime for LICHeE's output
  # conversion.
  #mutrel.check_posterior_sanity(mrel)

  return mutrel.Mutrel(
    vids = vids,
    rels = mrel,
  )
