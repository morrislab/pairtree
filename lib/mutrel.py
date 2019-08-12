from collections import namedtuple
import numpy as np
import common
from common import Models
import util

Mutrel = namedtuple('Mutrel', ('vids', 'rels'))

def init_mutrel(vids):
  M = len(vids)
  mrel = Mutrel(vids=list(vids), rels=np.nan*np.ones((M, M, len(Models._all))))
  return mrel

def remove_variants_by_vidx(mrel, vidxs):
  # Make set for efficient `in`.
  vidxs = set(vidxs)
  return Mutrel(
    vids = [vid for vidx, vid in enumerate(mrel.vids) if vidx not in vidxs],
    rels = util.remove_rowcol(mrel.rels, vidxs),
  )

def sort_mutrel_by_vids(mrel):
  sorted_vids = common.sort_vids(mrel.vids)
  if mrel.vids == sorted_vids:
    return mrel

  sort_map = {vid: vidx for vidx, vid in enumerate(mrel.vids)}
  order = [sort_map[vid] for vid in sorted_vids]
  assert sorted_vids == [mrel.vids[idx] for idx in order]
  return Mutrel(
    vids = sorted_vids,
    rels = reorder_array(mrel.rels, order),
  )

def check_mutrel_sanity(mrel):
  '''Check properties that should be true of all mutrel arrays.'''
  assert not np.any(np.isnan(mrel))
  for model in ('garbage', 'cocluster', 'diff_branches'):
    # These should be symmetric.
    midx = getattr(Models, model)
    mat = mrel[:,:,midx]
    assert np.allclose(mat, mat.T)

def check_posterior_sanity(posterior):
  check_mutrel_sanity(posterior)
  assert np.all(0 <= posterior) and np.all(posterior <= 1)
  assert np.allclose(1, np.sum(posterior, axis=2))

  diag = range(len(posterior))
  noncocluster = [getattr(Models, M) for M in Models._all if M != 'cocluster']
  for M in noncocluster:
    assert np.allclose(0, posterior[diag,diag,M])
  assert np.allclose(1, posterior[diag,diag,Models.cocluster])

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

def reorder_array(arr, order):
  '''Reorder first two dimensions of `arr` according to `order`. This
  assumes that first two dimensions have same size. Further dimensions are
  left unmodified.'''
  # common.reorder_square_matrix does hierarchical clustering to determine
  # order. This is a much simpler and more elegant function that uses a
  # user-defined order.
  M = len(arr)
  assert arr.ndim >= 2
  assert arr.shape[:2] == (M, M)
  assert set(order) == set(range(M))

  arr = arr[order,:]
  arr = arr[:,order]
  return arr