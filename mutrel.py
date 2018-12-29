from collections import namedtuple
import numpy as np
import common
from common import Models

Mutrel = namedtuple('Mutrel', ('vids', 'rels'))

def init_mutrel(vids):
  M = len(vids)
  mrel = Mutrel(vids=list(vids), rels=np.nan*np.ones((M, M, len(Models._all))))
  return mrel

def _remove_rowcol(arr, indices):
  '''Remove rows and columns at `indices`.'''
  # Using a mask requires only creating a new array once (and then presumably
  # another array to apply the mask). Calling np.delete() multiple times would
  # creat separate arrays for every element in `indices`.
  shape = list(arr.shape)
  # Must convert to list *before* array, as indices is often a set, and
  # converting a set directly to an array yields a dtype of `object` that
  # contains the set. It's really weird.
  indices = np.array(list(indices))
  # Only check the first two dimensions, as we ignore all others.
  assert np.all(0 <= indices) and np.all(indices < min(shape[:2]))

  for axis in (0, 1):
    arr = np.delete(arr, indices, axis=axis)
    shape[axis] -= len(indices)

  assert np.array_equal(arr.shape, shape)
  return arr

def remove_variants(mrel, vidxs):
  # Make set for efficient `in`.
  vidxs = set(vidxs)
  return Mutrel(
    vids = [vid for vidx, vid in enumerate(mrel.vids) if vidx not in vidxs],
    rels = _remove_rowcol(mrel.rels, vidxs),
  )

def sort_mutrel_by_vids(mrel):
  sorted_vids = common.sort_vids(mrel.vids)
  sort_map = {vid: vidx for vidx, vid in enumerate(sorted_vids)}
  order = [sort_map[vid] for vid in mrel.vids]
  return Mutrel(
    vids = sorted_vids,
    rels = reorder_array(mrel.rels, order),
  )

def add_garbage(posterior, garb_svids):
  assert len(set(posterior.vids) & set(garb_svids)) == 0
  new_vids = posterior.vids + garb_svids
  new_posterior = init_mutrel(new_vids)
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

  check_posterior_sanity(new_posterior.rels)
  return new_posterior

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
  #
  M = len(arr)
  assert arr.ndim >= 2
  assert arr.shape[:2] == (M, M)
  assert set(order) == set(range(M))

  arr = arr[order,:]
  arr = arr[:,order]
  return arr
