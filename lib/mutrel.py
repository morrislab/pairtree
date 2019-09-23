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
