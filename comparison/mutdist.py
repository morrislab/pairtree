from collections import namedtuple
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import util
import evalutil
import common

Mutdist = namedtuple('Mutdist', ('vids', 'assays', 'dists'))

def _calc_dist(mphi, baseline):
  dist = np.abs(mphi - baseline)
  return dist

def calc_mutdist(cluster_phis, llhs, clusterings, baseline, counts):
  assert len(cluster_phis) == len(llhs) == len(clusterings) == len(counts)
  weights = util.softmax(llhs + np.log(counts))
  assert np.isclose(1, np.sum(weights))

  vids = None
  # TODO: make assays meaningful, rather than just always setting it to None.
  assays = None
  dists = None

  for (cluster_phi, clustering, weight) in zip(cluster_phis, clusterings, weights):
    cluster_phi = evalutil.fix_rounding_errors(cluster_phi)
    assert np.all(0 <= cluster_phi) and np.all(cluster_phi <= 1)
    V, membership = evalutil.make_membership_mat(clustering)
    mphi = np.dot(membership, cluster_phi)

    if vids is None:
      vids = V
    assert V == vids
    if dists is None:
      dists = np.zeros(mphi.shape)

    weighted = weight * _calc_dist(mphi, baseline['phi'])
    assert not np.any(np.isnan(weighted)) and not np.any(np.isinf(weighted))
    dists += weighted

  assert list(vids) == list(baseline['vids'])
  if assays is not None:
    assert list(assays) == list(baseline['assays'])
  return Mutdist(vids=vids, assays=assays, dists=dists)

def write_mutdist(mdist, mutdistfn):
  # calc_mutdist should have created `mutdist` with sorted vids, but double-check
  # this is true.
  assert list(mdist.vids) == common.sort_vids(mdist.vids)
  np.savez_compressed(mutdistfn, dists=mdist.dists, vids=mdist.vids, assays=mdist.assays)

def load_mutdist(mutdistfn):
  results = np.load(mutdistfn, allow_pickle=True)
  return Mutdist(vids=results['vids'], assays=results['assays'], dists=results['dists'])
