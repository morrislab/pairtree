from collections import namedtuple
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import util
import evalutil
import common
import mutstat

def _calc_dist(mphi, baseline):
  dist = np.abs(mphi - baseline)
  return dist

def calc_mutdist(cluster_phis, llhs, clusterings, baseline, counts):
  assert len(cluster_phis) == len(llhs) == len(clusterings) == len(counts)
  weights = util.softmax(llhs + np.log(counts))
  assert np.isclose(1, np.sum(weights))
  baseline_phis = baseline.stats

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

    weighted = weight * _calc_dist(mphi, baseline_phis)
    assert not np.any(np.isnan(weighted)) and not np.any(np.isinf(weighted))
    dists += weighted

  assert list(vids) == list(baseline.vids)
  if assays is not None:
    assert list(assays) == list(baseline.assays)
  return mutstat.Mutstat(vids=vids, assays=assays, stats=dists)
