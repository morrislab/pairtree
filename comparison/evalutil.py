import numpy as np
from collections import namedtuple

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import mutrel
import common
from common import Models

Mutphi = namedtuple('Mutphi', ('vids', 'phi'))

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

def _softmax(V):
  V = np.copy(V) - np.max(V)
  return np.exp(V) / np.sum(np.exp(V))

def _make_logweights(llhs, tree_weights):
  if tree_weights == 'llh':
    return np.copy(llhs)
  elif tree_weights == 'uniform':
    return np.zeros(len(llhs))
  else:
    raise Exception('Unknown tree_weights=%s' % tree_weights)

def calc_mutphi(cluster_phis, llhs, clusterings, tree_weights, fix_rounding=True):
  logweights = _make_logweights(llhs, tree_weights)
  weights = _softmax(logweights)
  assert len(cluster_phis) == len(llhs) == len(clusterings)

  vids = None

  for (cluster_phi, clustering, weight) in zip(cluster_phis, clusterings, weights):
    if fix_rounding:
      assert np.all(0 <= cluster_phi) and np.all(cluster_phi <= 1)
    V, membership = make_membership_mat(clustering)
    mutphi = np.dot(membership, cluster_phi)
    if vids is None:
      vids = V
      soft_mutphi = np.zeros(mutphi.shape)
    else:
      assert vids == V
    soft_mutphi += weight * mutphi

  if fix_rounding:
    soft_mutphi = _fix_rounding_errors(soft_mutphi)
  return Mutphi(vids=vids, phi=soft_mutphi)

def save_mutphi(mutphi, mutphifn):
  # calc_mutphi should have created `mutphi` with sorted vids, but double-check
  # this is true.
  assert list(mutphi.vids) == common.sort_vids(mutphi.vids)
  np.savez_compressed(mutphifn, phi=mutphi.phi, vids=mutphi.vids)

def load_mutphi(mutphifn):
  results = np.load(mutphifn)
  return Mutphi(phi=results['phi'], vids=results['vids'])

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
  weights = _softmax(uniq_logweights + np.log(counts))

  vids = None
  for adjm, clustering, weight in zip(uniq_adjms, uniq_clusterings, weights):
    mrel = mutrel.make_mutrel_tensor_from_cluster_adj(adjm, clustering)
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
