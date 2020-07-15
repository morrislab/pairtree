import numpy as np
import numpy.ma as ma
import util

def _fix_eta(eta):
  # Remove non-cancerous population from eta.
  eta = eta[1:]
  eta = eta / np.sum(eta, axis=0)
  assert np.allclose(1, np.sum(eta, axis=0))

  eta = ma.masked_equal(eta, 0)
  return eta

def calc_cdi(eta):
  'Compute the clone diversity index (CDI), which is the entropy of eta.'
  eta = _fix_eta(eta)
  K, S = eta.shape
  H = -ma.sum(eta * ma.log2(eta), axis=0)
  return H

def calc_cmdi(eta, clusters, struct):
  '''Compute the clone and mutation diversity index (CMDI), which is the joint
  entropy of eta and the mutations presnt in a clone (i.e., the mutations
  specific to it as well as the mutations inherited from its ancestors).'''
  K, S = eta.shape

  adj = util.convert_parents_to_adjmatrix(struct)
  anc = util.make_ancestral_from_adj(adj, check_validity=True)
  assert anc.shape == (K, K)

  vids, mutmem = util.make_membership_mat(clusters)
  M = len(vids)
  # Root node has no associated mutations.
  mutmem = np.insert(mutmem, 0, 0, axis=1)
  assert mutmem.shape == (M, K)
  assert np.sum(mutmem) == M
  # `mutanc[i,j] = 1` iff mutation `i` occurred in node `j` or a node ancestral
  # to it.
  mutanc = np.dot(mutmem, anc)
  # `mutanc_cnt[i]` = number of mutations that occurred in clone `i` and all
  # clones ancestral to it.
  mutanc_cnt = np.sum(mutanc, axis=0)

  assert mutanc_cnt[0] == 0 and np.all(mutanc_cnt[1:] > 0)
  M_k = np.repeat(mutanc_cnt[1:][:,np.newaxis], S, axis=1)
  eta = _fix_eta(eta)
  assert eta.shape == M_k.shape

  H_joint = -ma.sum(eta * (ma.log2(eta) - np.log2(M_k)), axis=0)
  assert H_joint.shape == (S,)
  return H_joint
