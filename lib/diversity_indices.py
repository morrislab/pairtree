import numpy as np
import numpy.ma as ma
import util

def _fix_eta(eta):
  assert np.allclose(1, np.sum(eta, axis=0))
  # Remove non-cancerous population from eta.
  eta = eta[1:]

  # prevents divide by zero error if a sample contains only the non-cancerous populuation
  eta = np.divide(eta, np.sum(eta, axis=0), out=np.zeros_like(eta), where=eta!=0.0)
  eta_sum = np.sum(eta, axis=0)

  # ensure that all samples which contain at least one subclone 
  # now have subclone population frequencies which sum to 1
  assert np.allclose(1, eta_sum[eta_sum.nonzero()[0]])

  eta = ma.masked_equal(eta, 0)
  return eta

def calc_cdi(eta):
  '''
  Compute the clone diversity index (CDI), which is the entropy of eta.

  >>> cdi = calc_cdi([[0.5], [0.3], [0.2]])
  >>> np.isclose(cdi[0], 0.9709505944546686)
  True
  '''
  eta = _fix_eta(eta)
  K, S = eta.shape
  H = -ma.sum(eta * ma.log2(eta), axis=0)
  return H

def calc_cmdi(eta, clusters, struct):
  '''
  Compute the clone and mutation diversity index (CMDI), which is the joint
  entropy of eta and the mutations presnt in a clone (i.e., the mutations
  specific to it as well as the mutations inherited from its ancestors).

  >>> eta = np.array([[0.5], [0.3], [0.2]])
  >>> clusters = [['s0', 's1'], ['s2', 's3', 's4']]
  >>> struct = [0, 1]
  >>> cmdi = calc_cmdi(eta, clusters, struct)
  >>> np.isclose(cmdi[0], 2.4997218324096133)
  True
  '''
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

def calc_cadi(eta, struct):
  '''
  Compute the clone and ancestor diversity index (CADI), which is the joint
  entropy of eta and the subclones ancestral to a clone.

  >>> eta = np.array([[0.5], [0.2], [0.2], [0.1]])
  >>> struct = [0, 1, 1]
  >>> cadi = calc_cadi(eta, struct)
  >>> np.isclose(cadi[0], 2.1219280948873624)
  True
  '''
  K, S = eta.shape

  adj = util.convert_parents_to_adjmatrix(struct)
  anc = util.make_ancestral_from_adj(adj, check_validity=True)
  assert anc.shape == (K, K)
  A = np.sum(anc, axis=0) - 1
  A = np.repeat(A[1:][:,np.newaxis], S, axis=1)
  assert np.all(A >= 1)

  eta = _fix_eta(eta)
  assert A.shape == eta.shape

  H_joint = -ma.sum(eta * (ma.log2(eta) - np.log2(A)), axis=0)
  assert H_joint.shape == (S,)
  return H_joint

def calc_sdi(eta, clusters, eta_threshold=0.01):
  '''
  Compute the Shannon diversity index as given in
  https://doi.org/10.1016/j.cell.2019.08.032, which is the entropy of the
  proportion of mutations in each cluster. Here we include a cluster in a
  sample only if the subclone's associated `eta` exceeds `eta_threshold`.

  >>> eta = np.array([[0.5], [0.4], [0.099], [0.001]])
  >>> clusters = [['s0', 's1'], ['s2', 's3', 's4'], ['s5']]
  >>> sdi = calc_sdi(eta, clusters)
  >>> np.isclose(sdi[0], 0.9709505944546686)
  True
  '''
  eta = _fix_eta(eta)
  K, S = eta.shape
  assert len(clusters) == K

  sdi = []
  for sidx in range(S):
    present = eta[:,sidx] >= eta_threshold
    present_clusters = [C for is_present, C in zip(present, clusters) if is_present]
    N_k = np.sum([len(C) for C in present_clusters])
    p_K = np.array([len(C) / N_k for C in present_clusters])
    sdi_s = -np.sum(p_K * np.log2(p_K))
    sdi.append(sdi_s)

  return np.array(sdi)

if __name__ == '__main__':
  import doctest
  doctest.testmod()
