import csv
import numpy as np
import sklearn.cluster
from collections import namedtuple
import vaf_correcter

class Models:
  _all = ('garbage', 'cocluster', 'A_B', 'B_A', 'diff_branches')
for idx, M in enumerate(Models._all):
  setattr(Models, M, idx)

def parse_ssms(sampid, ssmfn):
  vaf = []
  ssm_ids = []
  var_names = []
  variants = {}

  with open(ssmfn) as F:
    reader = csv.DictReader(F, delimiter='\t')
    for row in reader:
      if False and len(variants) >= 3:
        break
      variant = {
        'id': row['id'],
        'name': row['gene'],
        'ref_reads': np.array([float(V) for V in row['a'].split(',')]),
        'total_reads': np.array([float(V) for V in row['d'].split(',')]),
        'mu_v': float(row['mu_v']),
      }
      variant['var_reads'] = variant['total_reads'] - variant['ref_reads']
      variant['vaf'] = variant['var_reads'] / variant['total_reads']
      variant['chrom'], variant['pos'] = variant['name'].split('_')
      variant['pos'] = int(variant['pos'])
      variants[row['id']] = variant

  vaf_correcter.correct_vafs(sampid, variants)
  return variants

def make_ancestral_from_adj(adj):
  K = len(adj)
  Z = np.zeros((K,K))

  def _find_desc(I, vec):
    # Base case: if I have no children, my ancestor vec is just myself.
    if np.sum(vec) == 0:
      return vec
    else:
      children = np.array([_find_desc(idx, adj[idx]) for (idx, val) in enumerate(vec) if idx != I and val == 1])
      self_and_child = vec + np.sum(children, axis=0)
      self_and_child[self_and_child > 1] = 1
      return self_and_child

  for k in range(K):
    # If we know `adj` is topologically sorted, we can reduce the complexity of
    # this -- we would start at leaves and work our way upward, eliminating
    # need for recursive DFS. But since we don't expect `K` to be large, we can
    # write more general version that works for non-toposorted trees.
    Z[k] = _find_desc(k, adj[k])

  return Z

def make_cluster_supervars(clusters, variants):
  cluster_supervars = {}
  svid2svidx = {}

  for cidx, cluster in enumerate(clusters):
    if len(cluster) == 0:
      continue
    cvars = [variants['s%s' % vidx] for vidx in cluster]

    cluster_total_reads = np.array([V['total_reads'] for V in cvars])
    cluster_var_reads = np.array([V['var_reads'] for V in cvars])
    # Correct for sex variants.
    mu_v = np.array([V['mu_v'] for V in cvars])[:,np.newaxis]
    vaf_corrections = np.array([V['vaf_correction'] for V in cvars])[:,np.newaxis]
    cluster_var_reads = np.round(vaf_corrections * (cluster_var_reads / (2*(1 - mu_v))))

    S = {
      'gene': None,
      'id': 'C%s' % cidx,
      'name': 'C%s' % cidx,
      'chrom': None,
      'pos': None,
      'cluster': cidx,
      'mu_v': 0.499,
      'var_reads': np.sum(cluster_var_reads, axis=0),
      'total_reads': np.sum(cluster_total_reads, axis=0),
    }
    S['ref_reads'] = S['total_reads'] - S['var_reads']
    S['vaf'] = S['var_reads'] / S['total_reads']

    svid2svidx[S['id']] = len(cluster_supervars)
    cluster_supervars[S['id']] = S

  svidx2svid = {V: K for K, V in svid2svidx.items()}
  return (cluster_supervars, svid2svidx, svidx2svid)

def agglo_children_to_adjlist(children, nleaves):
  assert len(children) == nleaves - 1
  adjlist = {}
  for idx, C in enumerate(children):
    adjlist[idx + nleaves] = C
  root = nleaves + len(children) - 1
  assert max(adjlist.keys()) == root
  return (adjlist, root)

def dfs(adjlist, root):
  ordered = []
  def _dfs(A, parent):
    if parent not in A:
      ordered.append(parent)
      return
    for child in A[parent]:
      _dfs(A, child)
  _dfs(adjlist, root)
  return np.array(ordered)

def reorder_rows(mat, start=None, end=None):
  N = len(mat)
  if start is None:
    start = 0
  if end is None:
    end = N

  fullidxs = np.array(range(N))
  submat = mat[start:end]

  # n_clusters doesn't matter, as we're only interested in the linkage tree
  # between data points.
  agglo = sklearn.cluster.AgglomerativeClustering(
    n_clusters = len(submat),
    affinity = 'l2',
    linkage = 'average',
    compute_full_tree = True,
  )
  labels = agglo.fit_predict(submat)
  adjlist, root = agglo_children_to_adjlist(agglo.children_, agglo.n_leaves_)
  idxs = dfs(adjlist, root)
  submat = [submat[I] for I in idxs]

  fullidxs[start:end] = idxs + start
  mat = np.vstack((mat[:start], submat, mat[end:]))
  return (mat, fullidxs)

def reorder_cols(mat, start=None, end=None):
  mat = mat.T
  mat, idxs = reorder_rows(mat, start, end)
  return (mat.T, idxs)

def reorder_square_matrix(mat):
  # Reorders the rows, then reorders the columns using the same indices. This
  # works for symmetric matrices. It can also be used for non-symmetric
  # matrices, as it doesn't require symmetry.
  assert mat.shape[0] == mat.shape[1]
  mat, idxs = reorder_rows(mat)
  mat = mat.T[idxs,:].T
  return (mat, idxs)

def is_xeno(samp):
  return 'xeno' in samp.lower()

def extract_patient_samples(variants, sampnames):
  munged = {}
  patient_mask = np.array([not is_xeno(S) for S in sampnames])
  for vid in variants.keys():
    munged[vid] = dict(variants[vid])
    for K in ('total_reads', 'ref_reads', 'var_reads', 'vaf'):
      munged[vid][K] = variants[vid][K][patient_mask]
  variants = munged
  sampnames = [S for (S, is_patient) in zip(sampnames, patient_mask) if is_patient]
  return (variants, sampnames)
