import numpy as np
from collections import namedtuple

_LOGEPSILON = -30
_EPSILON    = np.exp(_LOGEPSILON)

Variant = namedtuple('Variant', (
  'id',
  'var_reads',
  'ref_reads',
  'total_reads',
  'vaf',
  'omega_v',
))

def sort_vids(vids):
  return sorted(vids, key = lambda V: int(V[1:]))

def extract_vids(variants):
  return sort_vids(variants.keys())

def convert_variant_dict_to_tuple(V):
  return Variant(**{K: V[K] for K in Variant._fields})

# An `IntEnum` would be cleaner, but it creates Numba problems, so use
# `namedtuple` instead.
_ModelChoice = namedtuple('_ModelChoice', (
  'garbage',
  'cocluster',
  'A_B',
  'B_A',
  'diff_branches',
))
Models = _ModelChoice(garbage=0, cocluster=1, A_B=2, B_A=3, diff_branches=4)
NUM_MODELS = len(Models)
ALL_MODELS = _ModelChoice._fields

def ensure_valid_tree(adj):
  # I had several issues with subtle bugs in my tree initialization algorithm
  # creating invalid trees. This function is useful to ensure that `adj`
  # corresponds to a valid tree.
  adj = np.copy(adj)
  K = len(adj)
  assert np.all(np.diag(adj) == 1)
  np.fill_diagonal(adj, 0)
  visited = set()

  stack = [0]
  while len(stack) > 0:
    P = stack.pop()
    assert P not in visited
    visited.add(P)
    C = list(np.flatnonzero(adj[P]))
    stack += C
  assert visited == set(range(K))

def convert_adjlist_to_adjmatrix(adjlist):
  all_children = [child for C in adjlist.values() for child in C]
  root = 0
  assert root not in all_children and root in adjlist.keys()

  N = max(all_children) + 1
  adjm = np.eye(N)

  for parent in adjlist.keys():
    children = adjlist[parent]
    adjm[parent, children] = 1

  return adjm

def convert_adj_matrix_to_json_adjlist(adjm):
  adjm = np.copy(adjm)
  np.fill_diagonal(adjm, 0)
  adjl = {}

  for parent, child in zip(*np.nonzero(adjm)):
    # JSON keys must be strings.
    parent = str(parent)
    # Must convert from NumPy ints to Python ints. Ugh.
    child = int(child)
    if parent not in adjl:
      adjl[parent] = []
    adjl[parent].append(child)

  return adjl

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
  # Avoid importing sklearn unless necessary.
  import warnings
  with warnings.catch_warnings():
    warnings.simplefilter('ignore', category=DeprecationWarning)
    import sklearn.cluster

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
  # Disable reordering for Pieter's data, since he has a defined sample order
  # he needs.
  #return (mat, np.array(range(len(mat.T))))
  mat = mat.T
  mat, idxs = reorder_rows(mat, start, end)
  return (mat.T, idxs)

def reorder_square_matrix(mat):
  # Reorders the rows, then reorders the columns using the same indices. This
  # works for symmetric matrices, which would have the same result whether I
  # order the rows and then columns, or vice versa. It can also be used for
  # non-symmetric matrices, as it doesn't require symmetry, but there the
  # results will obviously be dependent on whether you order the original
  # matrix or its transpose.
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

def debug(*args, **kwargs):
  if hasattr(debug, 'DEBUG') and debug.DEBUG:
    print(*args, **kwargs)
