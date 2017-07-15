import plotly
import plotly.graph_objs as go
import numpy as np
import sklearn.cluster

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
  # between data points, and sklearn will calculate the full tree for us
  # regardless. (Though, with more data points, we might need to force this by
  # passing the `complete_full_tree` to the clusterer constructor.
  agglo = sklearn.cluster.AgglomerativeClustering(n_clusters=3)
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

def reorder_samples(eta, sampnames):
  eta, idxs = reorder_cols(eta)
  sampnames = [sampnames[I] for I in idxs]
  return (eta, sampnames)

def plot_eta(eta, sampnames, outf):
  eta, sampnames = reorder_samples(eta, sampnames)

  # eta: KxS, K = # of clusters, S = # of samples
  etacum = np.cumsum(eta, axis=0)
  traces = [go.Scatter(
    name = 'Population %s' % cidx,
    x = sampnames,
    y = cumrow,
    text = ['Population %s: %.0f%%' % (cidx, 100*E) for E in origrow],
    hoverinfo = 'x+text',
    mode = 'lines',
    fill = 'tonexty'
  ) for cidx, cumrow, origrow in zip(range(len(etacum)), etacum, eta)]
  layout = go.Layout(
    xaxis = {'showgrid': False}
  )
  fig = go.Figure(data=traces, layout=layout)

  plot = plotly.offline.plot(
    fig,
    output_type = 'div',
  )
  print(plot, file=outf)
