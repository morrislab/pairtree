import plotly
import plotly.graph_objs as go
import numpy as np
import sklearn.cluster
import common

def find_xeno_ranges(sampnames):
  assert not common.is_xeno(sampnames[0])

  last_was_xeno = False
  xeno_range_start = None
  xeno_ranges = []
  for idx, S in enumerate(sampnames[1:]):
    idx += 1
    if common.is_xeno(S) and not last_was_xeno:
      assert xeno_range_start is None
      xeno_range_start = idx
    elif not common.is_xeno(S) and last_was_xeno:
      xeno_ranges.append((xeno_range_start, idx))
      xeno_range_start = None
    last_was_xeno = common.is_xeno(S)

  if xeno_range_start is not None:
    xeno_ranges.append((xeno_range_start, len(sampnames)))
    xeno_range_start = None

  proposed_xenos = len(sampnames) * [False]
  for start, end in xeno_ranges:
    for idx in range(start, end):
      proposed_xenos[idx] = True
  truth_xenos = [common.is_xeno(S) for S in sampnames]
  assert np.array_equal(np.array(proposed_xenos), np.array(truth_xenos))

  return xeno_ranges

def reorder_samples(eta, sampnames):
  xeno_ranges = find_xeno_ranges(sampnames)
  for start, end in xeno_ranges:
    eta, idxs = common.reorder_cols(eta, start, end)
    sampnames = [sampnames[I] for I in idxs]
  return (eta, sampnames)

def hide_unwanted(eta, sampnames):
  def _is_unwanted(samp):
    if 'cns' in samp.lower() or 'spleen' in samp.lower():
      return True
    if samp in ('D2', 'D3', 'R2'):
      return True
    if samp.startswith('Pre-Diagnosis Xeno'):
      return True
    return False
  sampnames, eta = zip(*[(S, col) for S, col in zip(sampnames, eta.T) if not _is_unwanted(S)])
  return (np.array(eta).T, sampnames)

def plot_eta(eta, sampnames, outf):
  eta, sampnames = hide_unwanted(eta, sampnames)
  eta, sampnames = reorder_samples(eta, sampnames)
  short_sampnames = [S.replace('Diagnosis Xeno', 'DX').replace('Relapse Xeno', 'RX')  for S in sampnames]
  widths = np.array([1.0 if not common.is_xeno(S) else 0.4 for S in sampnames])

  # eta: KxS, K = # of clusters, S = # of samples
  traces = [go.Bar(
    name = 'Population %s' % cidx,
    x = short_sampnames,
    y = row,
    text = ['Population %s: %.0f%%' % (cidx, 100*E) for E in row],
    hoverinfo = 'x+text',
    width = widths,
  ) for cidx, row in zip(range(len(eta)), eta)]
  layout = go.Layout(
    barmode = 'stack',
    bargap = 0.1,
  )
  fig = go.Figure(data=traces, layout=layout)

  plot = plotly.offline.plot(
    fig,
    output_type = 'div',
  )
  print(plot, file=outf)
