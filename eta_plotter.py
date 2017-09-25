import plotly
import plotly.graph_objs as go
import numpy as np
import common
import json

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

def reorder_samples(mat, sampnames):
  xeno_ranges = find_xeno_ranges(sampnames)
  for start, end in xeno_ranges:
    mat, idxs = common.reorder_cols(mat, start, end)
    sampnames = [sampnames[I] for I in idxs]
  return (mat, sampnames)

def hide_unwanted(mat, sampnames):
  def _is_unwanted(samp):
    if 'cns' in samp.lower() or 'spleen' in samp.lower():
      return True
    if samp.startswith('Pre-Diagnosis Xeno'):
      return True
    return False
  sampnames, mat = zip(*[(S, col) for S, col in zip(sampnames, mat.T) if not _is_unwanted(S)])
  return (np.array(mat).T, sampnames)

def rename_samples(sampnames):
  return [S.replace('Diagnosis Xeno', 'DX').replace('Relapse Xeno', 'RX')  for S in sampnames]

def plot_eta(eta, sampnames, outf):
  eta, sampnames = hide_unwanted(eta, sampnames)
  eta, sampnames = reorder_samples(eta, sampnames)
  short_sampnames = rename_samples(sampnames)
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

def write_phi_json(phi, sampnames, jsonfn):
  phi, sampnames = hide_unwanted(phi, sampnames)
  phi, sampnames = reorder_samples(phi, sampnames)
  short_sampnames = rename_samples(sampnames)
  with open(jsonfn, 'w') as outf:
    json.dump({
      'samples': short_sampnames,
      'phi': phi.tolist()
    }, outf)
