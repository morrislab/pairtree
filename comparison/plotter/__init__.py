import numpy as np
import pandas as pd
import plotly.io as pio
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from collections import defaultdict
import sys
import os
import re

MISSING = -1
HAPPY_METHOD_NAMES = [
  ('truth', 'Truth'),
  ('pairtree_handbuilt', 'Pairtree (manually constructed trees)'),
  ('pairtree_single', 'Pairtree (single-chain)'),
  ('pairtree_multi', 'Pairtree'),
  ('citup', 'CITUP'),
  ('lichee', 'LICHeE'),
  ('pastri', 'PASTRI'),
  ('pplus_supervars', 'PhyloPlus'),
  ('pplus_allvars', 'PhyloPlus (unclustered)'),
  ('pwgs_supervars', 'PhyloWGS'),
  ('pwgs_allvars', 'PhyloWGS (unclustered)'),
  ('pairtree_clustrel', 'Pairs tensor'),
  ('mle_unconstrained', 'Unconstrained lineage frequencies'),
]
SORTED_METHODS = [M for M, M_full in HAPPY_METHOD_NAMES]
HAPPY_METHOD_NAMES = {M: M_full for (M, M_full) in HAPPY_METHOD_NAMES}

def augment(results, param_names):
  if len(param_names) < 1:
    return results
  params = defaultdict(list)
  for rid in results.runid:
    for token in re.findall(rf'([{param_names}]\d+)', rid):
      K, V = token[0], token[1:]
      params[K].append(int(V))

  lengths = np.array([len(V) for V in params.values()])
  if not np.all(lengths == len(results)):
    raise Exception('Some params missing from some points')

  for K in params.keys():
    results[K] = params[K]
  return results

def load_results(resultsfn):
  results = pd.read_csv(resultsfn)
  #results = results.drop(['pairtree_clustrel'], axis=1)
  methods = get_method_names(results)

  for M in methods:
    inf_idxs = np.isinf(results[M])
    if np.any(inf_idxs):
      # Sometimes the LLH may be -inf. (Thanks, PASTRI.) Consider this to be a
      # failed run.
      print('%s has %s infs' % (M, np.sum(inf_idxs)), file=sys.stderr)
      results.loc[inf_idxs,M] = MISSING
  return (results, methods)

def get_method_names(results):
  methnames = [K for K in results.keys() if K not in ('runid', 'truth')]
  return methnames

def partition(results, methods, key):
  if key is not None:
    V_vals = sorted(pd.unique(results[key]))
    V_filters = [results[key] == V for V in V_vals]
  else:
    V_vals = [None]
    V_filters = [len(results)*[True]]

  partitioned = {}
  for V, V_filter in zip(V_vals, V_filters):
    partitioned[V] = {}
    V_results = results[V_filter]

    for M in methods:
      R = V_results[M]
      failed = R == MISSING
      R_succ = R[np.logical_not(failed)]
      partitioned[V][M] = {
        'scores': R_succ,
        'complete': len(R_succ),
        'total': len(R),
      }
    if len(partitioned[V]) == 0:
      print('No method has output for %s=%s' % (key, V))
      del partitioned[V]

  return partitioned

def sort_methods(methods):
  methods = set(methods)
  M_sorted = [M for M in SORTED_METHODS if M in methods]
  unspecified = methods - set(M_sorted)
  M_sorted += sorted(unspecified)
  return M_sorted

def make_bar_trace(methods, complete, total, name=None):
  methods = sort_methods(methods)
  complete = np.array([complete[M] for M in methods])
  total = np.array([total[M] for M in methods])
  missing = total - complete
  frac_missing = 100*(missing / total)

  trace = {
    'type': 'bar',
    'x': [HAPPY_METHOD_NAMES.get(M, M) for M in methods],
    'y': frac_missing,
  }
  if name is not None:
    trace['name'] = name
  return trace

def make_pie_fig(methods, complete, total, name=None):
  methods = sort_methods(methods)
  happy_methods = [HAPPY_METHOD_NAMES.get(M, M) for M in methods]
  complete = np.array([complete[M] for M in methods])
  total = np.array([total[M] for M in methods])
  missing = total - complete

  traces = [go.Pie(
    labels = ('Percent succeeded', 'Percent failed'),
    values = (C, D),
    name = HAPPY_METHOD_NAMES.get(M, M),
    textinfo = 'none',
    sort = False,
    marker = {
      'colors': ['rgb(255,255,255)', 'rgb(228,26,28)'],
      'line': {'color': 'rgba(0,0,0,0.5)', 'width': 2},
    },
  ) for (C, D, M) in zip(complete, missing, methods)]

  pie_fig = make_subplots(
    rows = 1,
    cols = len(traces),
    specs = [len(traces)*[{'type': 'domain'}]],
    subplot_titles=happy_methods,
  )
  for idx, T in enumerate(traces):
    pie_fig.add_trace(T, row=1, col=(idx+1))
  if name is not None:
    pie_fig.update_layout(title_text=name)
  return pie_fig

def make_points_for_methods(scores):
  methods = sort_methods(scores.keys())
  points = [(HAPPY_METHOD_NAMES.get(M, M), R) for M in methods for R in scores[M]]
  X, Y = zip(*points)
  return (X, Y)

def make_points_for_single_method(results, method, X_key):
  points = [(row[X_key], row[method]) for idx, row in results.iterrows() if row[method] != MISSING]
  return points

def make_violin_trace(X, Y, side='both', group=None, name=None, bandwidth=None, colour=None):
  trace = {
    'type': 'violin',
    'x': X,
    'y': Y,
    'spanmode': 'hard',
    'side': side,

    #'box': {'visible': True},
    'points': False,
    #'pointpos': -1,
    #'jitter': 0.5,
    #'marker': {'size': 4, 'opacity': 0.7},
    #'opacity': 0.7,
    'bandwidth': bandwidth,
    'width': 1.0,
  }
  if group is not None:
    trace['legendgroup'] = group
  if name is not None:
    trace['name'] = name
  if colour is not None:
    assert colour.startswith('rgb(')
    assert colour.endswith(')')
    triplet = 'rgba(' + colour[4:-1]
    trace['fillcolor'] = triplet + ',0.5)'
    trace['line_color'] = triplet + ',1.0)'

  return trace

def make_box_trace(X, Y, group=None, name=None, colour=None):
  trace = {
    'type': 'box',
    'x': X,
    'y': Y,
  }
  if group is not None:
    trace['legendgroup'] = group
  if name is not None:
    trace['name'] = name
  if colour is not None:
    trace['marker_color'] = colour

  return trace

def make_ticks(traces):
  Y = np.array([val for T in traces for val in T['y']])
  minY = np.floor(np.min(Y))
  maxY = np.ceil(np.max(Y))
  N = int(maxY - minY)
  tickvals = np.linspace(minY, maxY, num=(N+1))
  assert np.allclose(tickvals, tickvals.astype(np.int))
  ticktext = ['10<sup>%.0f</sup>' % T for T in tickvals]
  return (tickvals, ticktext)

def make_fig(traces, template, ytitle, max_y=None, layout_options=None, log_y_axis=False):
  yaxis = {'title': ytitle}

  if log_y_axis:
    for T in traces:
      T['y'] = np.log10(T['y'])
    yaxis['tickmode'] = 'array'
    yaxis['tickvals'], yaxis['ticktext'] = make_ticks(traces)

  if max_y is not None:
    yaxis['range'] = (0, max_y)
  fig = {
    'data': traces,
    'layout': {
      'template': template,
      'yaxis': yaxis,
    },
  }
  if layout_options is not None:
    fig['layout'] = {**fig['layout'], **layout_options}
  return fig

def read_plotly():
  plotly_path = os.path.join(os.path.dirname(__file__), 'plotly', 'plotly.min.js')
  with open(plotly_path) as F:
    return F.read()

def write_figs(figs, outfn, export_dims=None):
  plot = ''
  if export_dims is None:
    export_dims = {}
  else:
    # Duplicate, give that we modify it
    export_dims = dict(export_dims)

  for idx, fig in enumerate(figs):
    imgfn = os.path.basename(outfn)
    if imgfn.endswith('.html'):
      imgfn = imgfn[:-5]
    imgfn = imgfn.replace('.', '_')
    imgfn = '%s_%s' % (imgfn, idx + 1)

    if idx not in export_dims:
      export_dims[idx] = (750, 450)

    #print(pio.to_json(fig))
    plot += pio.to_html(
      fig,
      include_plotlyjs = 'cdn',
      config = {
        'showLink': False,
        'toImageButtonOptions': {
          'format': 'svg',
          'width': export_dims[idx][0],
          'height': export_dims[idx][1],
          'filename': imgfn,
        },
      },
    )

  with open(outfn, 'w') as outf:
    print(plot, file=outf)

def set_missing_to(results, methods, val):
  for M in methods:
    missing = results[M] == MISSING
    results.loc[missing,M] = val
  return results

def partition_by_threshold(parted, threshold):
  def _combine(vals):
    methods = set([M for V in vals for M in parted[V].keys()])
    combined = {M: {
      'scores': np.concatenate([parted[V][M]['scores'] if M in parted[V] else [] for V in vals]),
      'complete': sum([parted[V][M]['complete'] if M in parted[V] else 0 for V in vals]),
      'total': sum([parted[V][M]['total'] if M in parted[V] else 0 for V in vals]),
    } for M in methods}
    return combined

  V_lte = [V for V in parted.keys() if V <= threshold]
  V_gt  = [V for V in parted.keys() if V >  threshold]
  threshed = {
    'lte': _combine(V_lte),
    'gt': _combine(V_gt),
  }
  return threshed

def munge(results, methods, baseline, score_type, plot_type):
  if baseline is not None:
    assert baseline in results
    for M in methods:
      present = results[M] != MISSING
      results.loc[present,M] -= results.loc[present,baseline]

  if score_type == 'mutrel':
    for M in methods:
      present = results[M] != MISSING
      results.loc[present,M] *= 100
    if plot_type == 'violin':
      set_missing_to(results, methods, 1)

def make_score_ytitle(score_type, plot_fn):
  if score_type == 'mutrel':
    if 'steph' in plot_fn:
      ytitle = 'Tree reconstruction differences from expert baseline'
    else:
      ytitle = 'Tree reconstruction error'
    ytitle += '<br>(% relations)'
  elif score_type == 'mutphi':
    if 'steph' in plot_fn:
      ytitle = 'Frequency reconstruction differences from expert baseline'
    else:
      ytitle = 'Frequency reconstruction error'
    ytitle += '<br>(Δbits / mutation / tissue sample)'
  elif score_type in ('mutdistl1', 'mutdistl2'):
    ytitle = 'Mutation frequency distance '
    if score_type.endswith('l1'):
      ytitle += '(L1)'
    else:
      ytitle += '(L2)'
  elif score_type in ('cputime', 'walltime'):
    ytitle = 'Runtime (seconds)'
  else:
    raise Exception('Unknown score type %s' % score_type)
  return ytitle

