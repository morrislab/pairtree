import pandas as pd
import argparse
import plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from collections import defaultdict
import re
import numpy as np
import sys
import os

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
  ('pairtree_tensor', 'Pairs tensor'),
  ('mle_unconstrained', 'MLE lineage frequencies'),
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

def make_pie_traces(methods, complete, total):
  methods = sort_methods(methods)
  complete = np.array([complete[M] for M in methods])
  total = np.array([total[M] for M in methods])
  missing = total - complete

  traces = [go.Pie(
    labels = ('Succeeded', 'Failed'),
    values = (C, D),
    name = M,
    sort = False,
  ) for (C, D, M) in zip(complete, missing, methods)]
  return traces

def make_pie_fig(methods, complete, total, name=None):
  methods = sort_methods(methods)
  complete = np.array([complete[M] for M in methods])
  total = np.array([total[M] for M in methods])
  missing = total - complete

  traces = [go.Pie(
    labels=('Succeeded', 'Failed'),
    values=(C, D),
    name=M,
    textinfo='none',
  ) for (C, D, M) in zip(complete, missing, methods)]

  pie_fig = make_subplots(
    rows = 1,
    cols = len(traces),
    specs = [len(traces)*[{'type': 'domain'}]],
    subplot_titles=methods,
  )
  for idx, T in enumerate(traces):
    pie_fig.add_trace(T, row=1, col=(idx+1))
  if name is not None:
    pie_fig.update_layout(title_text=name)
  return pie_fig

def _make_points(vals):
  methods = sort_methods(vals.keys())
  points = [(HAPPY_METHOD_NAMES.get(M, M), R) for M in methods for R in vals[M]]
  X, Y = zip(*points)
  return (X, Y)

def make_violin_trace(vals, side='both', group=None, name=None, bandwidth=None):
  X, Y = _make_points(vals)
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

  return trace

def make_box_trace(vals, group=None, name=None):
  X, Y = _make_points(vals)
  trace = {
    'type': 'box',
    'x': X,
    'y': Y,
  }
  if group is not None:
    trace['legendgroup'] = group
  if name is not None:
    trace['name'] = name

  return trace

def make_ticks(traces):
  Y = np.array([val for T in traces for val in T['y']])
  minY = np.floor(np.min(Y))
  maxY = np.ceil(np.max(Y))
  N = int(maxY - minY)
  tickvals = np.linspace(minY, maxY, num=(N+1))
  ticktext = 10**tickvals
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

def write_figs(figs, outfn):
  plot = ''
  for idx, fig in enumerate(figs):
    imgfn = os.path.basename(outfn)
    if imgfn.endswith('.html'):
      imgfn = imgfn[:-5]
    imgfn = imgfn.replace('.', '_')
    imgfn = '%s_%s' % (imgfn, idx + 1)

    plot += plotly.offline.plot(
      fig,
      output_type = 'div',
      include_plotlyjs = False,
      config = {
        'showLink': True,
        'toImageButtonOptions': {
          'format': 'svg',
          'width': 750,
          'height': 450,
          'filename': imgfn,
        },
      },
    )

  with open(outfn, 'w') as outf:
    #print('<script type="text/javascript">%s</script>' % read_plotly(), file=outf)
    print( '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>', file=outf)
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
  elif score_type in ('cputime', 'walltime'):
    ytitle = 'Runtime (seconds)'
  else:
    raise Exception('Unknown score type %s' % score_type)
  return ytitle

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--S-threshold', type=int)
  parser.add_argument('--template', default='seaborn')
  parser.add_argument('--max-y')
  parser.add_argument('--score-type', required=True, choices=('mutrel', 'mutphi', 'cputime', 'walltime'))
  parser.add_argument('--baseline')
  parser.add_argument('--bandwidth', type=float)
  parser.add_argument('--plot-type', choices=('box', 'violin'), default='box')
  parser.add_argument('--log-y-axis', action='store_true')
  parser.add_argument('--hide-method', nargs='+')
  parser.add_argument('results_fn')
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  results, methods = load_results(args.results_fn)
  if args.hide_method is not None:
    methods = [M for M in methods if M not in args.hide_method]

  if args.baseline is not None:
    assert args.baseline in results
    for M in methods:
      present = results[M] != MISSING
      results.loc[present,M] -= results.loc[present,args.baseline]

  if args.score_type == 'mutrel':
    set_missing_to(results, methods, 1)
    for M in methods:
      results[M] *= 100

  sides = {
    'lte': 'negative',
    'gt': 'positive',
  }
  names = {
    'lte': '1, 3, or 10 tissue samples',
    'gt': '30 or 100 tissue samples',
  }

  score_traces = []
  completion_bar_traces = []
  completion_pie_figs = []

  if args.S_threshold is not None:
    results = augment(results, 'S')
    parted_by_S = partition(results, methods, key='S')
    parted_by_thresh = partition_by_threshold(parted_by_S, args.S_threshold)
    for group in parted_by_thresh.keys():
      vals = parted_by_thresh[group]
      methods = [M for M in vals.keys()]
      scores = {M: vals[M]['scores'] for M in methods}
      complete = {M: vals[M]['complete'] for M in methods}
      total = {M: vals[M]['total'] for M in methods}

      if args.plot_type == 'box':
        score_trace = make_box_trace(scores, group=group, name=names[group])
      elif args.plot_type == 'violin':
        score_trace = make_violin_trace(scores, side=sides[group], group=group, name=names[group], bandwidth=args.bandwidth)
      else:
        raise Exception('Unknown plot type %s' % args.plot_type)

      score_traces.append(score_trace)
      completion_bar_traces.append(make_bar_trace(methods, complete, total, name=names[group]))
      completion_pie_figs.append(make_pie_fig(methods, complete, total, name=names[group]))
  else:
    scores = {M: np.array(results[M]) for M in methods}
    total = {M: len(scores[M]) for M in methods}
    complete = {M: np.sum(scores[M] != MISSING) for M in methods}
    if args.plot_type == 'box':
      score_trace = make_box_trace(scores)
    elif args.plot_type == 'violin':
      score_trace = make_violin_trace(scores, side='both', bandwidth=args.bandwidth)
    score_traces.append(score_trace)
    completion_bar_traces.append(make_bar_trace(methods, complete, total))
    completion_pie_figs.append(make_pie_fig(methods, complete, total))

  figs = [
    make_fig(
      score_traces,
      args.template,
      make_score_ytitle(args.score_type, args.plot_fn),
      args.max_y,
      log_y_axis = args.log_y_axis,
      layout_options = {
        'boxmode': 'group',
        'violinmode': 'overlay',
        'violingap': 0.0,
        'violingroupgap': 0.0,
      },
    ),

    make_fig(
      completion_bar_traces,
      args.template,
      'Proportion of failed runs',
      max_y = None,
      layout_options = {'barmode': 'group'},
    ),
  ]
  figs += completion_pie_figs
  write_figs(figs, args.plot_fn)

if __name__ == '__main__':
  main()
