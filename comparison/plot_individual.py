import pandas as pd
import argparse
import plotly
import plotly.figure_factory as ff
from collections import defaultdict
import re
import numpy as np
import sys

HAPPY_METHOD_NAMES = {
  'mle_unconstrained': 'MLE lineage frequencies',
  'pairtree_clustrel': 'Pairwise cluster relations',
  'pairtree_trees_llh': 'Pairtree (automated tree search)',
  'pairtree_handbuilt': 'Pairtree (manually constructed trees)',
  'pastri_trees_llh': 'PASTRI',
  'pwgs_allvars_single_llh': 'PhyloWGS (no clustering enforced)',
  'pwgs_allvars_multi_llh': 'Multi-chain PhyloWGS (no clustering enforced)',
  'pwgs_supervars_single_llh': 'PhyloWGS (enforced clustering)',
  'pwgs_supervars_multi_llh': 'Multi-chain PhyloWGS (enforced clustering)',
  'truth': 'Truth',
}

for name, synonym in (
  (('pwgs_trees_single_llh', 'pwgs_supervars_single_llh')),
  (('pwgs_trees_multi_llh', 'pwgs_supervars_multi_llh')),
):
  HAPPY_METHOD_NAMES[name] = HAPPY_METHOD_NAMES[synonym]

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
      results.loc[inf_idxs,M] = -1
    #missing_idxs = results[M] == -1
    #results.loc[missing_idxs,M] = 0

  return (results, methods)

def get_method_names(results):
  methnames = [K for K in results.keys() if K != 'runid']
  methnames = [M for M in methnames if not M.endswith('uniform')]
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
      failed = R == -1
      R_succ = R[np.logical_not(failed)]
      if len(R_succ) == 0:
        continue
      partitioned[V][M] = {
        'scores': R_succ,
        'frac_complete': len(R_succ) / len(R),
      }
    if len(partitioned[V]) == 0:
      print('No method has output for %s=%s' % (key, V))
      del partitioned[V]

  return partitioned

def make_bar_traces(parted, _make_legend_label):
  traces = []
  for V in parted.keys():
    methods = sorted(parted[V].keys(), key=lambda M: HAPPY_METHOD_NAMES[M])
    trace = {
      'type': 'bar',
      'x': [HAPPY_METHOD_NAMES[M] for M in methods],
      'y': [parted[V][M] for M in methods],
    }
    if V is None:
      assert len(parted) == 1
    else:
      trace['name'] = _make_legend_label(V)
    traces.append(trace)
  return traces

def make_box_traces(parted, _make_legend_label):
  traces = []
  for V in parted.keys():
    methods = sorted(parted[V].keys(), key=lambda M: HAPPY_METHOD_NAMES[M])
    points = [(HAPPY_METHOD_NAMES[M], R) for M in methods for R in parted[V][M]]
    X, Y = zip(*points)
    trace = {
      'type': 'box',
      'boxpoints': False,
      'x': X,
      'y': Y,
      'boxmean': True
    }
    if V is None:
      assert len(parted) == 1
    else:
      trace['name'] = _make_legend_label(V)
    traces.append(trace)
  return traces

def make_fig(traces, template, ytitle, max_y=None, layout_options=None):
  yaxis = {'title': ytitle}
  if max_y is not None:
    yaxis['range'] = (0, max_y)
  fig = {
    'data': traces,
    'layout': {
      'template': template,
      'yaxis': yaxis,
    }
  }
  if layout_options is not None:
    fig['layout'] = {**fig['layout'], **layout_options}
  return fig

def write_figs(figs, outfn):
  plot = ''
  for fig in figs:
    plot += plotly.offline.plot(
      fig,
      output_type = 'div',
      include_plotlyjs = 'cdn',
    )
  with open(outfn, 'w') as outf:
    print(plot, file=outf)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--partition-by-samples', action='store_true')
  parser.add_argument('--template', default='seaborn')
  parser.add_argument('--max-y')
  parser.add_argument('results_fn')
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  results, methods = load_results(args.results_fn)
  if args.partition_by_samples:
    results = augment(results, 'S')
    parted = partition(results, methods, key='S')
  else:
    parted = partition(results, methods, key=None)

  figs = []
  results_score = {}
  results_frac_complete = {}
  for V in parted:
    results_score[V]         = {M: parted[V][M]['scores']        for M in parted[V]}
    results_frac_complete[V] = {M: parted[V][M]['frac_complete'] for M in parted[V]}

  def _make_legend_label(V):
    suffix = 'sample' if V == 1 else 'samples'
    return '%s %s' % (V, suffix)
  box_traces = make_box_traces(results_score, _make_legend_label)
  bar_traces = make_bar_traces(results_frac_complete, _make_legend_label)

  figs = [
    make_fig(
      box_traces,
      args.template,
      'Error (bits)',
      args.max_y,
      {'boxmode': 'group'},
    ),
    make_fig(
      bar_traces,
      args.template,
      'Proportion of successful runs',
      args.max_y,
      {'barmode': 'group'},
    ),
  ]
  write_figs(figs, args.plot_fn)

if __name__ == '__main__':
  main()
