import argparse

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import plotter
from plotter import MISSING
import plotly.subplots
import pandas as pd
import numpy as np

from plotly.figure_factory._violin import calc_stats as calc_violin_stats

AXIS_TITLES = {}

def _make_axis_titles():
  AXIS_TITLES['mutphi'] = {
    'all_scores': 'VAF reconstruction loss<br>(bits / mutation / tissue sample)',
    'S_comparison': 'Median VAF reconstruction loss<br>(bits / mutation / tissue sample)',
    'single_vs_others': 'VAF reconstruction loss<br>(bits / mutation / tissue sample)<br>relative to %s',
  }

  AXIS_TITLES['mutrel'] = {
    'all_scores': 'Relationship reconstruction error<br>(bits / mutation pair)',
    'S_comparison': 'Median relationship reconstruction error<br>(bits / mutation pair)',
    'single_vs_others': 'Relationship reconstruction loss<br>(bits / mutation pair)r<br>relative to %s',
  }

  AXIS_TITLES['walltime'] = {
    'all_scores': 'Wall-clock time (seconds)',
    'S_comparison': 'Median wall-clock time (seconds)',
    'single_vs_others': 'Extra wall-clock time (seconds)<br>relative to %s',
  }

  AXIS_TITLES['cputime'] = {
    'all_scores': 'CPU time (seconds)',
    'S_comparison': 'Median CPU time (seconds)',
    'single_vs_others': 'Extra CPU time (seconds)<br>relative to %s',
  }

  for lpdist in ('l1', 'l2'):
    AXIS_TITLES[f'mutdist{lpdist}'] = {
      'all_scores': f'Relationship reconstruction error<br>({lpdist.upper()} distance / mutation / tissue sample)',
      'S_comparison': f'Median relationship reconstruction error<br>{lpdist.upper()} distance / mutation / tissue sample)',
      'single_vs_others': f'Relationship reconstruction error<br>({lpdist.upper()} distance / mutation / tissue sample)<br>relative to %s',
    }

def _plot_single_vs_others(results, single, methods, method_colours, score_type):
  assert single in methods
  others = plotter.sort_methods(set(methods) - set((single,)))

  traces = []
  for M in others:
    points = [(row['runid'], row[M] - row[single]) for idx, row in results.iterrows() \
      if MISSING not in (row[single], row[M])]
    if len(points) == 0:
      continue
    runids, Y = zip(*points)
    traces.append({
      'type': 'box',
      'y': Y,
      'text': runids,
      'name': plotter.HAPPY_METHOD_NAMES.get(M, M),
      #'spanmode': 'hard',
      #'bandwidth': 0.07,
      'boxpoints': 'all',
      'marker': {
        'outliercolor': plotter.format_colour(method_colours[M], 0.5),
        'color': plotter.format_colour(method_colours[M], 0.5),
      },
      'line': {
        'color': plotter.format_colour(method_colours[M]),
      },
      'jitter': 0.6,
      'pointpos': 2.0,
    })
  fig = {
    'data': traces,
    'layout': {
      'yaxis': {
        'title': AXIS_TITLES[score_type]['single_vs_others'] % plotter.HAPPY_METHOD_NAMES.get(single, single),
        'zeroline': True,
        'zerolinewidth': 1,
        'zerolinecolor': 'rgba(0,0,0,0.3)',
      },
    },
  }

  return fig

def _plot_scores(results, methods, method_colours, score_type, use_same_x_limit=True):
  score_traces = []
  success_traces = []
  K_vals = sorted(pd.unique(results['K']))
  N = len(K_vals)
  M_sorted = list(reversed(plotter.sort_methods(methods)))
  M_happy = {M: plotter.HAPPY_METHOD_NAMES.get(M, M) for M in M_sorted}

  for kidx, K in enumerate(K_vals):
    K_rows = [row for idx, row in results.iterrows() if row['K'] == K]
    missing_fracs = {M: len([row for row in K_rows if row[M] == MISSING]) / len(K_rows) for M in M_sorted}
    points = {M: [row[M] for row in K_rows if row[M] != MISSING] for M in M_sorted}
    runids = {M: [row['runid'] for row in K_rows if row[M] != MISSING] for M in M_sorted}

    score_traces.append([{
      'type': 'box',
      'boxmean': False,
      'x': points[M],
      'text': runids[M],
      'name': M_happy[M],
      'marker_color': plotter.format_colour(method_colours[M]),
      'orientation': 'h',
      'width': 0.8,
      # Hack: only display the legend for the first trace. Otherwise, it will
      # be duplicated for each trace. This assumes that the methods on each
      # plot and their corresponding colours are invariant; if this is
      # violated, the legend will be wrong.
      'showlegend': kidx == 0,
    } for M in M_sorted])
    success_traces.append({
      'type': 'bar',
      'x': [1 - missing_fracs[M] for M in M_sorted],
      'y': [M_happy[M] for M in M_sorted],
      'marker': {
        'color': [plotter.format_colour(method_colours[M], 0.5) for M in M_sorted],
        'line': {
          'width': 2,
          'color': [plotter.format_colour(method_colours[M]) for M in M_sorted],
        },
      },
      'orientation': 'h',
      'showlegend': False,
    })

  fig = plotly.subplots.make_subplots(
    rows=N,
    cols=2,
    column_widths=[0.2, 0.8],
    shared_xaxes = True,
    shared_yaxes = True,
    row_titles = ['%s subclones' % K for K in K_vals],
    horizontal_spacing = 0.03,
    vertical_spacing = 0.03,
  )
  for idx, (st, ft) in enumerate(zip(score_traces, success_traces)):
    fig.add_trace(ft, col=1, row=idx+1)
    for T in st:
      fig.add_trace(T, col=2, row=idx+1)
  fig.update_xaxes(range=(1, 0), col=1)

  # Whiskers run from `d1` (lower) to `d2` (higher).
  max_x = np.array([max([calc_violin_stats(T['x'])['q2'] for T in st if len(T['x']) > 0]) for st in score_traces])
  min_x = np.array([min([np.min(T['x']) for T in st if len(T['x']) > 0]) for st in score_traces])
  if use_same_x_limit:
    max_x[:] = np.max(max_x)
    min_x[:] = np.min(min_x)
  for idx, (A, B) in enumerate(zip(min_x, max_x)):
    assert B >= A
    rng = B - A
    fig.update_xaxes(
      range=(A - 0.05*rng, B + 0.05*rng),
      col=2,
      row=idx+1
    )

  fig.update_xaxes(
    title_text='Success rate',
    row=N,
    col=1,
  )
  fig.update_xaxes(
    tickformat='.0%',
    col=1,
  )
  fig.update_xaxes(
    zeroline=True,
    zerolinewidth=1,
    zerolinecolor='rgba(0,0,0,0.3)',
    col=2,
  )
  fig.update_xaxes(
    title_text = AXIS_TITLES[score_type]['all_scores'],
    row = N,
    col = 2,
  )

  fig.update_layout(showlegend=True, legend={'traceorder': 'reversed'})
  return fig

def _plot_success_rates(results, methods, method_colours, K_vals, S_vals):
  M_sorted = plotter.sort_methods(methods)
  M_happy = {M: plotter.HAPPY_METHOD_NAMES.get(M, M) for M in M_sorted}
  fig = plotly.subplots.make_subplots(
    rows = 1,
    cols = len(K_vals),
    subplot_titles = [plotter.pluralize(K, 'subclone') for K in K_vals],
    x_title = 'Tissue samples',
  )

  for Kidx, K in enumerate(K_vals):
    for M in M_sorted:
      M_successes = {}
      for S in S_vals:
        KS_rows = [row for idx, row in results.iterrows() if row['S'] == S and row['K'] == K]
        M_successes[S] = len([row for row in KS_rows if row[M] != MISSING]) / len(KS_rows)

      Y = np.array([M_successes[S] for S in S_vals])
      fig.add_trace({
        'type': 'scatter',
        'x': [str(S) for S in S_vals],
        'y': Y,
        'name': M_happy[M],
        'line': {'color': plotter.format_colour(method_colours[M]), 'width': 4,},
        'marker': {'size': 14},
        'showlegend': Kidx == 0,
      }, row=1, col=Kidx+1)

  fig.update_xaxes(
    tickangle = 0,
    type = 'category',
  )
  fig.update_yaxes(
    showticklabels = False,
  )
  fig.update_yaxes(
    title = 'Success rate',
    tickformat = '%s%%',
    showticklabels = True,
    col = 1,
  )
  return fig

def _plot_S_comparison(results, methods, method_colours, K_vals, S_vals, score_type):
  M_sorted = plotter.sort_methods(methods)
  M_happy = {M: plotter.HAPPY_METHOD_NAMES.get(M, M) for M in M_sorted}
  fig = plotly.subplots.make_subplots(
    rows = 1,
    cols = len(K_vals),
    subplot_titles = [plotter.pluralize(K, 'subclone') for K in K_vals],
    x_title = 'Tissue samples',
    shared_yaxes = True,
  )

  for Kidx, K in enumerate(K_vals):
    for M in M_sorted:
      M_upper = {}
      M_error = {}
      M_lower = {}

      for S in S_vals:
        scores = np.array([row[M] for idx, row in results.iterrows() if row['S'] == S and row['K'] == K and row[M] != MISSING])
        if len(scores) == 0:
          continue
        M_lower[S] = np.quantile(scores, 0.25)
        M_upper[S] = np.quantile(scores, 0.75)
        M_error[S] = np.median(scores)

      S_present = sorted(M_error.keys())
      fig.add_trace({
        'type': 'scatter',
        'x': [str(S) for S in S_present],
        'y': [M_error[S] for S in S_present],
        'error_y': {
          'type': 'data',
          'symmetric': False,
          'array': [M_upper[S] - M_error[S] for S in S_present],
          'arrayminus': [M_error[S] - M_lower[S] for S in S_present],
        },
        'name': M_happy[M],
        'line': {'color': plotter.format_colour(method_colours[M]), 'width': 4,},
        'marker': {'size': 14},
        'showlegend': Kidx == 0,
      }, row=1, col=Kidx+1)

  fig.update_xaxes(
    tickangle = 45,
    type = 'category',
  )
  fig.update_yaxes(
    title = AXIS_TITLES[score_type]['S_comparison'],
    col = 1,
  )
  return fig

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--template', default='seaborn')
  parser.add_argument('--score-type', required=True, choices=('mutrel', 'mutphi', 'mutdistl1', 'mutdistl2', 'cputime', 'walltime'))
  parser.add_argument('--baseline')
  parser.add_argument('--hide-methods', default='')
  parser.add_argument('results_fn')
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  _make_axis_titles()

  results, methods = plotter.load_results(args.results_fn)
  plotter.munge(results, methods, args.baseline)
  for key in ('K', 'S'):
    results = plotter.augment(results, key)
  methods = set(methods)
  methods -= set(args.hide_methods.split(','))
  method_colours = plotter.choose_method_colours(methods)

  figs = {
    'scores': _plot_scores(
      results,
      methods - set(('mle_unconstrained',)),
      method_colours,
      args.score_type
    ),
    'success_rate': _plot_success_rates(
      results,
      methods - set(('mle_unconstrained', 'pwgs_supervars', 'lichee', 'pairtree_clustrel')),
      method_colours,
      (3, 10),
      (1, 3, 10),
    ),
    'S_comparison': _plot_S_comparison(
      results,
      methods - set(('mle_unconstrained', 'citup', 'pastri')),
      method_colours,
      (3, 10, 30),
      (1, 3, 10, 30, 100),
      args.score_type
    ),
  }

  export_dims = {
    'success_rate': (400, 485),
  }
  export_dims['scores'] = (700, 850)
  export_dims['S_comparison'] = (500, 485)

  for M in plotter.sort_methods(methods):
    if M == 'mle_unconstrained':
      continue
    name = f'{M}_vs_others'
    figs[name] = _plot_single_vs_others(
      results,
      M,
      methods,# - set(('mle_unconstrained',)),
      method_colours,
      args.score_type
    )
    export_dims[name] = (600, 485)

  for F in figs.values():
    if 'layout' not in F:
      F['layout'] = {}
    F['layout']['template'] = args.template
  plotter.write_figs(figs, args.plot_fn, export_dims)

if __name__ == '__main__':
  main()
