import argparse

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import plotter
from plotter import MISSING
import plotly.subplots
import pandas as pd
import plotly.express as px
import numpy as np

from plotly.figure_factory._violin import calc_stats as calc_violin_stats


def _plot_scores(results, methods, method_colours):
  score_traces = []
  failure_traces = []
  K_vals = sorted(pd.unique(results['K']))
  N = len(K_vals)
  M_sorted = list(reversed(plotter.sort_methods(methods)))
  M_happy = {M: plotter.HAPPY_METHOD_NAMES.get(M, M) for M in M_sorted}


  for kidx, K in enumerate(K_vals):
    K_rows = [row for idx, row in results.iterrows() if row['K'] == K]
    missing_fracs = {M: len([row for row in K_rows if row[M] == MISSING]) / len(K_rows) for M in M_sorted}
    points = {M: [row[M] for row in K_rows if row[M] != MISSING] for M in M_sorted}

    score_traces.append([{
      'type': 'box',
      'x': points[M],
      'name': M_happy[M],
      'marker_color': _format_colour(method_colours[M]),
      'orientation': 'h',
      'width': 0.8,
      # Hack: only display the legend for the first trace. Otherwise, it will
      # be duplicated for each trace. This assumes that the methods on each
      # plot and their corresponding colours are invariant; if this is
      # violated, the legend will be wrong.
      'showlegend': kidx == 0,
    } for M in M_sorted])
    failure_traces.append({
      'type': 'bar',
      'x': [missing_fracs[M] for M in M_sorted],
      'y': [M_happy[M] for M in M_sorted],
      'marker': {
        'color': [_format_colour(method_colours[M], 0.5) for M in M_sorted],
        'line': {
          'width': 2,
          'color': [_format_colour(method_colours[M]) for M in M_sorted],
        },
      },
      'orientation': 'h',
      'showlegend': False,
    })

  fig = plotly.subplots.make_subplots(
    rows=N,
    cols=2,
    column_widths=[0.2, 0.8],
    shared_yaxes=True,
    row_titles = ['%s subclones' % K for K in K_vals],
    horizontal_spacing = 0.03,
    vertical_spacing = 0.03,
  )
  for idx, (st, ft) in enumerate(zip(score_traces, failure_traces)):
    fig.add_trace(ft, col=1, row=idx+1)
    for T in st:
      fig.add_trace(T, col=2, row=idx+1)
    # Whiskers run from `d1` (lower) to `d2` (higher).
    max_x = max([calc_violin_stats(T['x'])['q2'] for T in st if len(T['x']) > 0])
    fig.update_xaxes(
      range=(-0.15*max_x, 1.05*max_x),
      col=2,
      row=idx+1
    )
  fig.update_xaxes(range=(1, 0), col=1)

  fig.update_xaxes(
    title_text='Failure rate',
    row=N,
    col=1,
  )
  fig.update_xaxes(
    tickformat='.0%',
    col=1,
  )
  fig.update_xaxes(
    zeroline=True,
    zerolinewidth=2,
    zerolinecolor='rgba(0,0,0,0.5)',
    col=2,
  )
  fig.update_xaxes(
    title_text = 'Frequency reconstruction error<br>(Δbits / mutation / tissue sample)',
    row = N,
    col = 2,
  )
  fig.update_layout(showlegend=True, legend={'traceorder': 'reversed'})
  return fig

def _pluralize(N, unit):
  lbl = '%s %s' % (N, unit)
  if N != 1:
    lbl += 's'
  return lbl

def _plot_failure_rates(results, methods, method_colours, K_vals, S_vals):
  M_sorted = plotter.sort_methods(methods)
  M_happy = {M: plotter.HAPPY_METHOD_NAMES.get(M, M) for M in M_sorted}
  fig = plotly.subplots.make_subplots(
    rows = 1,
    cols = len(K_vals),
    subplot_titles = [_pluralize(K, 'subclone') for K in K_vals],
  )

  for Kidx, K in enumerate(K_vals):
    traces = {}
    for M in M_sorted:
      M_failures = {}
      for S in S_vals:
        KS_rows = [row for idx, row in results.iterrows() if row['S'] == S and row['K'] == K]
        M_failures[S] = len([row for row in KS_rows if row[M] == MISSING]) / len(KS_rows)

      Y = np.array([M_failures[S] for S in S_vals])
      # Don't bother plotting if failure rate is 100% for all `S`.
      #if np.all(Y == 1):
      #  continue
      fig.add_trace({
        'type': 'scatter',
        'x': [_pluralize(S, 'sample') for S in S_vals],
        'y': Y,
        'name': M_happy[M],
        'line': {'color': _format_colour(method_colours[M]), 'width': 4,},
        'marker': {'size': 14},
        'showlegend': Kidx == 0,
      }, row=1, col=Kidx+1)

  fig.update_yaxes(
    showticklabels = False,
  )
  fig.update_yaxes(
    title = 'Failure rate',
    tickformat = '%s%%',
    showticklabels = True,
    col = 1,
  )
  return fig

def _plot_error_rate(results, methods, method_colours, K_vals, S_vals):
  M_sorted = plotter.sort_methods(methods)
  M_happy = {M: plotter.HAPPY_METHOD_NAMES.get(M, M) for M in M_sorted}
  fig = plotly.subplots.make_subplots(
    rows = 1,
    cols = len(K_vals),
    subplot_titles = [_pluralize(K, 'subclone') for K in K_vals],
    shared_yaxes = True,
  )

  for Kidx, K in enumerate(K_vals):
    traces = {}
    for M in M_sorted:
      M_error = {}
      for S in S_vals:
        scores = np.array([row[M] for idx, row in results.iterrows() if row['S'] == S and row['K'] == K and row[M] != MISSING])
        if len(scores) == 0:
          continue
        M_error[S] = np.median(scores)

      S_present = sorted(M_error.keys())
      fig.add_trace({
        'type': 'scatter',
        'x': [_pluralize(S, 'sample') for S in S_present],
        'y': [M_error[S] for S in S_present],
        'name': M_happy[M],
        'line': {'color': _format_colour(method_colours[M]), 'width': 4,},
        'marker': {'size': 14},
        'showlegend': Kidx == 0,
      }, row=1, col=Kidx+1)

  fig.update_yaxes(
    title = 'Median frequency reconstruction error<br>(Δbits / mutation / tissue sample)',
    col = 1,
  )
  return fig

def _format_colour(triplet, opacity=1):
  assert 0 <= opacity <= 1
  return 'rgba(%s,%s,%s,%s)' % (*triplet, opacity)

def _choose_method_colours(methods):
  colour_scale = px.colors.qualitative.Plotly
  method_colours = {M: colour_scale[idx] for idx, M in enumerate(plotter.sort_methods(methods))}

  rgb = {}
  for meth, C in method_colours.items():
    assert len(C) == 7 and C.startswith('#')
    rgb[meth] = [int(C[1:][idx:idx+2], 16) for idx in range(0, 6, 2)]
  return rgb

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--template', default='seaborn')
  parser.add_argument('--score-type', required=True, choices=('mutrel', 'mutphi', 'mutdistl1', 'mutdistl2', 'cputime', 'walltime'))
  parser.add_argument('--baseline')
  parser.add_argument('results_fn')
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  results, methods = plotter.load_results(args.results_fn)
  plotter.munge(results, methods, args.baseline, args.score_type, 'box')
  for key in ('K', 'S'):
    results = plotter.augment(results, key)
  methods = set(methods)
  method_colours = _choose_method_colours(methods)

  figs = [
    _plot_scores(
      results,
      methods - set(('mle_unconstrained',)),
      method_colours,
    ),
    _plot_failure_rates(
      results,
      methods - set(('mle_unconstrained', 'pwgs_supervars', 'lichee')),
      method_colours,
      (3, 10),
      (1, 3, 10),
    ),
    _plot_error_rate(
      results,
      methods - set(('mle_unconstrained', 'citup', 'pastri')),
      method_colours,
      (3, 10, 30),
      (1, 3, 10, 30, 100),
    ),
  ]

  for F in figs:
    if 'layout' not in F:
      F['layout'] = {}
    F['layout']['template'] = args.template
  plotter.write_figs(figs, args.plot_fn, export_dims = {
    0: (700, 850),
    1: (400, 485),
    2: (400, 485),
  })

if __name__ == '__main__':
  main()
