import argparse

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import plotter
from plotter import MISSING
import plotly.subplots
import pandas as pd
import plotly.express as px

from plotly.figure_factory._violin import calc_stats as calc_violin_stats

def plot(results, methods, template):
  score_traces = []
  failure_traces = []
  K_vals = sorted(pd.unique(results['K']))
  N = len(K_vals)
  M_sorted = list(reversed(plotter.sort_methods(methods)))
  M_happy = {M: plotter.HAPPY_METHOD_NAMES.get(M, M) for M in M_sorted}

  colour_scale = px.colors.qualitative.Plotly
  colours = {M: colour_scale[idx] for idx, M in enumerate(M_sorted)}

  for K in K_vals:
    K_rows = [row for idx, row in results.iterrows() if row['K'] == K]
    missing_fracs = {M: len([row for row in K_rows if row[M] == MISSING]) / len(K_rows) for M in M_sorted}
    points = {M: [row[M] for row in K_rows if row[M] != MISSING] for M in M_sorted}

    score_traces.append([{
      'type': 'box',
      'x': points[M],
      'name': M_happy[M],
      'marker_color': colours[M],
      'orientation': 'h',
      'width': 0.8,
    } for M in M_sorted])
    failure_traces.append({
      'type': 'bar',
      'x': [100*missing_fracs[M] for M in M_sorted],
      'y': [M_happy[M] for M in M_sorted],
      'marker_color': [colours[M] for M in M_sorted],
      'orientation': 'h',
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
    max_x = max([calc_violin_stats(T['x'])['q3'] for T in st if len(T['x']) > 0])
    fig.update_xaxes(range=(-0.5, 1.05*max_x), col=2, row=idx+1)
  fig.update_xaxes(range=(100, 0), col=1)

  fig.update_xaxes(title_text='Failure rate', row=N, col=1)
  fig.update_xaxes(title_text='Score', row=N, col=2)
  fig.update_layout(showlegend=False, template=template)
  return [fig]

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--template', default='seaborn')
  parser.add_argument('--score-type', required=True, choices=('mutrel', 'mutphi', 'mutdistl1', 'mutdistl2', 'cputime', 'walltime'))
  parser.add_argument('--baseline')
  parser.add_argument('--hide-method', nargs='+')
  parser.add_argument('results_fn')
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  results, methods = plotter.load_results(args.results_fn)
  if args.hide_method is not None:
    methods = [M for M in methods if M not in args.hide_method]
  plotter.munge(results, methods, args.baseline, args.score_type, 'box')
  for key in ('K', 'S'):
    results = plotter.augment(results, key)

  figs = plot(results, methods, args.template)
  plotter.write_figs(figs, args.plot_fn)

if __name__ == '__main__':
  main()
