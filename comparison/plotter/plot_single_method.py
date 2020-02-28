import argparse
import numpy as np
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import plotter
from plotter import MISSING
import pandas as pd

def make_score_traces(results, method):
  S_vals = sorted(pd.unique(results['S']))
  traces = []

  for S in S_vals:
    points = [(row['K'], row[method]) for idx, row in results.iterrows() \
      if row[method] != MISSING and row['S'] == S]
    if len(points) == 0:
      continue
    points = sorted(points, key = lambda V: V[0])
    X, Y = zip(*points)
    X = ['%s subclones' % K for K in X]
    trace = {
      'type': 'box',
      'x': X,
      'y': Y,
      'legendgroup': str(S),
      'name': '%s samples' % S,
      'boxmean': True,
    }
    traces.append(trace)

  return traces

def make_completion_traces(results, method):
  S_vals = sorted(pd.unique(results['S']))
  K_vals = sorted(pd.unique(results['K']))
  traces = []

  for S in S_vals:
    prop_missing = {}

    for K in K_vals:
      total = len([
        row for idx, row in results.iterrows() \
        if row['K'] == K \
        and row['S'] == S
      ])
      if total == 0:
        continue
      complete = len([
        row for idx, row in results.iterrows() \
        if row[method] != MISSING \
        and row['K'] == K \
        and row['S'] == S
      ])
      missing = total - complete
      prop_missing[K] = missing / total

    if len(prop_missing) == 0:
      continue
    K_sorted = sorted(prop_missing.keys())
    X = ['%s subclones' % K for K in K_sorted]
    Y = [100*prop_missing[k] for k in K_sorted]
    traces.append({
      'type': 'bar',
      'x': X,
      'y': Y,
      'name': '%s samples' % S,
      })
  return traces

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--template', default='seaborn')
  parser.add_argument('--max-y')
  parser.add_argument('--score-type', required=True, choices=('mutrel', 'mutphi', 'mutdistl1', 'mutdistl2'))
  parser.add_argument('--baseline')
  parser.add_argument('results_fn')
  parser.add_argument('method')
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  results, methods = plotter.load_results(args.results_fn)
  plot_type = 'box'
  plotter.munge(results, methods, args.baseline)

  for key in ('K', 'S'):
    results = plotter.augment(results, key)
  score_traces = make_score_traces(results, args.method)
  completion_traces = make_completion_traces(results, args.method)

  figs = {
    f'scores_{args.method}': plotter.make_fig(
      score_traces,
      args.template,
      plotter.make_score_ytitle(args.score_type, args.plot_fn),
      args.max_y,
      log_y_axis = False,
      layout_options = {
        'boxmode': 'group',
        'violinmode': 'overlay',
        'violingap': 0.0,
        'violingroupgap': 0.0,
      },
    ),
    f'completion_{args.method}': plotter.make_fig(
      completion_traces,
      args.template,
      'Failure rate',
      100,
      log_y_axis = False,
      layout_options = {
        'barmode': 'group',
      },
    ),
  }
  plotter.write_figs(figs, args.plot_fn)

if __name__ == '__main__':
  main()
