import argparse
import numpy as np
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import plotter
from plotter import MISSING
import pandas as pd

def make_traces(results, method):
  S_vals = sorted(pd.unique(results['S']))
  traces = []

  for S in S_vals:
    points = [(row['K'], row[method]) for idx, row in results.iterrows() \
      if row[method] != MISSING and row['S'] == S]
    points = sorted(points, key = lambda V: V[0])
    X, Y = zip(*points)
    X = ['%s subclones' % K for K in X]
    trace = plotter.make_box_trace(X, Y, group=str(S), name='%s samples' % S)
    traces.append(trace)

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
  plotter.munge(results, methods, args.baseline, args.score_type, plot_type)

  for key in ('K', 'S'):
    results = plotter.augment(results, key)
  traces = make_traces(results, args.method)

  figs = [
    plotter.make_fig(
      traces,
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
  ]
  plotter.write_figs(figs, args.plot_fn)

if __name__ == '__main__':
  main()
