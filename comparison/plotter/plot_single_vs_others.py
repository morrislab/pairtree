import argparse
import numpy as np
import plotly.graph_objs as go

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import plotter
from plotter import MISSING
import plotly

def make_boxes(results, methods, single):
  assert single in methods
  others = plotter.sort_methods(set(methods) - set((single,)))

  traces = []
  for M in others:
    points = [(row['runid'], row[M] - row[single]) for idx, row in results.iterrows() \
      if MISSING not in (row[single], row[M])]
    if len(points) == 0:
      continue
    runids, Y = zip(*points)
    traces.append(go.Box(
      y = Y,
      text = runids,
      name = '%s (%s runs)' % (plotter.HAPPY_METHOD_NAMES.get(M, M), len(points)),
      boxpoints = 'all',
      jitter = 0.3,
      pointpos = 1.8,
    ))

  assert len(traces) > 0
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
  parser.add_argument('single_method')
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  results, methods = plotter.load_results(args.results_fn)
  plot_type = 'box'
  plotter.munge(results, methods, args.baseline, args.score_type, plot_type)
  for key in ('K', 'S'):
    results = plotter.augment(results, key)

  boxes = make_boxes(results, methods, args.single_method)
  figs = {
      f'scores_{args.single_method}_vs_others': plotter.make_fig(
      boxes,
      args.template,
      plotter.make_score_ytitle(args.score_type, args.plot_fn),
      args.max_y,
      log_y_axis = False,
      layout_options = {
      },
    ),
  }
  plotter.write_figs(figs, args.plot_fn)

if __name__ == '__main__':
  main()
