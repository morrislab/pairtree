import argparse

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import plotter
from plotter import MISSING
import numpy as np

from plotly.figure_factory._violin import calc_stats as calc_violin_stats

def _plot_single_vs_others(results, single, methods, method_colours):
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
        'title': 'Frequency reconstruction error<br>relative to %s' % plotter.HAPPY_METHOD_NAMES.get(M, M),
        'zeroline': True,
        'zerolinewidth': 2,
        'zerolinecolor': 'rgba(0,0,0,0.5)',
      },
    },
  }
  return fig

def _plot_scores(results, methods, method_colours):
  M_sorted = list(plotter.sort_methods(methods))
  M_happy = {M: plotter.HAPPY_METHOD_NAMES.get(M, M) for M in M_sorted}

  points = {}
  runids = {}

  for M in M_sorted:
    vals = np.array([row[M] for idx, row in results.iterrows()])
    assert not np.any(vals == MISSING)
    points[M] = vals
    runids[M] = [row['runid'] for idx, row in results.iterrows()]

  opacity = 0.5
  score_traces = [{
    'type': 'box',
    'y': points[M],
    'text': runids[M],
    'name': M_happy[M],
    'orientation': 'v',
    'boxpoints': 'all',
    'marker': {
      'outliercolor': plotter.format_colour(method_colours[M], opacity),
      'color': plotter.format_colour(method_colours[M], opacity),
      'size': 12,
    },
    'line': {
      'color': plotter.format_colour(method_colours[M]),
    },
    'jitter': 0.6,
    'pointpos': 2.0,
  } for M in points]

  max_x = max([calc_violin_stats(T['y'])['q2'] for T in score_traces if len(T['y']) > 0])
  fig = {
    'data': score_traces,
    'layout': {
      'yaxis': dict(
        zeroline = True,
        zerolinewidth = 2,
        zerolinecolor = 'rgba(0,0,0,0.5)',
        range = (-0.15*max_x, 1.05*max_x),
        title_text = 'Frequency reconstruction error<br>(Δbits / mutation / tissue sample)',
      ),
      'showlegend': False,
    },
  }
  return fig

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--template', default='seaborn')
  parser.add_argument('--score-type', required=True, choices=('mutphi', 'mutdistl1', 'mutdistl2'))
  parser.add_argument('--hide-methods')
  parser.add_argument('--baseline')
  parser.add_argument('results_fn')
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  results, methods = plotter.load_results(args.results_fn)
  plotter.munge(results, methods, args.baseline)
  methods = set(methods)
  methods -= set(args.hide_methods.split(','))
  method_colours = plotter.choose_method_colours(methods)

  figs = {
    'scores': _plot_scores(
      results,
      methods - set(('mle_unconstrained',)),
      method_colours,
    ),
  }
  export_dims = {
    'scores': (600, 500),
  }

  for M in methods:
    if M == 'mle_unconstrained':
      continue
    name = f'{M}_vs_others'
    figs[name] = _plot_single_vs_others(
      results,
      M,
      methods,# - set(('mle_unconstrained',)),
      method_colours,
    )
    export_dims[name] = (600, 485)

  for F in figs.values():
    if 'layout' not in F:
      F['layout'] = {}
    F['layout']['template'] = args.template
  plotter.write_figs(figs, args.plot_fn, export_dims)

if __name__ == '__main__':
  main()
