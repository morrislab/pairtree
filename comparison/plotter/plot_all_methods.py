import argparse
import numpy as np
import plotter
from plotter import MISSING

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--S-threshold', type=int)
  parser.add_argument('--template', default='seaborn')
  parser.add_argument('--max-y')
  parser.add_argument('--score-type', required=True, choices=('mutrel', 'mutphi', 'mutdistl1', 'mutdistl2', 'cputime', 'walltime'))
  parser.add_argument('--baseline')
  parser.add_argument('--bandwidth', type=float)
  parser.add_argument('--plot-type', choices=('box', 'violin'), default='box')
  parser.add_argument('--log-y-axis', action='store_true')
  parser.add_argument('--hide-method', nargs='+')
  parser.add_argument('results_fn')
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  results, methods = plotter.load_results(args.results_fn)
  if args.hide_method is not None:
    methods = [M for M in methods if M not in args.hide_method]
  plotter.munge(results, methods, args.baseline, args.score_type, args.plot_type)

  sides = {
    'lte': 'negative',
    'gt': 'positive',
  }
  names = {
    'lte': '1 or 3 tissue samples',
    'gt': '10, 30, or 100 tissue samples',
  }
  colours = {
    'lte':'rgb(55,126,184)',
    'gt': 'rgb(152,78,163)',
  }

  score_traces = []
  completion_bar_traces = []
  completion_pie_figs = []

  if args.S_threshold is not None:
    results = plotter.augment(results, 'S')
    parted_by_S = plotter.partition(results, methods, key='S')
    parted_by_thresh = plotter.partition_by_threshold(parted_by_S, args.S_threshold)
    for group in parted_by_thresh.keys():
      vals = parted_by_thresh[group]
      if len(vals) == 0:
        continue
      methods = [M for M in vals.keys()]
      scores = {M: vals[M]['scores'] for M in methods}
      complete = {M: vals[M]['complete'] for M in methods}
      total = {M: vals[M]['total'] for M in methods}
      X, Y = plotter.make_points_for_methods(scores)

      if args.plot_type == 'box':
        score_trace = plotter.make_box_trace(X, Y, group=group, name=names[group], colour=colours[group])
      elif args.plot_type == 'violin':
        score_trace = plotter.make_violin_trace(X, Y, side=sides[group], group=group, name=names[group], bandwidth=args.bandwidth, colour=colours[group])
      else:
        raise Exception('Unknown plot type %s' % args.plot_type)

      score_traces.append(score_trace)
      completion_bar_traces.append(plotter.make_bar_trace(methods, complete, total, name=names[group]))
      completion_pie_figs.append(plotter.make_pie_fig(methods, complete, total, name=names[group]))
  else:
    scores = {M: np.array(results[M]) for M in methods}
    total = {M: len(scores[M]) for M in methods}
    complete = {M: np.sum(scores[M] != MISSING) for M in methods}
    X, Y = plotter.make_points_for_methods(scores)
    if args.plot_type == 'box':
      score_trace = plotter.make_box_trace(X, Y, colour=colours['lte'])
    elif args.plot_type == 'violin':
      score_trace = plotter.make_violin_trace(X, Y, side='both', bandwidth=args.bandwidth)
    score_traces.append(score_trace)
    completion_bar_traces.append(plotter.make_bar_trace(methods, complete, total))
    completion_pie_figs.append(plotter.make_pie_fig(methods, complete, total))

  figs = [
    plotter.make_fig(
      score_traces,
      args.template,
      plotter.make_score_ytitle(args.score_type, args.plot_fn),
      args.max_y,
      log_y_axis = args.log_y_axis,
      layout_options = {
        'boxmode': 'group',
        'violinmode': 'overlay',
        'violingap': 0.0,
        'violingroupgap': 0.0,
      },
    ),

    plotter.make_fig(
      completion_bar_traces,
      args.template,
      'Proportion of failed runs',
      max_y = None,
      layout_options = {'barmode': 'group'},
    ),
  ]
  figs += completion_pie_figs
  plotter.write_figs(figs, args.plot_fn)

if __name__ == '__main__':
  main()
