import pandas as pd
import argparse
import numpy as np

import plotly
import re

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import plotter
from plotter import MISSING

S_COLOURS = {S: plotter.hex2rgb(C) for S, C in {
  1: '#ef553b',
  3: '#636efa',
  10: '#00cc96',
  30: '#ab63fa',
  100: '#ffa15a',
}.items()}

def partition(results, keys):
  keys = set(keys)
  key_vals = {K: [] for K in keys}

  for runid in results.runid:
    params = re.findall('[A-Z]\d+', runid)
    for P in params:
      key, val = P[0], int(P[1:])
      if key not in keys:
        continue
      key_vals[key].append(val)

  for K in key_vals.keys():
    results[K] = key_vals[K]
  unique = {K: sorted(set(V)) for K, V in key_vals.items()}
  return (results, unique)

def plot(results, unique_keys, result_key, ytitle, shared_y=False, log_y=False):
  K_vals = unique_keys['K']
  S_vals = unique_keys['S']

  fig = plotly.subplots.make_subplots(
    rows = 1,
    cols = len(K_vals),
    subplot_titles = [plotter.pluralize(K, 'subclone') for K in K_vals],
    shared_yaxes = shared_y,
    x_title = 'Tissue samples',
    y_title = ytitle,
  )
  min_y, max_y = np.inf, -np.inf

  for Kidx, K in enumerate(K_vals):
    for S in S_vals:
      KS_rows = [row for idx, row in results.iterrows() if row['S'] == S and row['K'] == K]
      if len(KS_rows) == 0:
        continue
      trace = {
        'type': 'box',
        'y': [row[result_key] for row in KS_rows],
        'text': [row['runid'] for row in KS_rows],
        'name': '%s' % S,
        'marker': {
          'outliercolor': plotter.format_colour(S_COLOURS[S], 0.5),
          'color': plotter.format_colour(S_COLOURS[S], 0.5),
        },
        'line': {
          'color': plotter.format_colour(S_COLOURS[S]),
        },
      }

      #min_y = np.min([min_y, np.min(trace['y'])])
      #max_y = np.max([max_y, np.max(trace['y'])])
      #if log_y:
      #  trace['y'] = np.log10(trace['y'])
      fig.add_trace(trace, row=1, col=Kidx+1)

  fig.update_layout(
    showlegend = False,
    title_text = ytitle,
  )
  fig.update_xaxes(
    tickangle = 0,
    type = 'category',
  )
  if log_y:
    fig.update_yaxes(type = 'log')
    #floor, ceil = np.floor(np.log10(min_y)), np.ceil(np.log10(max_y)) + 1
    #N = int(ceil - floor)
    #tickvals = np.linspace(floor, ceil, num=(N+1)).astype(np.int)
    #print(tickvals, floor, ceil)
    #assert np.allclose(tickvals, tickvals.astype(np.int))
    #fig.update_yaxes(
    #  tickmode = 'array',
    #  tickvals = tickvals,
    #  ticktext = ['%s' % T for T in tickvals],
    #)

  return fig

def main():
  parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--mutrels', dest='mutrel_fn')
  parser.add_argument('entropy_fn')
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  results = pd.read_csv(args.entropy_fn)
  if args.mutrel_fn is not None:
    mrel_meth, mrel_fn = args.mutrel_fn.split('=', 1)
    mrel = pd.read_csv(mrel_fn)
    mrel = mrel.filter(('runid', mrel_meth))
    mrel = mrel.rename(columns = {mrel_meth: 'mutrel'})
    assert MISSING not in mrel['mutrel']
    assert set(results['runid']) == set(mrel['runid'])
    results = pd.merge(results, mrel, on='runid', how='outer')

  results['H_trees_pairtree_3_minus_H_trees_truth'] = results['H_trees_pairtree_3'] - results['H_trees_truth']
  results, unique_keys = partition(results, ('K', 'S'))

  figs = {}
  export_dims = {}

  for name, title, shared_y, log_y in (
    ('true_trees', 'Trees consistent with<br>true lineage frequencies', False, True), 
    ('jsd_parents_mean', 'Parent JSD between tree parents<br>and Pairtree parents (bits)', True, False),
    ('H_trees_pairtree_3_minus_H_trees_truth', 'Difference in tree entropy distribution<br>between Pairtree and truth (bits)', False, False),
    ('mutrel', 'Pairwise relation error (bits)', True, False),
  ):
    figs[name] = plot(results, unique_keys, name, title, shared_y, log_y)
    export_dims[name] = (700, 400)

  figs['true_trees'].update_yaxes(rangemode = 'tozero')
  #for idx, K in enumerate(unique_keys['K']):
  #  logmax = np.log10(np.max(results['true_trees'][results['K'] == K]))
  #  figs['true_trees'].update_yaxes(range = [-0.1, np.ceil(logmax)], row=1, col=idx+1)
  plotter.write_figs(figs, args.plot_fn, export_dims)

if __name__ == '__main__':
  main()
