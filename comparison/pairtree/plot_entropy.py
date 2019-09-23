import pandas as pd
import argparse
import numpy as np

import plotly
import plotly.graph_objs as go
import re

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

def plot(results, result_key, logy=False):
  results, unique_keys = partition(results, ('K', 'S'))
  unique_keys['K'] = list(reversed(unique_keys['K']))
  rows = len(unique_keys['K'])
  cols = len(unique_keys['S'])

  fig = plotly.subplots.make_subplots(rows=rows, cols=cols, shared_yaxes=True, horizontal_spacing=0, vertical_spacing=0.05)
  for idx, kval in enumerate(unique_keys['K']):
    fig.update_yaxes(title_text=f'K={kval}', row=idx+1, col=1)
  for idx, sval in enumerate(unique_keys['S']):
    fig.update_xaxes(title_text=f'S={sval}', row=rows, col=idx+1)
  fig.update_layout(showlegend=False)
  if logy:
    fig.update_layout(yaxis_type='log')

  for ridx in range(rows):
    for cidx in range(cols):
      K = unique_keys['K'][ridx]
      S = unique_keys['S'][cidx]
      rows = results[(results.K == K) & (results.S == S)]
      fig.add_trace(go.Box(
        y = rows[result_key],
        name = f'(K,S) = ({K},{S}) ({len(rows)} runs)',
      ), row = ridx+1, col = cidx+1)
      fig.update_xaxes(showticklabels=False, row=ridx+1, col=cidx+1)
  fig.update_layout(title_text=result_key)

  html = plotly.offline.plot(
    fig,
    output_type = 'div',
    include_plotlyjs = False,
    config = {
      'showLink': True,
      'toImageButtonOptions': {
        'format': 'svg',
        'width': 750,
        'height': 450,
        'filename': 'entropy',
      },
    },
  )
  return html

def main():
  parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('entropy_fn')
  parser.add_argument('out_fn')
  args = parser.parse_args()

  results = pd.read_csv(args.entropy_fn)
  results['H_trees_pairtree_3_minus_H_trees_truth'] = results['H_trees_pairtree_3'] - results['H_trees_truth']
  results['H_parents_pairtree_minus_H_trees_pairtree_3'] = results['H_parents_pairtree'] - results['H_trees_pairtree_3']
  results['H_parents_truth_minus_H_trees_truth'] = results['H_parents_truth'] - results['H_trees_truth']
  results['H_parents_pairtree_minus_H_parents_truth'] = results['H_parents_pairtree'] - results['H_parents_truth']

  html = ''
  for key, logy in (
    (('true_trees', True)), 
    (('sampled_unique_trees', True)), 
    (('prop_good', False)), 
    (('H_trees_truth', False)), 
    (('H_trees_pairtree_3', False)), 
    (('H_parents_truth', False)), 
    (('H_parents_pairtree', False)), 
    (('H_trees_pairtree_3_minus_H_trees_truth', False)),
    (('H_parents_pairtree_minus_H_parents_truth', False)),
    (('H_parents_pairtree_minus_H_trees_pairtree_3', False)),
    (('H_parents_truth_minus_H_trees_truth', False)),
    (('prop_truth_recovered', False)),
    (('jaccard', False)),
    (('jsd_trees', False)),
    (('jsd_parents_sum', False)),
    (('jsd_parents_mean', False)),
    (('jsd_parents_max', False)),
    (('jsd_parents_phi_mean', False)),
    (('jsd_parents_phi_max', False)),
  ):
    html += plot(results, key, logy)

  with open(args.out_fn, 'w') as outf:
    print( '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>', file=outf)
    print(html, file=outf)

if __name__ == '__main__':
  main()
