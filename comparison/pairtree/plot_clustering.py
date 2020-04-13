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

  fig = plotly.subplots.make_subplots(rows=rows, cols=cols, shared_yaxes=False, horizontal_spacing=0.05, vertical_spacing=0.05)
  for idx, kval in enumerate(unique_keys['K']):
    fig.update_yaxes(title_text=f'K={kval}', row=idx+1, col=1)
  for idx, sval in enumerate(unique_keys['S']):
    fig.update_xaxes(title_text=f'S={sval}', row=rows, col=idx+1)
  fig.update_layout(showlegend=False)

  for ridx in range(rows):
    for cidx in range(cols):
      K = unique_keys['K'][ridx]
      S = unique_keys['S'][cidx]
      rows = results[(results.K == K) & (results.S == S)]
      fig.add_trace(go.Box(
        y = rows[result_key],
        name = f'(K,S) = ({K},{S}) ({len(rows)} runs)',
        boxmean = True,
      ), row = ridx+1, col = cidx+1)
      fig.update_xaxes(showticklabels=False, row=ridx+1, col=cidx+1)
  fig.update_layout(title_text=result_key)
  if logy:
    fig.update_yaxes(type='log')

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
        'filename': result_key,
      },
    },
  )
  return html

def main():
  parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('clusterstats_fn')
  parser.add_argument('out_fn')
  args = parser.parse_args()

  results = pd.read_csv(args.clusterstats_fn)
  results['found_nlglh_minus_true_nlglh'] = results['found_nlglh'] - results['true_nlglh']
  results['found_clusters_minus_true_clusters'] = results['found_clusters'] - results['true_clusters']

  html = ''
  for key, logy in (
    ('found_nlglh_minus_true_nlglh', False),
    ('found_clusters_minus_true_clusters', False),
    ('homogeneity', False),
    ('completeness', False),
    ('vmeasure', False),
    ('ami', False),
  ):
    html += plot(results, key, logy)

  with open(args.out_fn, 'w') as outf:
    print( '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>', file=outf)
    print(html, file=outf)

if __name__ == '__main__':
  main()
