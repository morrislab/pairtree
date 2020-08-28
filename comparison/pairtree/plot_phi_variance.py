import argparse
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import util
import resultserializer

def _calc_depth(parents):
  K = len(parents) + 1
  depth = np.nan + np.zeros(K, dtype=np.int)
  depth[0] = 0
  stack = [0]

  while len(stack) > 0:
    node = stack.pop()
    children = np.flatnonzero(parents == node) + 1
    depth[children] = depth[node] + 1
    stack += children.tolist()

  assert not np.any(np.isnan(depth))
  return depth[1:].astype(np.int)

def _calc_num_pops(parents):
  # Calculate number of populations in each subclone.
  adj = util.convert_parents_to_adjmatrix(parents)
  K = len(adj)
  assert adj.shape == (K,K)
  anc = util.make_ancestral_from_adj(adj)
  C = np.sum(anc, axis=1)
  assert C[0] == K
  return C[1:].astype(np.int)

def _calc_tree_stats(truth_fn):
  truth = resultserializer.Results(truth_fn)
  eta = truth.get('eta')
  phi = truth.get('phi')
  struct = truth.get('structure')

  phi_std = np.std(phi, axis=1)
  phi_mean = np.mean(phi, axis=1)
  depth = _calc_depth(struct)
  num_pops = _calc_num_pops(struct)

  df = pd.DataFrame({
    'phi_std': phi_std[1:],
    'phi_mean': phi_mean[1:],
    'largest_eta': np.max(eta, axis=1)[1:],
    'depth': depth,
    'num_pops': num_pops,
  })
  polyprimary = np.sum(struct == 0) > 1
  return (df, polyprimary)

def _process(truth_fns):
  tree_stats = pd.DataFrame()
  polyprimary, monoprimary = 0, 0
  for T in truth_fns:
    ts, is_pp = _calc_tree_stats(T)
    tree_stats = tree_stats.append(ts)
    if is_pp:
      polyprimary += 1
    else:
      monoprimary += 1
  return (tree_stats, polyprimary, monoprimary)

def to_html(fig):
  return pio.to_html(
    fig,
    full_html = False,
    include_plotlyjs = 'cdn',
    include_mathjax = 'cdn',
    config = {
      'showLink': True,
      'toImageButtonOptions': {
        'format': 'svg',
        'width': 750,
        'height': 450,
      },
    },
  )

def count_num_pops(df):
  counts = df['num_pops'].to_numpy()
  assert 0 not in counts
  vals = list(range(1, np.max(counts) + 1))
  cnt = [np.sum(counts == val) for val in vals]
  return (vals, cnt)

def make_box(df, x, y, title, xtitle, ytitle):
  # Similar to px.box(), but with categorical axis that's properly ordered.
  # px.box() seems to place categorical labels in random order.
  rng = list(range(df[x].min(), df[x].max() + 1))
  boxes = [go.Box(y=df[df[x] == N][y], name=N, boxmean=True, marker_color='#6773fa') for N in rng]
  fig = go.Figure(
    data=boxes,
    layout={
      'xaxis': {'title': xtitle, 'type': 'category'},
      'yaxis': {'title': ytitle},
      'title': title,
      'template': 'plotly_white',
    },
  )
  return fig

def plot_tree_stats(df, alpha, K):
  figs = []

  figs.append(make_box(
    df,
    x='num_pops',
    y='phi_std',
    title = f'φ standard deviation for {K} populations ({len(df):,} samples)',
    xtitle='Number of populations in subclone',
   ytitle='φ standard deviation',
  ))
  figs.append(make_box(
    df,
    x='num_pops',
    y='phi_mean',
    title = f'Mean φ for {K} populations ({len(df):,} samples)',
    xtitle='Number of populations in subclone',
    ytitle='φ mean',
  ))
  for F in figs:
    F.update_layout(showlegend = False)
    F.update_xaxes(type = 'category')

  K += 1
  C = np.arange(1, K)
  sum_var = np.sqrt( (C/K) * (1 - C/K) / (K*alpha + 1) )
  sum_mean = C/K
  analytical = pd.DataFrame({'C': C, 'sum_var': sum_var, 'sum_mean': sum_mean})
  figs[0].add_scatter(x=analytical.C, y=analytical.sum_var,  mode='markers', marker={'size': 22, 'opacity': 0.8, 'color': '#ff483b'})
  figs[1].add_scatter(x=analytical.C, y=analytical.sum_mean, mode='markers', marker={'size': 22, 'opacity': 0.8, 'color': '#ff483b'})

  subclone_sizes = count_num_pops(df)
  figs.append({
    'data': go.Bar(x=subclone_sizes[0], y=subclone_sizes[1]),
    'layout': {
      'xaxis': {'title': 'Number of populations in subclone', 'type': 'category'},
      'yaxis': {'title': 'Count'},
    }
  })
  figs.append(px.histogram(
    df,
    x = 'largest_eta',
    labels = {
      'largest_eta': 'Largest η for subpopulation across cancer samples',
    },
  ))
  figs[-1].update_xaxes(range = (0, 1))
  figs[-1].update_yaxes(title = 'Count')

  html = ''
  for idx, F in enumerate(figs):
    html += to_html(F)
  return html

def plot_polyprimary(polyprimary, monoprimary):
  fig = {
    'data': go.Pie(labels=['Monoprimary', 'Polyprimary'], values=[monoprimary, polyprimary], sort=False),
    'layout': { },
  }
  return to_html(fig)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--alpha', type=float, default=0.1)
  parser.add_argument('-K', type=int, required=True)
  parser.add_argument('plot_fn')
  parser.add_argument('truth_fns', nargs='+')
  args = parser.parse_args()

  tree_stats, polyprimary, monoprimary = _process(args.truth_fns)
  ts = plot_tree_stats(tree_stats, args.alpha, args.K)
  pp = plot_polyprimary(polyprimary, monoprimary)
  with open(args.plot_fn, 'w') as F:
    print(ts, file=F)
    print(pp, file=F)

if __name__ == '__main__':
  main()
