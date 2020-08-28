import argparse
import numpy as np
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import scipy.special

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import resultserializer

def _calc_eta_stats(truth_fn):
  truth = resultserializer.Results(truth_fn)
  eta = truth.get('eta')

  K, S = eta[1:].shape
  biggest_eta = np.max(eta, axis=1)[1:]
  return (K, S, biggest_eta)

def _load_eta_stats(truth_fns):
  eta_stats = {}
  for T in truth_fns:
    K, S, biggest_eta = _calc_eta_stats(T)
    if (K, S) not in eta_stats:
      eta_stats[(K, S)] = []
    eta_stats[(K, S)].append(biggest_eta)
  for key in list(eta_stats.keys()):
    eta_stats[key] = np.array(eta_stats[key])
  return eta_stats

def plot_eta_stats(eta_stats, alpha, threshold=0.01):
  K_vals = list(sorted(set([K for (K, S) in eta_stats.keys()])))
  S_vals = list(sorted(set([S for (K, S) in eta_stats.keys()])))
  K_map = {K: kidx for kidx, K in enumerate(K_vals)}
  S_map = {S: sidx for sidx, S in enumerate(S_vals)}

  below_thresh_empirical   = np.full((len(K_vals), len(S_vals)), np.nan)
  below_thresh_lol = {}
  below_thresh_theoretical = np.copy(below_thresh_empirical)

  for (K, S), biggest_eta in eta_stats.items():
    num_trees = len(biggest_eta)
    assert biggest_eta.shape == (num_trees, K)
    bt = np.sum(biggest_eta < threshold, axis=1)
    frac = bt/K
    assert len(frac) == num_trees
    below_thresh_empirical[K_map[K],S_map[S]] = np.mean(frac)
    below_thresh_lol[(K,S)] = frac

  for K in K_vals:
    for S in S_vals:
      prob = scipy.special.betainc(alpha, (K-1)*alpha, threshold)**S
      below_thresh_theoretical[K_map[K],S_map[S]] = prob

  fig_emp_hm = make_heatmap(
    S_vals,
    K_vals,
    below_thresh_empirical,
    'Cancer samples',
    'Populations',
    'Mean proportion of populations with η < 0.01<br>across all cancer samples (empirical)',
  )
  fig_thr_hm = make_heatmap(
    S_vals,
    K_vals,
    below_thresh_theoretical,
    'Cancer samples',
    'Populations',
    'Probability that η < 0.01<br>across all cancer samples (theoretical)',
  )

  fig_emp_sp = make_subplots(
    rows = len(K_vals),
    cols = 1,
    shared_xaxes = True,
    vertical_spacing = 0.05,
    row_titles = ['%s subclones' % K for K in K_vals],
    x_title = 'Cancer samples',
    y_title = 'Proportion of populations with η < 0.01<br>across all cancer samples',
  )
  for K in K_vals:
    for S in S_vals:
      if (K, S) not in below_thresh_lol:
        continue
      fig_emp_sp.add_trace(go.Box(y=below_thresh_lol[(K,S)], name=S, boxmean=True), row=K_map[K]+1, col=1)
    fig_emp_sp.add_trace(
      go.Scatter(
        x = S_vals,
        y = below_thresh_theoretical[K_map[K]],
        mode = 'markers',
        marker = {'size': 22, 'opacity': 0.65, 'color': '#000000'},
      ),
      row=K_map[K]+1,
      col=1
    )
  fig_emp_sp.update_xaxes(type = 'category')
  fig_emp_sp.update_yaxes(tickformat = ',.0%', range = (0,1))
  fig_emp_sp.update_layout(showlegend = False)

  fig_thr_sp = make_subplots(
    rows = len(K_vals),
    cols = 1,
    shared_xaxes = True,
    vertical_spacing = 0.09,
    x_title = 'Cancer samples',
    y_title = 'Proportion of populations with η < 0.01 (theoretical)',
  )
  for K in K_vals:
    fig_thr_sp.add_trace(
      go.Bar(x=S_vals, y=below_thresh_theoretical[K_map[K]]),
      row = len(K_vals) - K_map[K],
      col = 1,
    )
  fig_thr_sp.update_xaxes(type = 'category')
  fig_thr_sp.update_yaxes(tickformat = ',.0%', range = (0,1))

  return fig_emp_hm + fig_thr_hm + to_html(fig_emp_sp) + to_html(fig_thr_sp)

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
        'height': 1000,
      },
    },
  )

def make_heatmap(X, Y, Z, xtitle, ytitle, title):
  fig = {
    'data': go.Heatmap(
      x = X,
      y = Y,
      z = Z,
      type = 'heatmap',
      colorscale = 'Viridis',
      zmin = 0,
      zmax = 1,
    ),
    'layout': {
      'xaxis': {'title': xtitle, 'type': 'category'},
      'yaxis': {'title': ytitle, 'type': 'category'},
      'title': title,
      'template': 'plotly_white',
    },
  }

  html = to_html(fig)
  return html

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--alpha', type=float, default=0.1)
  parser.add_argument('plot_fn')
  parser.add_argument('truth_fns', nargs='+')
  args = parser.parse_args()

  eta_stats = _load_eta_stats(args.truth_fns)
  html = plot_eta_stats(eta_stats, args.alpha)
  with open(args.plot_fn, 'w') as F:
    print(html, file=F)

if __name__ == '__main__':
  main()
