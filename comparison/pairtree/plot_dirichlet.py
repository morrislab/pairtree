import numpy as np
import scipy.special
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.figure_factory as ff
import argparse

def compute_empirical_prob(alpha, K, S):
  eta = np.random.dirichlet(alpha = K*[alpha], size = S).T
  mean = alpha / (K*alpha)
  bigger = np.max(eta, axis=1) >= mean
  return np.sum(bigger) / K

def compute_empirical_repeatedly(alpha, K, S, trials):
  probs = []
  for _ in range(trials):
    result = compute_empirical_prob(alpha, K, S)
    probs.append(result)
  probs = np.array(probs)
  return probs

def compute_analytic_prob_single(K):
  log_alpha_range = np.linspace(-5, 5, 11)
  alpha = 10**log_alpha_range
  #S = np.array([1] + list(np.linspace(10, 100, steps - 1)))
  S = np.array([1, 3, 10, 30, 100])

  mean = alpha / (K*alpha)
  # Prob for single sample. Relies on the beta being the marginal for a
  # Dirichlet component.
  prob_lt_ss = scipy.special.betainc(
    alpha,
    (K-1)*alpha,
    mean,
  )
  # prob_lt[i,j] = probability that component `i` of Dirichlet (parameterized
  # by `alpha_i`) will have value <= `mean` in every Dirichlet sample, if we
  # take `S_j` (independent) samples from that Dirichlet
  S_rows = np.dot(np.ones(len(log_alpha_range))[:,None], S[None,:])
  prob_cols = np.dot(prob_lt_ss[:,None], np.ones(len(S))[None,:])
  assert S_rows.shape == prob_cols.shape
  prob_lt = prob_cols**S_rows
  #prob_lt = np.hstack((alpha[:,None], (K-1)*alpha[:,None], prob_lt))
  return (alpha, S, prob_lt)

def compute_analytic_prob_many(K, C, alpha, S, threshold=0.01):
  # prob_lt[i,j] = probability that component `i` of Dirichlet (parameterized
  # by `alpha_i`) will have value <= `mean` in every Dirichlet sample, if we
  # take `S_j` (independent) samples from that Dirichlet
  prob_lt = np.full((len(C), len(S)), np.nan)
  for cidx, c in enumerate(C):
    for sidx, s in enumerate(S):
      prob_lt[cidx,sidx] = scipy.special.betainc(c*alpha, (K - c)*alpha, threshold)**s

  return prob_lt

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
    },
  }

  return to_html(fig)

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

def plot_prob_single(alpha, K, S, prob_lt):
  html = make_heatmap(
    ['S=%s' % s for s in S],
    ['10<sup>%d</sup>' % b for b in np.log10(alpha)],
    prob_lt,
    'Samples',
    r'$\alpha$',
    r'$\text{Probability that all samples in }K=%s\text{ tree will have an } \eta \lt \frac{1}{K} = \frac{1}{%s}$' % (K, K),
  )
  return html

def _make_c_range(K):
  rng = {
    # Formula: to span 1 to K-1, there are K-2 gaps between integers. For `num`
    # values that yield integers, specify `num = [divisor of K-2] + 1`.
    3:   np.linspace(1, K-1, num=2).astype(np.int),
    10:  np.linspace(1, K-1, num=5).astype(np.int),
    30:  np.linspace(1, K-1, num=8).astype(np.int),
    100: np.linspace(1, K-1, num=8).astype(np.int),
  }
  return rng[K]

def plot_prob_many(alpha, K, threshold=0.05):
  C = np.arange(1, K)
  S = [1, 3, 10, 30, 100]
  prob = compute_analytic_prob_many(K, C, alpha, S)

  visible = np.any(prob >= threshold, axis=1)

  Z = prob[visible]
  labels = [ [f'{100*V:.0f}%' for V in row] for row in Z ]
  fig = ff.create_annotated_heatmap(
    Z,
    annotation_text=labels,
    x = ['%s' % x for x in S],
    y = ['%s' % y for y in C[visible]],
    colorscale = 'Viridis',
    zmin = 0,
    zmax = 1,
    showscale = True,
    colorbar = {'tickformat': ',.0%'},
  )
  fig.update_xaxes(title = 'Cancer samples', type='category')
  fig.update_yaxes(title = 'Populations in sub-tree', type='category')
  fig.update_layout(title = 'Probability that all populations in sub-tree will have η < 0.01<br>across all cancer samples')
  html = to_html(fig)

  rows = int(np.sum(visible))
  bars = make_subplots(
    rows = rows,
    cols = 1,
    shared_xaxes = True,
    vertical_spacing = 0.05,
    row_titles = [f'{c} pop{"" if c == 1 else "s"}.' % c for c in C[visible]],
    x_title = 'Cancer samples',
    y_title = 'Probability that all populations in sub-tree<br>will have η < 0.01 across all cancer samples',
  )
  for idx in range(rows):
    bars.add_trace(
      go.Bar(x=S, y=prob[visible][idx]),
      row = rows - idx,
      col = 1
    )
  bars.update_layout(showlegend = False)
  bars.update_xaxes(type = 'category')
  bars.update_yaxes(tickformat = ',.0%', range = (0,1))
  html += to_html(bars)

  return html

def plot_var(alpha, N, K, eta_sum_var):
  html = make_heatmap(
    ['N=%s' % n for n in N],
    ['10<sup>%d</sup>' % b for b in np.log10(alpha)],
    eta_sum_var,
    '$N$',
    r'$\alpha$',
    'Stdev of sum of N components in %s-element Dirichlet' % K,
  )
  return html

def compute_var(alpha, K):
  N = np.arange(1, K+1)

  N_tiled = np.dot(np.ones(len(alpha))[:,None], N[None,:])
  alpha_sum = K*alpha
  alpha_sum_tiled = np.dot(alpha_sum[:,None], np.ones(K)[None,:])
  assert N_tiled.shape == alpha_sum_tiled.shape
  # prob_lt[i,j] = probability that component `i` of Dirichlet (parameterized
  # by `alpha_i`) will have value <= `mean` in every Dirichlet sample, if we
  # take `S_j` (independent) samples from that Dirichlet
  eta_sum_var = np.sqrt( (N_tiled/K) * (1 - N_tiled/K) / (alpha_sum_tiled + 1) )
  return (N, eta_sum_var)

def main():
  import sys
  np.set_printoptions(linewidth=400, precision=3, threshold=sys.maxsize, suppress=True)
  np.seterr(divide='raise', invalid='raise')

  parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('-K', default=100, type=int)
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  #compute_empirical_repeatedly(alpha=0.001, K=args.K, S=10, trials=100)
  alpha, S, prob_lt = compute_analytic_prob_single(args.K)
  N, eta_sum_var = compute_var(alpha, args.K)

  html = '<h1><marquee>I FUCKING HATE TREES</marquee></h1>'
  html += plot_prob_single(alpha, args.K, S, prob_lt)
  html += plot_var(alpha, N, args.K, eta_sum_var)
  html += plot_prob_many(alpha=0.1, K=args.K)
  with open(args.plot_fn, 'w') as F:
    print(html, file=F)

if __name__ == '__main__':
  main()
