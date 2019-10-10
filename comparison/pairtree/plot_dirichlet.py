import numpy as np
import scipy.special
import plotly
import plotly.graph_objects as go
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

def compute_analytic_prob(K):
  log_alpha_range = np.linspace(-5, 5, 11)
  alpha = 10**log_alpha_range
  #S = np.array([1] + list(np.linspace(10, 100, steps - 1)))
  S = np.array([1, 3, 10, 30, 100])

  mean = alpha / (K*alpha)
  # Prob for single sample
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

def make_heatmap(X, Y, Z, xtitle, ytitle, title):
  fig = {
    'data': go.Heatmap(
      x = X,
      y = Y,
      z = Z,
      type = 'heatmap',
      colorscale = 'Viridis',
    ),
    'layout': {
      'xaxis': {'title': xtitle},
      'yaxis': {'title': ytitle},
      'title': title,
    },
  }

  html = plotly.offline.plot(
    fig,
    output_type = 'div',
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
  return html

def plot_prob(alpha, K, S, prob_lt):
  print(S)
  html = make_heatmap(
    ['S=%s' % s for s in S],
    ['10<sup>%d</sup>' % b for b in np.log10(alpha)],
    prob_lt,
    'Samples',
    r'$\alpha$',
    r'$\text{Probability that all samples in }K=%s\text{ tree will have } \eta \lt \frac{1}{K} = \frac{1}{%s}$' % (K, K),
  )
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
  alpha, S, prob_lt = compute_analytic_prob(args.K)
  N, eta_sum_var = compute_var(alpha, args.K)

  html = '<h1><marquee>I LOVE TREES</marquee></h1>'
  html += plot_prob(alpha, args.K, S, prob_lt)
  html += plot_var(alpha, N, args.K, eta_sum_var)
  with open(args.plot_fn, 'w') as F:
    print(html, file=F)

if __name__ == '__main__':
  main()
