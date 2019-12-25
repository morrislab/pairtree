from scipy.special import gammaln, betaln
import scipy.stats
import numpy as np
import plotly.express as px
import math

def logfactorial(X):
  result = np.empty(X.shape)
  for idx in range(len(X)):
    result[idx] = math.lgamma(X[idx] + 1)
  return np.array(result)

def log_N_choose_K(N, K):
  return logfactorial(N) - (logfactorial(K) + logfactorial(N - K))

#def logpmf(k, n, a, b):
#  logkx = (gammaln(n + 1) + gammaln(k + a) + gammaln(n - k + b) + gammaln(a + b)) - (gammaln(k + 1) + gammaln(n - k + 1) + gammaln(n + a + b) + gammaln(a) + gammaln(b))
#  return logkx

def logpmf(x, n, a, b):
  assert np.all(a > 0) and np.all(b > 0)
  logpx = log_N_choose_K(n, x) + betaln(x + a, n - x + b) - betaln(a, b)
  return logpx

def logpx_betabinom(var_reads, total_reads, omega, phi):
  delta = 0.5
  alpha = var_reads + delta
  beta = np.maximum(delta, omega*total_reads - var_reads + delta)
  return logpmf(var_reads, total_reads, alpha, beta)

def logpx_binom(var_reads, total_reads, omega, phi, epsilon=1e-5):
  P = omega*phi
  P = np.maximum(P, epsilon)
  P = np.minimum(P, 1 - epsilon)
  return scipy.stats.binom.logpmf(var_reads, total_reads, P) 

def main():
  import sys
  np.set_printoptions(linewidth=400, precision=3, threshold=sys.maxsize, suppress=True)
  np.seterr(divide='raise', invalid='raise')

  T = np.array([9000])
  x = np.arange(0, T[0] + 1)
  omega = np.array([0.5])
  phi = np.array([0])

  print('total_reads', T, sep='\t')
  print('omega', omega, sep='\t')
  print('phi', phi, sep='\t')
  limit = 20

  result_betabinom = logpx_betabinom(x, T, omega, phi)
  figs = []
  figs.append(px.line(x=x[:limit], y=np.exp(result_betabinom[:limit]), title='Beta binomial'))
  print(np.exp(scipy.special.logsumexp(result_betabinom)), np.exp(result_betabinom[:limit]))

  for epsilon in (1e-5, 1e-30):
    result = logpx_binom(x, T, omega, phi, epsilon=epsilon)
    fig = px.line(x=x[:limit], y=np.exp(result[:limit]), title=f'Binomial with epsilon={epsilon}')
    figs.append(fig)
    print(np.exp(scipy.special.logsumexp(result)), np.exp(result[:limit]))

  for F in figs:
    F.update_layout(yaxis = {'range': [0,1]})
    #F.show()


main()
