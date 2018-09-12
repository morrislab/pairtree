import scipy.special

def logfactorial(X):
  return scipy.special.gammaln(X + 1)

def log_N_choose_K(N, K):
  return logfactorial(N) - (logfactorial(K) + logfactorial(N - K))
