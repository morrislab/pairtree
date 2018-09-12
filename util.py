import scipy.special
import time

def logfactorial(X):
  return scipy.special.gammaln(X + 1)

def log_N_choose_K(N, K):
  return logfactorial(N) - (logfactorial(K) + logfactorial(N - K))

def time_exec(f):
  def wrap(*args):
    time1 = time.time()
    ret = f(*args)
    time2 = time.time()
    ms = (time2-time1)*1000.0
    time_exec._ms = ms
    return ret
  return wrap
time_exec._ms = None
