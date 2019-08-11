import numpy as np
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

def softmax(V):
  log_softmax = V - scipy.special.logsumexp(V)
  smax = np.exp(log_softmax)
  assert np.isclose(np.sum(smax), 1)
  # The vector sum will be close to 1 at this point, but not close enough to
  # make np.random.choice happy -- sometimes it will issue the error
  # "probabilities do not sum to 1" from mtrand.RandomState.choice.
  # To fix this, renormalize the vector.
  smax /= np.sum(smax)
  return smax

def remove_rowcol(arr, indices):
  '''Remove rows and columns at `indices`.'''
  # Using a mask requires only creating a new array once (and then presumably
  # another array to apply the mask). Calling np.delete() multiple times would
  # create separate arrays for every element in `indices`.
  shape = list(arr.shape)
  # Must convert to list *before* array, as indices is often a set, and
  # converting a set directly to an array yields a dtype of `object` that
  # contains the set. It's really weird.
  indices = np.array(list(indices))
  # Only check the first two dimensions, as we ignore all others.
  assert np.all(0 <= indices) and np.all(indices < min(shape[:2]))

  for axis in (0, 1):
    arr = np.delete(arr, indices, axis=axis)
    shape[axis] -= len(indices)

  assert np.array_equal(arr.shape, shape)
  return arr
