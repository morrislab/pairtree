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
  # The vector sum will be close to 1 at this point, but not close enough to
  # make np.random.choice happy -- sometimes it will issue the error
  # "probabilities do not sum to 1" from mtrand.RandomState.choice.
  # To fix this, renormalize the vector.
  smax /= np.sum(smax)
  assert np.isclose(np.sum(smax), 1)
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

def find_parents(adj):
  adj = np.copy(adj)
  np.fill_diagonal(adj, 0)
  assert np.all(np.sum(adj[:,1:], axis=0) == 1)
  return np.argmax(adj[:,1:], axis=0)

def convert_parents_to_adjmatrix(parents):
  K = len(parents) + 1
  adjm = np.eye(K)
  adjm[parents,np.arange(1, K)] = 1
  return adjm

def _calc_phi_hat(variants):
  _extract_mat = lambda key: np.array([var[key] for var in variants])

  V = _extract_mat('var_reads')
  T = _extract_mat('total_reads')
  omega = _extract_mat('omega_v')
  phi_hat = (V / T) / omega
  phi_hat = np.minimum(1, phi_hat)
  phi_hat = np.insert(phi_hat, 0, 1, axis=0)
  assert np.all(phi_hat >= 0)
  return phi_hat

def calc_nlglh(llh, K, S):
  return -llh / (np.log(2) * (K-1) * S)

def make_tree_struct(struct, count, llh, prob, phi, variants, sampnames):
  K, S = phi.shape
  phi_hat = _calc_phi_hat(variants)
  tree = {
    'phi': phi.tolist(),
    'phi_hat': phi_hat.tolist(),
    'llh': float(llh),
    'nlglh': float(calc_nlglh(llh, K, S)),
    'prob': float(prob),
    'count': int(count),
    'parents': struct.tolist(),
    'samples': sampnames,
  }
  return tree
