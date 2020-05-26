import numpy as np
import time
from numba import njit, vectorize
from common import Models
import math

@vectorize
def lgamma(X):
  return math.lgamma(X)

@njit
def logfactorial(X):
  result = lgamma(X + 1)
  return result

@njit
def log_N_choose_K(N, K):
  return logfactorial(N) - (logfactorial(K) + logfactorial(N - K))

@njit
def lbeta(X, Y):
  return lgamma(X) + lgamma(Y) - lgamma(X + Y)

@njit
def beta_binom_logpmf(x, n, a, b):
  # Numba doesn't like this assertion. Bleh.
  #assert np.all(a > 0.) and np.all(b > 0.)
  logpx = log_N_choose_K(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b)
  return logpx

@njit
def softmax(V):
  B = np.max(V)
  log_sum = B + np.log(np.sum(np.exp(V - B)))
  log_softmax = V - log_sum
  smax = np.exp(log_softmax)
  # The vector sum will be close to 1 at this point, but not close enough to
  # make np.random.choice happy -- sometimes it will issue the error
  # "probabilities do not sum to 1" from mtrand.RandomState.choice.
  # To fix this, renormalize the vector.
  smax /= np.sum(smax)
  #assert np.isclose(np.sum(smax), 1)
  return smax

@njit
def sample_multinom(probs):
  S = 0.
  R = np.random.rand()
  for idx, P in enumerate(probs):
    S += P
    if S >= R:
      return idx
  # Should not reach this point.
  assert False

@njit
def isclose(arr, V, tol=1e-6):
  # This is a poor person's version of `np.isclose`. Once Numba implements
  # this, I should remove it.
  return np.abs(arr - V) <= tol

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

def convert_adjmatrix_to_parents(adj):
  adj = np.copy(adj)
  np.fill_diagonal(adj, 0)
  return np.argmax(adj[:,1:], axis=0)

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
  eta = calc_eta(struct, phi)
  tree = {
    'phi': phi.tolist(),
    'phi_hat': phi_hat.tolist(),
    'eta': eta.tolist(),
    'llh': float(llh),
    'nlglh': float(calc_nlglh(llh, K, S)),
    'prob': float(prob),
    'count': int(count),
    'parents': struct.tolist(),
    'samples': sampnames,
  }
  return tree

def lpdist(A, B, p=1):
  dist = np.sum(np.abs(A - B)**p)**(1/p)
  assert dist >= 0
  return dist

@njit
def make_ancestral_from_adj(adj, check_validity=False):
  K = len(adj)
  root = 0

  if check_validity:
    # By default, disable checks to improve performance.
    assert np.all(1 == np.diag(adj))
    expected_sum = 2 * np.ones(K)
    expected_sum[root] = 1
    assert np.array_equal(expected_sum, np.sum(adj, axis=0))

  Z = np.copy(adj)
  np.fill_diagonal(Z, 0)
  stack = [root]
  while len(stack) > 0:
    P = stack.pop()
    C = np.flatnonzero(Z[P])
    if len(C) == 0:
      continue
    # Set ancestors of `C` to those of their parent `P`.
    C_anc = np.copy(Z[:,P])
    C_anc[P] = 1
    # Turn `C_anc` into column vector.
    Z[:,C] = np.expand_dims(C_anc, 1)
    stack += list(C)
  np.fill_diagonal(Z, 1)

  if check_validity:
    assert np.array_equal(Z[root], np.ones(K))
  return Z

@njit
def compute_node_relations(adj, check_validity=False):
  K = len(adj)
  anc = make_ancestral_from_adj(adj, check_validity)
  np.fill_diagonal(anc, 0)

  R = np.full((K, K), Models.diff_branches, dtype=np.int8)
  for idx in range(K):
    R[idx][anc[idx]   == 1] = Models.A_B
    R[idx][anc[:,idx] == 1] = Models.B_A
  np.fill_diagonal(R, Models.cocluster)

  if check_validity:
    assert np.all(R[0]   == Models.A_B)
    assert np.all(R[:,0] == Models.B_A)
  return R

def calc_eta(struct, phi):
  K, M = phi.shape
  assert len(struct) == K - 1
  adj = convert_parents_to_adjmatrix(struct)
  Z = make_ancestral_from_adj(adj)
  Z_inv = np.linalg.inv(Z)
  eta = np.dot(Z_inv, phi)

  assert np.allclose(0, eta[eta < 0])
  eta = np.abs(eta)
  # Renormalize.
  eta = eta / np.sum(eta, axis=0)
  return eta
