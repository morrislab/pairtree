import numpy.ma as ma
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import util
import evalutil

def compute_parent_dist(adjms, llhs):
  K = len(adjms[0])
  parent_dist = np.zeros((K, K))
  weights = util.softmax(llhs)

  for adjm, weight in zip(adjms, weights):
    np.fill_diagonal(adjm, 0)
    parent_dist += weight * adjm

  assert np.all(0 == parent_dist[:,0])
  parent_dist = parent_dist[:,1:]
  parent_dist = evalutil._fix_rounding_errors(parent_dist)

  assert np.all(0 <= parent_dist) and np.all(parent_dist <= 1)
  assert np.allclose(1, np.sum(parent_dist, axis=0))
  return parent_dist

def compute_entropy(parent_dist):
  K = len(parent_dist)
  assert parent_dist.shape == (K, K-1)
  parent_dist = ma.masked_equal(parent_dist, 0)
  parent_entropy = -ma.sum(parent_dist * ma.log2(parent_dist), axis=0)
  total_entropy = np.sum(parent_entropy)
  entropy_per_node = total_entropy / (K-1)
  return (total_entropy, entropy_per_node, parent_entropy)

# TODO: remove the `load` and `save` methods?
def load_pardist(pardist_path):
  return np.load(pardist_path)['parent_dist']

def save_pardist(pardist, pardist_path):
  np.savez_compressed(pardist_path, parent_dist = pardist)
