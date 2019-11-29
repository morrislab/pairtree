import pickle
from collections import namedtuple
import numpy as np

Neutree = namedtuple('Neutree', ('structs', 'phis', 'counts', 'logscores', 'clusterings', 'garbage'))

def save(ntree, outfn):
  # `phis`, `structs`, and `clusterings` can be NumPy arrays (if all of same
  # dimension) or lists of lists/arrays (if different dimensions).
  N = len(ntree.structs)
  for K in ('structs', 'phis', 'counts', 'logscores', 'clusterings'):
    assert len(getattr(ntree, K)) == N, '%s has length %s instead of %s' % (K, len(getattr(ntree, K)), N)

  # Ensure these values are definitely NumPy arrays.
  arr_vals = {K: np.array(getattr(ntree, K)) for K in ('counts', 'logscores')}
  ntree = ntree._replace(**arr_vals)

  with open(outfn, 'wb') as F:
    pickle.dump(ntree, F)

def load(fn):
  with open(fn, 'rb') as F:
    return pickle.load(F)
