import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import mutrel
from common import Models

def fix_rounding_errors(mrel):
  # Floating point error means that some entires can slightly exceed 1. Ensure
  # this doesn't happen by much.
  #
  # This issue typically arises when performing evaluations, using an "average"
  # mutrel created by summing many constituent weighted mutrels.
  assert np.allclose(1, mrel[mrel > 1])
  mrel = np.minimum(1, mrel)
  return mrel

def add_garbage(posterior, garb_svids):
  if len(garb_svids) == 0:
    return posterior
  assert len(set(posterior.vids) & set(garb_svids)) == 0
  new_vids = posterior.vids + garb_svids
  new_posterior = mutrel.init_mutrel(new_vids)
  G = len(garb_svids)
  M = len(new_posterior.vids)

  # Rather than carefully slicing and dicing the array to set it, just use a
  # series of carefully ordered overwrite operations to put it in the correct
  # state.
  new_posterior.rels[:] = 0
  new_posterior.rels[:,:,Models.garbage] = 1
  diag = range(M)
  new_posterior.rels[diag,diag,:] = 0
  new_posterior.rels[diag,diag,Models.cocluster] = 1
  new_posterior.rels[:-G,:-G,:] = posterior.rels

  mutrel.check_posterior_sanity(new_posterior.rels)
  return new_posterior

def softmax(V):
  V = np.copy(V) - np.max(V)
  return np.exp(V) / np.sum(np.exp(V))

def save_sorted_mutrel(mrel, mrelfn):
  mrel = mutrel.sort_mutrel_by_vids(mrel)
  mutrel.check_posterior_sanity(mrel.rels)
  np.savez_compressed(mrelfn, mutrel=mrel.rels)

