import numpy as np
import numpy.ma as ma
import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import mutrel
import inputparser
from common import Models

MISSING = -1

def load_mutrels(mutrel_args):
  mutrels = {}
  for mutrel_arg in mutrel_args:
    mutrel_name, mutrel_path = mutrel_arg.split('=', 1)
    assert mutrel_name not in mutrels, '%s is duplicate' % mutrel_name
    if os.path.exists(mutrel_path):
      mrel = np.load(mutrel_path)
      mutrels[mutrel_name] = mutrel.Mutrel(vids=mrel['vids'], rels=mrel['rels'])
    else:
      mutrels[mutrel_name] = None
  return mutrels

def discard_garbage(mutrels, clustered, garbage, ignore_garbage_for):
  for name in list(mutrels.keys()):
    if mutrels[name] is None:
      continue
    vids = mutrels[name].vids
    assert set(vids) == clustered | garbage, 'vids do not match expected set for %s' % name

    gidxs = [idx for idx, vid in enumerate(vids) if vid in garbage]
    if name not in ignore_garbage_for:
      for garbrels in (mutrels[name].rels[gidxs,:,Models.garbage], mutrels[name].rels[:,gidxs,Models.garbage].T):
        # Garbage mutations should have posterior that coclusters with
        # themselves, but that renders them garbage to every other mutation
        # (including other garbage)
        G, M = garbrels.shape
        # `expected` shape, given `G` garbage variants and `M` total variants: `GxM`
        expected = np.ones((G, M))
        expected[np.arange(G),gidxs] = 0
        assert np.allclose(expected, garbrels), '%s garbage relations are wrong' % name
    mutrels[name] = mutrel.remove_variants_by_vidx(mutrels[name], gidxs)
    assert set(mutrels[name].vids) == clustered

def _score_distl1(rels, truth):
  # Compute mean L1 distance.
  M = len(rels)
  dist = np.sum(np.abs(rels - truth), axis=2)
  assert np.allclose(dist, dist.T)
  assert np.allclose(0, np.diag(dist))
  # A pairwise mutrel vector can have a maximum L1 distance of 2, when no
  # elements have overlap with each other. Normalize this so scores can be
  # interpreted as "mean proportion of miscalled relations", given hard 0/1
  # calls for relation types.
  dist /= 2
  # Distances may be slightly higher than 1 because of floating point error.
  # Set these to exactly 1.
  dist[np.logical_and(dist > 1, np.isclose(1, dist))] = 1
  assert np.all(0 <= dist) and np.all(dist <= 1)

  # Take entries below main diagonal.
  dist_lower = np.tril(dist, 1)
  assert dist_lower.shape == (M, M)
  # There are (M choose 2) elements below the diagonal, so divide by this
  # when computing mean.
  score = np.sum(dist_lower) / (0.5*M*(M - 1))
  return score

def _compute_kld(P, Q):
  for A in (P, Q):
    assert np.all(A >= 0)
    assert np.allclose(1, np.sum(A, axis=2))

  logP = ma.log2(ma.masked_equal(P, 0))
  logQ = ma.log2(ma.masked_equal(Q, 0))
  kld = np.sum(P * (logP - logQ), axis=2)

  assert np.allclose(0, kld[kld < 0])
  kld = np.abs(kld)
  assert np.all(kld >= 0)
  return kld

def _score_jsd(rels, truth):
  M = len(rels)
  R = 0.5*(rels + truth)

  kld1 = _compute_kld(rels, R)
  kld2 = _compute_kld(truth, R)
  jsd = 0.5*(kld1 + kld2)
  assert np.allclose(1, jsd[jsd > 1])
  jsd = np.minimum(1, jsd)

  assert np.allclose(jsd, jsd.T)
  joint_jsd = np.sum(np.tril(jsd, 1))
  mean_jsd = joint_jsd / (0.5*M*(M-1))
  return mean_jsd

def compare(mutrels):
  assert 'truth' in mutrels
  M, _, num_models = mutrels['truth'].rels.shape
  assert mutrels['truth'].rels.shape == (M, M, num_models)
  correct = np.argmax(mutrels['truth'].rels, axis=2)

  names = sorted(mutrels.keys())
  scores = {}

  for name in names:
    mrel = mutrels[name]
    if mrel is None:
      scores[name] = MISSING
      continue

    assert mrel.rels.shape == (M, M, num_models)
    assert np.array_equal(mrel.vids, mutrels['truth'].vids)
    mutrel.check_posterior_sanity(mrel.rels)

    scores[name] = _score_jsd(mrel.rels, mutrels['truth'].rels)
    if name == 'truth':
      assert np.isclose(0, scores[name])
    
  print(*names, sep=',')
  print(*[scores[name] for name in names], sep=',')

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--discard-garbage', dest='discard_garbage', action='store_true')
  # Normally, we check to ensure that garbage mutations have all their
  # posterior mass on the garbage relation. However, if a method has been
  # allowed to cluster mutations itself, this won't necessarily hold true, so
  # permit this behaviour to be disabled for those methods.
  parser.add_argument('--ignore-garbage-for', action='append')
  parser.add_argument('--params', dest='paramsfn', required=True)
  parser.add_argument('mutrels', nargs='+')
  args = parser.parse_args()

  params = inputparser.load_params(args.paramsfn)
  garbage = set(params['garbage'])
  clustered = set([vid for C in params['clusters'] for vid in C])

  mutrels = load_mutrels(args.mutrels)
  if args.discard_garbage:
    ignore_garbage_for = set(args.ignore_garbage_for) if args.ignore_garbage_for is not None else set()
    discard_garbage(mutrels, clustered, garbage, ignore_garbage_for)
  compare(mutrels)

if __name__ == '__main__':
  main()
