from collections import namedtuple
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import util
import evalutil
import common

Mutstat = namedtuple('Mutstat', ('vids', 'assays', 'stats'))
MISSING = -1

def write(mstat, mutstatfn):
  # calc_mutstat should have created `mutstat` with sorted vids, but double-check
  # this is true.
  assert list(mstat.vids) == common.sort_vids(mstat.vids)
  np.savez_compressed(mutstatfn, stats=mstat.stats, vids=mstat.vids, assays=mstat.assays)

def load(mutstatfn):
  results = np.load(mutstatfn, allow_pickle=True)
  return Mutstat(vids=results['vids'], assays=results['assays'], stats=results['stats'])

def load_mutstats(mstat_args, check_inf=True, check_nan=True):
  mutstats = {}
  for mstat_arg in mstat_args:
    name, mstat_path = mstat_arg.split('=', 1)
    assert name not in mutstats
    if os.path.exists(mstat_path):
      mutstats[name] = load(mstat_path)
      if check_inf:
        assert not np.any(np.isinf(mutstats[name].stats)) 
      if check_nan:
        assert not np.any(np.isnan(mutstats[name].stats))
    else:
      mutstats[name] = None
  return mutstats

def remove_garbage(mutstats, garbage):
  # Remove garbage if present. Some mutstats have it (e.g., PWGS run on all
  # variants, without benefit of handbuilt clusters) while others don't.
  garbage = set(garbage)
  revised = {}

  for name, mstat in mutstats.items():
    if mstat is None:
      revised[name] = None
      continue
    garb_idxs = [idx for idx, vid in enumerate(mstat.vids) if vid in garbage]
    new_vids = [vid for idx, vid in enumerate(mstat.vids) if idx not in set(garb_idxs)]
    new_stats = np.delete(mstat.stats, garb_idxs, axis=0)
    revised[name] = Mutstat(stats=new_stats, vids=new_vids, assays=mstat.assays)

  return revised

def check_incomplete(mutstats, clusters):
  clustered = set([vid for C in clusters for vid in C])
  names = list(mutstats.keys())
  for name in names:
    if mutstats[name] is None:
      continue
    vids = set(mutstats[name].vids)
    assert vids.issubset(clustered)
    if vids != clustered:
      missing = clustered - vids
      msg = '%s lacks fraction=%s variants (%s)' % (name, len(missing) / len(clustered), missing)
      mutstats[name] = None
      raise Exception(msg)

def score_mutstats(mutstats, _score):
  names = list(mutstats.keys())
  scores = {N: MISSING for N in names}
  present = [N for N in names if mutstats[N] is not None]
  if len(present) == 0:
    return (names, scores)
  first_present = present[0]
  vids = mutstats[first_present].vids
  assays = mutstats[first_present].assays

  for name in present:
    mstat = mutstats[name]
    assert np.array_equal(mstat.vids, vids)
    assert np.array_equal(mstat.assays, assays)
    scores[name] = _score(mstat.stats)
  return (names, scores)

def impute_garbage(mstats, garbage, _impute):
  if len(garbage) == 0:
    return mstats
  garbage = set(garbage)
  old_vids = set(mstats.vids)
  assert len(old_vids & garbage) == 0

  new_vids = common.sort_vids(old_vids | garbage)
  M_old, S = mstats.stats.shape
  M_new = len(new_vids)
  new_stats = np.full((M_new, S), np.nan)
  new_assays = mstats.assays
  idxmap = {vid: idx for idx, vid in enumerate(mstats.vids)}

  for idx, vid in enumerate(new_vids):
    if vid in old_vids:
      new_stats[idx] = mstats.stats[idxmap[vid]]
    else:
      new_stats[idx] = _impute(vid)

  assert not np.any(np.isnan(new_stats))
  return Mutstat(vids=new_vids, stats=new_stats, assays=new_assays)
