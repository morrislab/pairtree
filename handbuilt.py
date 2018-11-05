import json
import collections
import numpy as np
import common

def _check_dupes(coll):
  if len(set(coll)) != len(coll):
    raise Exception('Duplicates: %s' % [E for E, cnt in collections.Counter(coll).items() if cnt > 1])

def _load_handbuilt(jsonfn, tree_type):
  with open(jsonfn) as F:
    H = json.load(F)
  return H[tree_type]

def _load_clusters(hb, variants):
  clusters = hb['clusters']

  clustered = [ssmid for C in hb['clusters'] for ssmid in C]
  garbage = hb['garbage']
  _check_dupes(clustered)
  _check_dupes(garbage)

  handbuilt_ssms = set(clustered) | set(garbage)
  handbuilt_ssms = set(['s%s' % S for S in handbuilt_ssms])

  allssms = set(variants.keys())
  if allssms  != handbuilt_ssms:
    raise Exception('all - handbuilt = %s, handbuilt - all = %s' % (str(allssms - handbuilt_ssms), str(handbuilt_ssms - allssms)))
  assert clusters[0] == []

  return clusters

def _load_tree(hb):
  struct = hb['structure']

  children = [node for C in struct.values() for node in C]
  _check_dupes(children)
  # Node 0 should be root, and so should not appear in children list.
  assert set(children) == set(range(1, max(children) + 1))
  assert len(children) + 1 == len(hb['clusters'])

  strkeys = list(struct.keys())
  for key in strkeys:
    struct[int(key)] = struct[key]
    del struct[key]
  return struct

def _calc_dfs_order(tstruct, mean_vafs):
  ordered = []
  def _dfs(node):
    ordered.append(node)
    if node in tstruct:
      for child in sorted(tstruct[node], key=lambda node_idx: -mean_vafs[node_idx]):
        _dfs(child)
  _dfs(0)
  return ordered

def _renumber_clusters(clusters, tstruct, variants, sampnames, colourings):
  vafs = [np.array([variants['s%s' % V]['vaf'] for V in C]) for C in clusters]
  vafs[0] = np.ones(shape=(1, len(sampnames)))
  mean_vafs = np.array([np.mean(V, axis=0) for V in vafs])
  left_sampidx = sampnames.index(colourings[0]['left'])
  left_mean_vafs = mean_vafs[:,left_sampidx]

  reordered = _calc_dfs_order(tstruct, left_mean_vafs)
  new_clusters = [clusters[cidx] for cidx in reordered]
  mapping = {old: new for new, old in enumerate(reordered)}
  new_tstruct = {mapping[old_parent]: [mapping[old_child] for old_child in tstruct[old_parent]] for old_parent in tstruct.keys()}
  return (new_clusters, new_tstruct)

def load_garbage(handbuilt_jsonfn, tree_type):
  hb = _load_handbuilt(handbuilt_jsonfn, tree_type)
  return hb['garbage']

def load_clusters_and_tree(handbuilt_jsonfn, variants, tree_type, sampnames):
  hbjson = _load_handbuilt(handbuilt_jsonfn, tree_type)
  clusters = _load_clusters(hbjson, variants)
  tstruct = _load_tree(hbjson)
  colourings = hbjson['colourings']
  clusters, tstruct = _renumber_clusters(clusters, tstruct, variants, sampnames, colourings)
  adjm = common.convert_adjlist_to_adjmatrix(tstruct)
  return (clusters, adjm, colourings)

def load_clusters(handbuilt_jsonfn, variants, tree_type, sampnames):
  hbjson = _load_handbuilt(handbuilt_jsonfn, tree_type)
  clusters = _load_clusters(hbjson, variants)
  return clusters

def load_tree(handbuilt_jsonfn, tree_type):
  hbjson = _load_handbuilt(handbuilt_jsonfn, tree_type)
  if 'structure' not in hbjson:
    return None
  tstruct = _load_tree(hbjson)
  adjm = common.convert_adjlist_to_adjmatrix(tstruct)
  return adjm

def load_samporders(handbuilt_jsonfn, tree_type):
  hbjson = _load_handbuilt(handbuilt_jsonfn, tree_type)
  samporders = hbjson['samporders']
  for S in samporders:
    assert len(S) == len(set(S)), 'Duplicate sample names: %s' % set([s for s in S if S.count(s) > 1])
  return samporders
