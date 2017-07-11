import json
import collections
import numpy as np

def _check_dupes(coll):
  if len(set(coll)) != len(coll):
    raise Exception('Duplicates: %s' % [E for E, cnt in collections.Counter(coll).items() if cnt > 1])

def _load_handbuilt(jsonfn):
  with open(jsonfn) as F:
    H = json.load(F)
  return H

def load_clusters(handbuilt_jsonfn, variants):
  hb = _load_handbuilt(handbuilt_jsonfn)
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

def load_garbage(handbuilt_jsonfn):
  hb = _load_handbuilt(handbuilt_jsonfn)
  return hb['garbage']

def convert_adjlist_to_adjmatrix(adjlist):
  # Convert keys from strings to ints.
  adjlist = {int(K): V for (K, V) in adjlist.items()}

  all_children = [child for C in adjlist.values() for child in C]
  root = 0
  assert root not in all_children and root in adjlist.keys()

  N = max(all_children) + 1
  adjm = np.eye(N)

  for parent in adjlist.keys():
    children = adjlist[parent]
    adjm[parent, children] = 1

  return adjm

def load_tree(handbuilt_jsonfn):
  hb = _load_handbuilt(handbuilt_jsonfn)
  struct = hb['structure']

  children = [node for C in struct.values() for node in C]
  _check_dupes(children)
  # Node 0 should be root, and so should not appear in children list.
  assert set(children) == set(range(1, max(children) + 1))
  assert len(children) + 1 == len(hb['clusters'])

  adjm = convert_adjlist_to_adjmatrix(struct)
  return adjm
