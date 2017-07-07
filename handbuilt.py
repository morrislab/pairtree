import json
import collections

def _check_dupes(coll):
  if len(set(coll)) != len(coll):
    raise Exception('Duplicates: %s' % [E for E, cnt in collections.Counter(coll).items() if cnt > 1])

def _load_handbuilt(jsonfn, allssms):
  with open(jsonfn) as F:
    H = json.load(F)

  clustered = [ssmid for C in H['clusters'] for ssmid in C]
  garbage = H['garbage']
  _check_dupes(clustered)
  _check_dupes(garbage)

  handbuilt_ssms = set(clustered) | set(garbage)
  handbuilt_ssms = set(['s%s' % S for S in handbuilt_ssms])

  allssms = set(allssms)
  if allssms  != handbuilt_ssms:
    raise Exception('SSMs missing from handbuilt: %s' % str(allssms - handbuilt_ssms))
    raise Exception('Extraneous SSMs in handbuilt: %s' % str(handbuilt_ssms - allssms))

  struct = H['structure']
  children = [node for C in struct.values() for node in C]
  print(children)
  _check_dupes(children)
  # Node 0 should be root, and so should not appear in children list.
  assert set(children) == set(range(1, max(children) + 1))

  return H

def cluster_variants(handbuilt_jsonfn, variants):
  hb = _load_handbuilt(handbuilt_jsonfn, variants.keys())
  clusters = hb['clusters']
  assert clusters[0] == []
  return clusters[1:]
