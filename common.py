import csv
import numpy as np
from collections import namedtuple

class Models:
  _all = ('cocluster', 'A_B', 'B_A', 'diff_branches')
for idx, M in enumerate(Models._all):
  setattr(Models, M, idx)

def parse_ssms(ssmfn):
  vaf = []
  ssm_ids = []
  var_names = []
  variants = {}

  with open(ssmfn) as F:
    reader = csv.DictReader(F, delimiter='\t')
    for row in reader:
      if False and len(variants) >= 3:
        break
      variant = {
        'id': row['id'],
        'name': row['gene'],
        'ref_reads': np.array([float(V) for V in row['a'].split(',')]),
        'total_reads': np.array([float(V) for V in row['d'].split(',')]),
      }
      variant['var_reads'] = variant['total_reads'] - variant['ref_reads']
      variants[row['id']] = variant

  return variants

def make_ancestral_from_adj(adj):
  K = len(adj)
  Z = np.zeros((K,K))

  def _find_desc(I, vec):
    # Base case: if I have no children, my ancestor vec is just myself.
    if np.sum(vec) == 0:
      return vec
    else:
      children = np.array([_find_desc(idx, adj[idx]) for (idx, val) in enumerate(vec) if idx != I and val == 1])
      self_and_child = vec + np.sum(children, axis=0)
      self_and_child[self_and_child > 1] = 1
      return self_and_child

  for k in range(K):
    # If we know `adj` is topologically sorted, we can reduce the complexity of
    # this -- we would start at leaves and work our way upward, eliminating
    # need for recursive DFS. But since we don't expect `K` to be large, we can
    # write more general version that works for non-toposorted trees.
    Z[k] = _find_desc(k, adj[k])

  return Z
