import csv
import numpy as np
from collections import namedtuple
import vaf_correcter

class Models:
  _all = ('garbage', 'cocluster', 'A_B', 'B_A', 'diff_branches')
for idx, M in enumerate(Models._all):
  setattr(Models, M, idx)

def parse_ssms(sampid, ssmfn):
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
        'mu_v': float(row['mu_v']),
      }
      variant['var_reads'] = variant['total_reads'] - variant['ref_reads']
      variants[row['id']] = variant

  vaf_correcter.correct_vafs(sampid, variants)
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

def make_cluster_supervars(clusters, variants):
  cluster_supervars = {}
  svid2svidx = {}

  for cidx, cluster in enumerate(clusters):
    if len(cluster) == 0:
      continue
    cvars = [variants['s%s' % vidx] for vidx in cluster]

    cluster_total_reads = np.array([V['total_reads'] for V in cvars])
    cluster_var_reads = np.array([V['var_reads'] for V in cvars])
    # Correct for sex variants.
    mu_v = np.array([V['mu_v'] for V in cvars])[:,np.newaxis]
    vaf_corrections = np.array([V['vaf_correction'] for V in cvars])[:,np.newaxis]
    cluster_var_reads = np.round(vaf_corrections * (cluster_var_reads / (2*(1 - mu_v))))

    S = {
      'gene': None,
      'id': 'C%s' % cidx,
      'name': 'C%s' % cidx,
      'chrom': None,
      'pos': None,
      'cluster': cidx,
      'mu_v': 0.499,
      'var_reads': np.sum(cluster_var_reads, axis=0),
      'total_reads': np.sum(cluster_total_reads, axis=0),
    }
    S['ref_reads'] = S['total_reads'] - S['var_reads']
    S['vaf'] = S['var_reads'] / S['total_reads']

    svid2svidx[S['id']] = len(cluster_supervars)
    cluster_supervars[S['id']] = S

  svidx2svid = {V: K for K, V in svid2svidx.items()}
  return (cluster_supervars, svid2svidx, svidx2svid)
