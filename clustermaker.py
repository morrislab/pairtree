import numpy as np
import pairwise

def _check_clusters(variants, clusters, garbage):
  for cidx, C in enumerate(clusters):
    if cidx == 0:
      assert len(C) == 0
    else:
      assert len(C) > 0

  vids = [int(V['id'][1:]) for V in variants]
  clustered = [child for C in clusters for child in C]
  garbage = set(garbage)
  clustered = set(clustered)
  assert len(clustered & garbage) == 0
  assert set(vids) == clustered | garbage == set(range(len(variants)))

def use_pre_existing(variants, mutrel, prior, parallel, clusters, garbage):
  vidmap = {V['id']: vidx for vidx, V in enumerate(variants)}
  # Convert from vids ("s0", "s1", "s3" -- gaps are possible, even though they
  # shouldn't occur, which we'll enforce elsewhere) to contiguous integer
  # indices.
  clusters = [[vidmap[vid] for vid in clust] for clust in clusters]
  garbage = [vidmap[vid] for vid in garbage]

  supervars = _make_cluster_supervars(clusters, variants)
  clust_posterior, clust_evidence = pairwise.calc_posterior(supervars, prior, 'supervariant', parallel)

  _check_clusters(variants, clusters, garbage)
  return (supervars, clust_posterior, clusters, garbage)

def cluster_and_discard_garbage(variants, mutrel, prior, parallel):
  raise Exception('Not implemented yet')

def _make_cluster_supervars(clusters, variants):
  supervars = []

  for cidx, cluster in enumerate(clusters):
    if len(cluster) == 0:
      continue
    cvars = [variants[vidx] for vidx in cluster]

    cluster_total_reads = np.array([V['total_reads'] for V in cvars])
    cluster_var_reads = np.array([V['var_reads'] for V in cvars])
    # Correct for sex variants.
    omega_v = np.array([V['omega_v'] for V in cvars])[:,np.newaxis]
    cluster_var_reads = np.round(cluster_var_reads / (2*omega_v))

    S_name = 'C%s' % len(supervars)
    S = {
      'id': S_name,
      'name': S_name,
      'chrom': None,
      'pos': None,
      'cluster': cidx,
      'omega_v': 0.5,
      'var_reads': np.sum(cluster_var_reads, axis=0, dtype=np.int),
      'total_reads': np.sum(cluster_total_reads, axis=0, dtype=np.int),
    }
    S['ref_reads'] = S['total_reads'] - S['var_reads']
    S['vaf'] = S['var_reads'] / S['total_reads']

    supervars.append(S)

  return supervars

def make_superclusters(supervars):
  N = len(supervars)
  superclusters = [[C] for C in range(N)]
  # Add empty initial cluster, which serves as tree root.
  superclusters.insert(0, [])
  return superclusters
