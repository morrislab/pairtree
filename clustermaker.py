import numpy as np
import pairwise
import common

def _check_clusters(variants, clusters, garbage):
  for C in clusters:
      assert len(C) > 0

  vids = common.extract_vids(variants)
  clustered = [child for C in clusters for child in C]
  garbage = set(garbage)
  clustered = set(clustered)
  assert len(clustered & garbage) == 0
  assert set(vids) == (clustered | garbage)

def use_pre_existing(variants, mutrel, prior, parallel, clusters, garbage):
  supervars = make_cluster_supervars(clusters, variants)
  clust_posterior, clust_evidence = pairwise.calc_posterior(supervars, prior, 'supervariant', parallel)
  _check_clusters(variants, clusters, garbage)
  return (supervars, clust_posterior, clust_evidence, clusters, garbage)

def cluster_and_discard_garbage(variants, mutrel, prior, parallel):
  raise Exception('Not implemented yet')

def make_cluster_supervars(clusters, variants):
  supervars = {}

  for cidx, cluster in enumerate(clusters):
    assert len(cluster) > 0
    cvars = [variants[vid] for vid in cluster]

    cluster_total_reads = np.array([V['total_reads'] for V in cvars])
    cluster_var_reads = np.array([V['var_reads'] for V in cvars])
    # Correct for sex variants.
    omega_v = np.array([V['omega_v'] for V in cvars])[:,np.newaxis]
    cluster_var_reads = np.round(cluster_var_reads / (2*omega_v))

    S_name = 'S%s' % len(supervars)
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

    supervars[S_name] = S

  return supervars

def make_superclusters(supervars):
  N = len(supervars)
  svids = common.extract_vids(supervars)
  superclusters = [[S] for S in svids]
  return superclusters
