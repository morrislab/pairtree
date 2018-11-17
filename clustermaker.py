import numpy as np
import pairwise
import common
from common import Models, Mutrel, debug
from tqdm import tqdm

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
  clust_posterior, clust_evidence = pairwise.calc_posterior(supervars, prior, rel_type='supervariant', parallel=parallel)
  _check_clusters(variants, clusters, garbage)
  return (supervars, clust_posterior, clust_evidence, clusters, garbage)

def _minkowski_dist(A, B, p=2):
  return np.sum(np.abs(A - B)**p)**(1/p)

def _calc_row_dist(mat):
  M = len(mat)
  dist = np.zeros((M, M))
  for I in range(M):
    for J in range(I):
      dist[I,J] = _minkowski_dist(mat[I], mat[J])
      dist[J,I] = dist[I,J]
  return dist

def _remove_rowcol(arr, indices):
  '''Remove rows and columns at `indices`.'''
  # Using a mask requires only creating a new array once (and then presumably
  # another array to apply the mask). Calling np.delete() multiple times would
  # creat separate arrays for every element in `indices`.
  shape = list(arr.shape)
  # Must convert to list *before* array, as indices is often a set, and
  # converting a set directly to an array yields a dtype of `object` that
  # contains the set. It's really weird.
  indices = np.array(list(indices))
  # Only check the first two dimensions, as we ignore all others.
  assert np.all(0 <= indices) and np.all(indices < min(shape[:2]))

  for axis in (0, 1):
    arr = np.delete(arr, indices, axis=axis)
    shape[axis] -= len(indices)

  assert np.array_equal(arr.shape, shape)
  return arr

def _remove_variants(mutrel, vidxs):
  # Make set for efficient `in`.
  vidxs = set(vidxs)
  return Mutrel(
    vids = [vid for vidx, vid in enumerate(mutrel.vids) if vidx not in vidxs],
    rels = _remove_rowcol(mutrel.rels, vidxs),
  )

def _cluster_variants(variants, mutrel_posterior, mutrel_evidence, prior, parallel):
  # Don't change original variants. We store all our intermediate supervariants
  # in the copy we create.
  variants = dict(variants)
  # Each variant begins in its own cluster.
  clusters = [[vid] for vid in mutrel_posterior.vids]
  # Ensure no name collisions by always increment SV names.
  svidx_max = len(variants)

  with tqdm(desc='Clustering variants', unit='iteration', dynamic_ncols=True, miniters=1, disable=None) as pbar:
    while True:
      pbar.update()
      M = len(mutrel_posterior.rels)
      cocluster_rel = np.argmax(mutrel_posterior.rels, axis=2) == Models.cocluster
      # Must ignore diagonal, which will always want to cocluster.
      cocluster_rel[np.tril_indices(M)] = False
      if not np.any(cocluster_rel):
        break

      debug('')
      debug('cocluster argmax = ', np.sum(np.argmax(mutrel_posterior.rels, axis=2) == Models.cocluster) / mutrel_posterior.rels.size)
      debug(np.argmax(mutrel_posterior.rels, axis=2))
      assert np.array_equal(mutrel_posterior.vids, mutrel_evidence.vids)
      assert mutrel_posterior.rels.shape == mutrel_evidence.rels.shape
      # We have four tasks:
      #
      # 0. Pick variants to merge into new cluster
      # 1. Update the clustering
      # 2. Create a new supervariant
      # 3. Compute the mutrels for that supervariant relative to other variants

      # 0. Pick variants to merge into new cluster
      cocluster_prob = mutrel_posterior.rels[:,:,Models.cocluster]
      row_dist = _calc_row_dist(cocluster_prob)
      # Distances on diagonal will always be zero. Plus, we want A < B, so fill
      # the lower triangular portion as well.
      row_dist[np.tril_indices(M)] = np.inf
      # A, B: indices such that row_dist[A,B] is the smallest
      # element in the array.
      #debug('row_dist', row_dist)
      A, B = np.unravel_index(np.argmin(row_dist, axis=None), row_dist.shape)
      debug('A =', A, 'B =', B, 'dist =', row_dist[A,B], 'post =', mutrel_posterior.rels[A,B])
      debug(np.vstack((mutrel_posterior.rels[A,:,Models.cocluster], mutrel_posterior.rels[B,:,Models.cocluster])))
      debug((np.sort(np.abs(variants[mutrel_posterior.vids[A]]['vaf'] - variants[mutrel_posterior.vids[B]]['vaf']))))
      debug(mutrel_posterior.rels[A] - mutrel_posterior.rels[B])
      # Since row_dist is upper triangular, we should have A < B.
      assert A < B
      if mutrel_posterior.rels[A,B,Models.cocluster] < 1 / len(Models._all):
        break

      # 1. Update the clustering
      clusters.append(clusters[A] + clusters[B])
      # Delete in reverse order, since B > A, and this will change indices.
      del clusters[B]
      del clusters[A]

      # 2. Create a new supervariant
      svid = 'S%s' % svidx_max
      svidx_max += 1
      supervar = _make_supervar(svid, (variants[mutrel_posterior.vids[A]], variants[mutrel_posterior.vids[B]]))
      assert svid not in variants
      variants[svid] = supervar

      # 3. Compute the mutrels for that supervariant relative to other variants
      mutrel_posterior = _remove_variants(mutrel_posterior, (A, B))
      mutrel_evidence = _remove_variants(mutrel_evidence, (A, B))
      mutrel_posterior, mutrel_evidence = pairwise.add_variant(svid, variants, mutrel_posterior, mutrel_evidence, prior=prior, parallel=parallel)

  supervars, mutrel_posterior, mutrel_evidence, clusters = _sort_clusters_by_vaf(variants, mutrel_posterior, mutrel_evidence, clusters)
  return (supervars, mutrel_posterior, mutrel_evidence, clusters)

def _sort_clusters_by_vaf(variants, mutrel_posterior, mutrel_evidence, clusters):
  supervars = [variants[svid] for svid in mutrel_posterior.vids]
  sv_vaf = np.array([S['vaf'] for S in supervars])
  mean_vaf = np.mean(sv_vaf, axis=1)
  order = np.argsort(-mean_vaf)
  assert len(order) == len(mutrel_posterior.vids) == len(mutrel_evidence.vids) == len(clusters)

  clusters = [clusters[idx] for idx in order]
  new_vids = [mutrel_posterior.vids[idx] for idx in order]
  mutrel_posterior = mutrel_posterior._replace(vids=new_vids)
  mutrel_evidence = mutrel_evidence._replace(vids=new_vids)

  # Collect supervars into single dict. Supervar `S<i>` should correspond to
  # cluster `i`.
  supervars = {}
  for vidx, V in enumerate(mutrel_posterior.vids):
    assert mutrel_evidence.vids[vidx] == V
    # Copy original variant -- if the variant is the lone variant in a cluster,
    # it won't have been put in an associated supervariant, and so we want to
    # ensure we avoid modifying it (specifically, we want to avoid changing its
    # name & ID).
    var = dict(variants[V])
    del variants[V]
    S_name = 'S%s' % vidx
    # Ensure we've never used this name before (though it shouldn't matter if
    # we have).
    assert S_name not in variants
    var['id'] = var['name'] = S_name
    supervars[S_name] = variants[S_name] = var
    for mutrel in (mutrel_posterior, mutrel_evidence):
      mutrel.vids[vidx] = S_name

  return supervars, mutrel_posterior, mutrel_evidence, clusters

def _discard_garbage(mutrel_posterior, mutrel_evidence):
  garbage = []

  with tqdm(desc='Discarding garbage', unit='variant', dynamic_ncols=True, miniters=1, disable=None) as pbar:
    while True:
      garbage_pairs = np.argmax(mutrel_posterior.rels, axis=2) == Models.garbage
      if not np.any(garbage_pairs):
        break
      num_garbage = np.sum(garbage_pairs, axis=1)
      most_garbage = np.argmax(num_garbage)
      garbage.append(mutrel_posterior.vids[most_garbage])
      mutrel_posterior = _remove_variants(mutrel_posterior, (most_garbage,))
      mutrel_evidence = _remove_variants(mutrel_evidence, (most_garbage,))
      pbar.update()

  return garbage

def cluster_and_discard_garbage(variants, mutrel_posterior, mutrel_evidence, prior, parallel):
  # Copy mutrels so we don't modify them.
  mutrel_posterior = Mutrel(vids=list(mutrel_posterior.vids), rels=np.copy(mutrel_posterior.rels))
  mutrel_evidence = Mutrel(vids=list(mutrel_evidence.vids), rels=np.copy(mutrel_evidence.rels))

  # As a side effect, _discard_garbage() will remove the garbage variants
  # from the mutrels.
  garbage = _discard_garbage(mutrel_posterior, mutrel_evidence)
  supervars, clust_posterior, clust_evidence, clusters = _cluster_variants(variants, mutrel_posterior, mutrel_evidence, prior, parallel)
  return (supervars, clust_posterior, clust_evidence, clusters, garbage)

def _make_supervar(name, variants):
  cluster_total_reads = np.array([V['total_reads'] for V in variants])
  cluster_var_reads = np.array([V['var_reads'] for V in variants])
  # Correct for sex variants.
  omega_v = np.array([V['omega_v'] for V in variants])[:,np.newaxis]
  cluster_var_reads = np.round(cluster_var_reads / (2*omega_v))

  S = {
    'id': name,
    'name': name,
    'chrom': None,
    'pos': None,
    'omega_v': 0.5,
    'var_reads': np.sum(cluster_var_reads, axis=0, dtype=np.int),
    'total_reads': np.sum(cluster_total_reads, axis=0, dtype=np.int),
  }
  S['ref_reads'] = S['total_reads'] - S['var_reads']
  S['vaf'] = S['var_reads'] / S['total_reads']

  return S

def make_cluster_supervars(clusters, variants):
  supervars = {}

  for cidx, cluster in enumerate(clusters):
    assert len(cluster) > 0
    cvars = [variants[vid] for vid in cluster]
    S_name = 'S%s' % len(supervars)
    supervars[S_name] = _make_supervar(S_name, cvars)

  return supervars

def make_superclusters(supervars):
  N = len(supervars)
  svids = common.extract_vids(supervars)
  superclusters = [[S] for S in svids]
  return superclusters
