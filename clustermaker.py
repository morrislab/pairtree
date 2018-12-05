import numpy as np
import pairwise
import common
from common import Models, Mutrel, debug
from progressbar import progressbar

def _check_clusters(variants, clusters, garbage):
  for C in clusters:
      assert len(C) > 0

  vids = common.extract_vids(variants)
  clustered = [child for C in clusters for child in C]
  garbage = set(garbage)
  clustered = set(clustered)
  assert len(clustered & garbage) == 0
  assert set(vids) == (clustered | garbage)

def use_pre_existing(variants, prior, parallel, clusters, garbage):
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

def _find_identical_mle(mutrel_posterior):
  ml_rels = np.argmax(mutrel_posterior.rels, axis=2)
  row_dist = _calc_row_dist(ml_rels)
  # Distances on diagonal will always be zero. Plus, we want A < B, so fill
  # the lower triangular portion as well.
  row_dist[np.tril_indices(len(row_dist))] = np.inf
  if not np.any(np.isclose(0, row_dist)):
    return []

  # We have four tasks:
  #
  # 0. Pick variants to merge into new clusters
  # 1. Update the clustering
  # 2. Create a new supervariant
  # 3. Compute the mutrels for that supervariant relative to other variants

  # 0. Pick variants to merge into new clusters
  to_merge = []
  used_cidxs = set()
  for I in range(len(row_dist)):
    J = np.where(np.isclose(0, row_dist[I]))[0]
    group = set([I] + J.tolist())
    group -= used_cidxs
    if len(group) <= 1:
      continue
    to_merge.append(group)
    used_cidxs |= set(group)

  return to_merge

def _merge_until_no_pairwise_mle_remain(mutrel_posterior):
  row_dist = _calc_row_dist(mutrel_posterior.rels)
  row_dist[np.tril_indices(len(row_dist))] = np.inf

  threshold = 1 / len(Models._all)
  # Only want to merge pairs whose cluster probability exceeds threshold.
  noncoclustered = mutrel_posterior.rels[:,:,Models.cocluster] < threshold
  row_dist[noncoclustered] = np.inf
  if np.all(np.isinf(row_dist)):
    return []

  A, B = np.unravel_index(np.argmin(row_dist, axis=None), row_dist.shape)
  #debug('row_dist', row_dist)
  #debug('A =', A, 'B =', B, 'dist =', row_dist[A,B], 'post =', mutrel_posterior.rels[A,B])
  #debug(np.vstack((mutrel_posterior.rels[A,:,Models.cocluster], mutrel_posterior.rels[B,:,Models.cocluster])))
  #debug(mutrel_posterior.rels[A] - mutrel_posterior.rels[B])

  assert A < B
  return [(A, B)]

def _merge_clusters(to_merge, clusters, variants, mutrel_posterior, mutrel_evidence, prior, pbar, parallel):
  debug('to_merge', [[mutrel_posterior.vids[I] for I in C] for C in to_merge])
  # 1. Update the clustering
  to_remove = []
  for group in to_merge:
    new_cluster = [member for C in group for member in clusters[C]]
    clusters.append(new_cluster)
    to_remove += group
  # Remove only after cluster addition complete to avoid changing indices.
  for C in sorted(to_remove, reverse=True):
    del clusters[C]

  # 2. Create new supervariants
  new_svids = []
  svidx_max = len(variants)
  for group in to_merge:
    svid = 'S%s' % svidx_max
    svidx_max += 1
    new_svids.append(svid)
    supervar = _make_supervar(svid, [variants[mutrel_posterior.vids[I]] for I in group])
    assert svid not in variants
    variants[svid] = supervar

  # 3. Compute the mutrels for the new supervars
  mutrel_posterior = _remove_variants(mutrel_posterior, to_remove)
  mutrel_evidence = _remove_variants(mutrel_evidence, to_remove)
  mutrel_posterior, mutrel_evidence = pairwise.add_variants(new_svids, variants, mutrel_posterior, mutrel_evidence, prior=prior, pbar=pbar, parallel=parallel)

  return (clusters, mutrel_posterior, mutrel_evidence)

def _reorder_matrix(mat, order):
  # common.reorder_square_matrix does hierarchical clustering to determine
  # order. This is a much simpler and more elegant function that uses a
  # user-defined order.
  M = len(mat)
  assert mat.shape[:2] == (M, M)
  assert set(order) == set(range(M))

  mat = mat[order,:]
  mat = mat[:,order]
  return mat

def _sort_clusters_by_vaf(variants, mutrel_posterior, mutrel_evidence, clusters):
  supervars = [variants[svid] for svid in mutrel_posterior.vids]
  sv_vaf = np.array([S['vaf'] for S in supervars])
  mean_vaf = np.mean(sv_vaf, axis=1)
  order = np.argsort(-mean_vaf)
  assert len(order) == len(mutrel_posterior.vids) == len(mutrel_evidence.vids) == len(clusters)

  clusters = [clusters[idx] for idx in order]
  new_vids = [mutrel_posterior.vids[idx] for idx in order]

  mutrel_posterior = Mutrel(
    vids = new_vids,
    rels = _reorder_matrix(mutrel_posterior.rels, order),
  )
  mutrel_evidence = Mutrel(
    vids = new_vids,
    rels = _reorder_matrix(mutrel_evidence.rels, order),
  )

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

def _discard_garbage(clusters, mutrel_posterior, mutrel_evidence):
  garbage = []

  while True:
    M = len(clusters)
    assert len(mutrel_posterior.vids) == len(mutrel_evidence.vids) == M
    assert mutrel_posterior.rels.shape == mutrel_evidence.rels.shape == (M, M, len(Models._all))

    garbage_pairs = np.argmax(mutrel_posterior.rels, axis=2) == Models.garbage
    if not np.any(garbage_pairs):
      break
    num_garbage = np.sum(garbage_pairs, axis=1)
    most_garbage = np.argmax(num_garbage)
    debug('most_garbage', mutrel_posterior.vids[most_garbage], clusters[most_garbage])

    garbage += clusters[most_garbage]
    del clusters[most_garbage]
    mutrel_posterior = _remove_variants(mutrel_posterior, (most_garbage,))
    mutrel_evidence = _remove_variants(mutrel_evidence, (most_garbage,))

  return (clusters, mutrel_posterior, mutrel_evidence, garbage)

def _plot(mutrel_posterior, clusters, variants, garbage):
  if _plot.prefix is None:
    return
  import plotter
  import relation_plotter
  import vaf_plotter
  with open('%s_%s.html' % (_plot.prefix, _plot.idx), 'w') as outf:
    plotter.write_header(_plot.prefix, outf)
    _plot.idx += 1
    relation_plotter.plot_ml_relations(mutrel_posterior, outf)
    relation_plotter.plot_separate_relations(mutrel_posterior, outf)
    samps = ['Samp %s' % (sidx + 1) for sidx in range(len(next(iter(variants.values()))['var_reads']))]
    supervars = [variants[S] for S in mutrel_posterior.vids]
    for V in variants.values():
      V['chrom'] = V['pos'] = None
    vaf_plotter.plot_vaf_matrix(
      clusters,
      variants,
      supervars,
      garbage,
      None,
      samps,
      should_correct_vaf=True,
      outf=outf,
    )
_plot.idx = 1
_plot.prefix = None

def _iterate_clustering(selector, desc, variants, clusters, clust_posterior, clust_evidence, prior, parallel):
  # Do initial round of garbage rejection.
  clusters, clust_posterior, clust_evidence, garbage = _discard_garbage(clusters, clust_posterior, clust_evidence)

  with progressbar(desc=desc, unit='pair', dynamic_ncols=True, miniters=1) as pbar:
    while True:
      #_plot(clust_posterior, clusters, variants, garbage)
      pbar.update()
      to_merge = selector(clust_posterior)
      if len(to_merge) == 0:
        break
      clusters, clust_posterior, clust_evidence = _merge_clusters(to_merge, clusters, variants, clust_posterior, clust_evidence, prior, pbar, parallel)

      #_plot(clust_posterior, clusters, variants, garbage)
      pbar.update()
      clusters, clust_posterior, clust_evidence, garb_vids = _discard_garbage(clusters, clust_posterior, clust_evidence)
      garbage += garb_vids

  return (clusters, garbage, clust_posterior, clust_evidence)

def cluster_and_discard_garbage(variants, mutrel_posterior, mutrel_evidence, prior, parallel):
  # Copy mutrels so we don't modify them.
  clust_posterior = Mutrel(vids=list(mutrel_posterior.vids), rels=np.copy(mutrel_posterior.rels))
  clust_evidence = Mutrel(vids=list(mutrel_evidence.vids), rels=np.copy(mutrel_evidence.rels))
  # Don't change original variants. We store all our intermediate supervariants
  # in the copy we create.
  variants = dict(variants)

  # Each variant begins in its own cluster.
  clusters = [[vid] for vid in mutrel_posterior.vids]
  garbage = []
  for selector, desc in ((_find_identical_mle, 'Merging identical relations'), (_merge_until_no_pairwise_mle_remain, 'Merging similar relations')):
    clusters, G, clust_posterior, clust_evidence = _iterate_clustering(
      selector,
      desc,
      variants,
      clusters,
      clust_posterior,
      clust_evidence,
      prior,
      parallel
    )
    garbage += G
    debug('V=%s C=%s' % (len(variants), len(clusters)))

  supervars, clust_posterior, clust_evidence, clusters = _sort_clusters_by_vaf(variants, clust_posterior, clust_evidence, clusters)
  return (supervars, clust_posterior, clust_evidence, clusters, garbage)

def _make_supervar(name, variants):
  cluster_total_reads = np.array([V['total_reads'] for V in variants])
  cluster_var_reads = np.array([V['var_reads'] for V in variants])
  # Correct for sex variants.
  omega_v = np.array([V['omega_v'] for V in variants])
  M, S = omega_v.shape

  cluster_var_reads = np.round(cluster_var_reads / (2*omega_v))
  # Don't allow var read count to exceed total read count, which can happen
  # when `omega_v` is small.
  cluster_var_reads = np.minimum(cluster_var_reads, cluster_total_reads)

  svar = {
    'id': name,
    'name': name,
    'chrom': None,
    'pos': None,
    'omega_v': 0.5 * np.ones(S),
    'var_reads': np.sum(cluster_var_reads, axis=0, dtype=np.int),
    'total_reads': np.sum(cluster_total_reads, axis=0, dtype=np.int),
  }
  svar['ref_reads'] = svar['total_reads'] - svar['var_reads']
  svar['vaf'] = svar['var_reads'] / svar['total_reads']

  return svar

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
