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

def _merge_clusters(to_merge, clusters, mutrel_evidence, prior):
  debug('to_merge', [[mutrel_evidence.vids[I] for I in C] for C in to_merge])
  # 1. Update the clustering
  to_remove = []
  for group in to_merge:
    new_cluster = [member for C in group for member in clusters[C]]
    clusters.append(new_cluster)
    to_remove += group
  # Remove only after cluster addition complete to avoid changing indices.
  for C in sorted(to_remove, reverse=True):
    del clusters[C]

  # 2. Update the mutrels
  mutrel_posterior, mutrel_evidence = pairwise.merge_variants(to_merge, mutrel_evidence, prior)

  M = len(clusters)
  assert mutrel_posterior.rels.shape == mutrel_evidence.rels.shape == (M, M, len(Models._all))
  assert len(mutrel_posterior.vids) == len(mutrel_evidence.vids) == M
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

def _sort_clusters_by_vaf(clusters, variants, mutrel_posterior, mutrel_evidence):
  supervars = make_cluster_supervars(clusters, variants)
  svids = common.extract_vids(supervars)

  sv_vaf = np.array([supervars[S]['vaf'] for S in svids])
  mean_vaf = np.mean(sv_vaf, axis=1)
  order = np.argsort(-mean_vaf)
  assert len(order) == len(mutrel_posterior.vids) == len(mutrel_evidence.vids) == len(clusters)

  clusters = [clusters[idx] for idx in order]
  # Recreate supervars so that the SV IDs are in order of descending cluster VAF.
  supervars = make_cluster_supervars(clusters, variants)
  # Recreate `svids` as well to be safe, even though it should be unchanged.
  svids = common.extract_vids(supervars)

  mutrel_posterior = Mutrel(
    vids = svids,
    rels = _reorder_matrix(mutrel_posterior.rels, order),
  )
  mutrel_evidence = Mutrel(
    vids = svids,
    rels = _reorder_matrix(mutrel_evidence.rels, order),
  )

  return (clusters, supervars, mutrel_posterior, mutrel_evidence)

def _discard_garbage(clusters, mutrel_posterior, mutrel_evidence):
  garbage = []

  while True:
    M = len(clusters)
    assert len(mutrel_posterior.vids) == M
    assert mutrel_posterior.rels.shape == (M, M, len(Models._all))

    garbage_pairs = np.argmax(mutrel_posterior.rels, axis=2) == Models.garbage
    if not np.any(garbage_pairs):
      break
    num_garbage = np.sum(garbage_pairs, axis=1)
    most_garbage = np.argmax(num_garbage)
    debug('most_garbage', mutrel_posterior.vids[most_garbage], clusters[most_garbage])

    garbage += clusters[most_garbage]
    del clusters[most_garbage]
    mutrel_posterior = pairwise.remove_variants(mutrel_posterior, (most_garbage,))
    mutrel_evidence = pairwise.remove_variants(mutrel_evidence, (most_garbage,))

  return (clusters, garbage, mutrel_posterior, mutrel_evidence)

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

def _iterate_clustering(selector, desc, clusters, clust_posterior, clust_evidence, prior):
  # Do initial round of garbage rejection.
  clusters, garbage, clust_posterior, clust_evidence = _discard_garbage(clusters, clust_posterior, clust_evidence)

  with progressbar(desc=desc, unit='step', dynamic_ncols=True, miniters=1) as pbar:
    while True:
      #_plot(clust_posterior, clusters, variants, garbage)
      pbar.update()
      to_merge = selector(clust_posterior)
      if len(to_merge) == 0:
        break
      clusters, clust_posterior, clust_evidence = _merge_clusters(to_merge, clusters, clust_evidence, prior)

      #_plot(clust_posterior, clusters, variants, garbage)
      pbar.update()
      clusters, garb_vids, clust_posterior, clust_evidence = _discard_garbage(clusters, clust_posterior, clust_evidence)
      garbage += garb_vids

  return (clusters, garbage, clust_posterior, clust_evidence)

def cluster_and_discard_garbage(variants, mutrel_posterior, mutrel_evidence, prior, parallel):
  assert np.all(mutrel_posterior.vids == mutrel_evidence.vids)
  # Copy mutrels so we don't modify them.
  clust_posterior = Mutrel(vids=list(mutrel_posterior.vids), rels=np.copy(mutrel_posterior.rels))
  clust_evidence = Mutrel(vids=list(mutrel_evidence.vids), rels=np.copy(mutrel_evidence.rels))

  # Each variant begins in its own cluster.
  clusters = [[vid] for vid in clust_posterior.vids]
  garbage = []
  for selector, desc in ((_find_identical_mle, 'Merging identical relations'), (_merge_until_no_pairwise_mle_remain, 'Merging similar relations')):
    clusters, G, clust_posterior, clust_evidence = _iterate_clustering(
      selector,
      desc,
      clusters,
      clust_posterior,
      clust_evidence,
      prior
    )
    garbage += G
    debug('selector=%s V=%s C=%s' % (selector.__name__, len(variants), len(clusters)))

  clusters, supervars, clust_posterior, clust_evidence = _sort_clusters_by_vaf(clusters, variants, clust_posterior, clust_evidence)
  # Note that the posterior and evidence for clusters *won't* be the same as if
  # we had computed them on the supervariants, since the evidence was just
  # computed by summing cluster members' evidences. Thus, recompute those now.
  clust_posterior, clust_evidence = pairwise.calc_posterior(supervars, prior, 'cluster', parallel)
  return (supervars, clust_posterior, clust_evidence, clusters, garbage)

def _make_supervar(name, variants):
  N = np.array([var['total_reads'] for var in variants])
  M, S = N.shape
  V = np.array([var['var_reads'] for var in variants])
  omega_v = np.array([var['omega_v'] for var in variants])

  # In creating a supervariant, we must rescale either `V` or `N` (the variant
  # or total read counts, respectively), as we're fixing the supervariant's
  # `omega_v`. We want the ratio `V / N*omega_v` to remain constant, as this
  # represents the "fraction of informative reads that are variant reads". That
  # is, `N * omega_v` indicates the number of reads that *could* be variant,
  # given what we know of the locus' copy-number state.
  #
  # As we will fix `omega_v` at 0.5 for the rescaled variant, and we want `V /
  # N*omega_v` to remain constant, we can either change `V` or `N` to maintain
  # the relationship `V / N*omega_v = V_hat / N_hat*omega_v`. We can solve this
  # equation two ways, depending on whether we modify `V` or `N`:
  #
  # 1. V_hat = V / 2*omega_v
  # 2. N_hat = 2*N*omega_v
  #
  # The advantage to approach 2 is that, when `omega_v` is small, it preserves
  # the binomial variance. We have `var = np(1 - p) =~ np = E[X]` when `p` is
  # small, meaning the variance scales with the mean when `p` is small. To
  # maintain the variance, we should try to maintain `np`, and so since we're
  # making `p` bigger, we should make `N` smaller, rather than changing `V`.
  #
  # Also, this means we can avoid any weird hackery when `omega_v = 0`, since
  # we don't have to add edge case checks to avoid division by zero.
  N_hat = 2*N*omega_v
  V_hat = np.minimum(V, N_hat)
  omega_v_hat = 0.5 * np.ones(S)

  svar = {
    'id': name,
    'name': name,
    'chrom': None,
    'pos': None,
    'omega_v': omega_v_hat,
    'var_reads': np.round(np.sum(V_hat, axis=0)).astype(np.int),
    'total_reads': np.round(np.sum(N_hat, axis=0)).astype(np.int),
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
