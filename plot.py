import csv
import argparse
import json
import numpy as np
import sklearn.cluster
import colorlover as cl
from common import parse_ssms, Models
from vaf_plotter import plot_vaf_matrix
from collections import defaultdict
import tree_sampler
import tree_builder
import json_writer
import phi_fitter

np.set_printoptions(threshold=np.nan)
np.random.seed(1)

def create_matrix(model, model_probs, var_names):
  N = len(var_names)
  mat = np.zeros((N, N))
  for vids, P in model_probs.items():
    assert 0 <= P <= 1
    vidx1, vidx2 = [int(V[1:]) for V in vids.split(',')]
    mat[(vidx1, vidx2)] = P

  if model in ('cocluster', 'diff_branches'):
    # These should be symmetric.
    assert np.allclose(mat, mat.T)

  return mat

def _reorder_row(R, idxs):
  return np.array([R[idx] for idx in idxs])

def cluster_rows(mat):
  affprop = sklearn.cluster.AffinityPropagation(damping=0.5)
  affprop.fit(mat)
  labels = affprop.predict(mat)

  N = len(mat)
  annotated = list(zip(mat, labels, range(N)))
  # Sort by cluster, then by row index.
  annotated = sorted(annotated, key = lambda R: R[1:])
  rows, idxs = zip(*[(R, I) for R, L, I in annotated])
  rows = [_reorder_row(R, idxs) for R in rows]

  return (np.array(rows), idxs)

def make_colour_from_intensity(intensity):
  val = np.int(np.round(255*(1 - intensity)))
  return 'rgb(255, %s, %s)' % (2*(val,))

def make_colour_from_category(cat):
  scalenum = len(Models._all)
  assert 0 <= cat < scalenum
  scale = cl.scales[str(scalenum)]['qual']['Set1']
  return scale[cat]

def make_colour_matrix(vals, colour_maker):
  N, M = vals.shape
  colours = N*[None]

  for i in range(N):
    colours[i] = M*[None]
    for j in range(M):
      colours[i][j] = colour_maker(vals[i][j])

  return colours

def make_table_row(entries, visibilities, colours):
  assert len(entries) == len(visibilities)
  entries = [E if V else ('<span>%s</span>' % E) for E, V in zip(entries, visibilities)]
  return '<tr>' + ''.join(['<td style="background-color: %s">%s</td>' % (C, E) for E, C in zip(entries, colours)]) + '</tr>'

def write_header(sampid, extra, outf):
  print('<script type="text/javascript" src="https://code.jquery.com/jquery-3.2.1.slim.min.js"></script>', file=outf)
  print('<script type="text/javascript" src="highlight_table_labels.js"></script>', file=outf)
  print('<link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet">', file=outf)
  print('<h1>%s (%s)</h1>' % (sampid, extra), file=outf)
  print('<style>td, table { padding: 5px; margin: 0; border-collapse: collapse; } span { visibility: hidden; } td:hover > span { visibility: visible; } .highlighted { background-color: black !important; color: white; font-weight: bold; }</style>', file=outf)

def write_table(model, mat, labels, colours, outf):
  print('<h2>%s</h2>' % model, file=outf)
  print('<table class="matrix"><thead></thead><tbody>', file=outf)

  N = len(mat)

  entries        = [''] + labels
  visibility     = len(entries)*[True]
  header_colours = len(entries)*['transparent']
  print(make_table_row(entries, visibility, header_colours), file=outf)

  for label, row, row_colours in zip(labels, mat, colours):
    entries     = [label]           + ['%.2f' % P for P in row]
    visibility  = [True]          + N*[False]
    row_colours = ['transparent'] + row_colours
    print(make_table_row(entries, visibility, row_colours), file=outf)

  print('</tbody></table>', file=outf)

def plot_individual(model_probs, should_cluster, outf):
  for model in Models._all:
    mat = create_matrix(model, model_probs['model_probs'][model], model_probs['var_names'])
    if should_cluster:
      mat, ssmidxs = cluster_rows(mat)
      if model in ('cocluster', 'diff_branches'):
        # These should be symmetric.
        assert np.allclose(mat, mat.T)
    else:
      ssmidxs = list(range(len(mat)))
    colours = make_colour_matrix(mat, make_colour_from_intensity)
    write_table(model, mat, ['s%s' % I for I in ssmidxs], colours, outf)

def write_legend(outf):
  print('<br><table class="table"><tbody>', file=outf)
  for midx, model in enumerate(Models._all):
    colour = make_colour_from_category(midx)
    print('<tr><td style="background-color: %s">%s</td></tr>' % (colour, model), file=outf)
  print('</tbody></table>', file=outf)

def write_cluster_map(cmap, outf):
  print('<br><table class="table"><thead><tr><th>Cluster</th><th>SSMs</th></tr></thead><tbody>', file=outf)
  for cidx, ssmidxs in enumerate(cmap):
    ssmidxs = ', '.join(['s%s' % I for I in sorted(ssmidxs)])
    print('<tr><td style="font-weight: bold">C%s</td><td>%s</td></tr>' % (cidx, ssmidxs), file=outf)
  print('</tbody></table>', file=outf)

def create_model_prob_tensor(model_probs):
  M = len(model_probs['var_names'])
  num_models = len(Models._all)
  tensor = np.zeros((M, M, num_models))
  for midx, mdl in enumerate(Models._all):
    tensor[:,:,midx] = create_matrix(mdl, model_probs['model_probs'][mdl], model_probs['var_names'])
  return tensor

def calc_relations(model_probs):
  M = len(model_probs)
  relations = np.argmax(model_probs, axis=2)
  assert relations.shape == (M, M)
  return relations

def plot_relations(relations, should_cluster, outf):
  if should_cluster:
    relations, ssmidxs = cluster_rows(relations)
  else:
    ssmidxs = list(range(len(relations)))

  colours = make_colour_matrix(relations, make_colour_from_category)
  write_table('relations', relations, ['s%s' % I for I in ssmidxs], colours, outf)

def combine_identical(mat):
  '''Combine identical rows & columns.'''
  rowmap = {}
  retained_idxs = []

  for idx, row in enumerate(mat):
    row = tuple(row)
    if row in rowmap:
      rowmap[row].append(idx)
    else:
      rowmap[row] = [idx]
      retained_idxs.append(idx)

  idxmap = []
  for row in mat[retained_idxs,:]:
    # idxmap lists what rows of `collapsed` correspond to what rows of `mat`
    # (one-to-many relation).
    idxmap.append(rowmap[tuple(row)])
  collapsed = mat[retained_idxs,:][:,retained_idxs]
  return (collapsed, idxmap)

def combine_identical_clusters(relations, clusters):
  combined, idxmap = combine_identical(relations)
  assert len(combined) <= len(relations)

  # idxmap maps newidx -> oldidxs. Reverse this.
  revidxmap = {}
  for newidx, oldidxs in enumerate(idxmap):
    for oldidx in oldidxs:
      assert oldidx not in revidxmap
      revidxmap[oldidx] = newidx

  new_clusters = defaultdict(list)
  for old_cidx, old_cluster in enumerate(clusters):
    new_cidx = revidxmap[old_cidx]
    new_clusters[new_cidx] += old_cluster

  assert set(new_clusters.keys()) == set(range(len(combined)))
  # Convert dictionary -> list.
  new_clusters = [sorted(new_clusters[idx]) for idx in sorted(new_clusters.keys())]

  return (combined, new_clusters)

def remove_small_clusters(mat, clusters, threshold=1):
  assert len(clusters) == len(mat)

  N = len(mat)
  to_remove = set([idx for idx, C in enumerate(clusters) if len(C) <= threshold])
  to_keep = [idx for idx in range(N) if idx not in to_remove]
  assert len(to_remove) + len(to_keep) == N

  filtered_mat = mat[to_keep,:][:,to_keep]
  filtered_clusters = [C for idx, C in enumerate(clusters) if idx not in to_remove]
  assert len(filtered_clusters) == len(to_keep)

  return (filtered_mat, filtered_clusters)

def cluster_relations(relations, remove_small):
  sorted_idxs = toposort(relations)
  # Sort rows by toposorted indexes.
  relations = np.array([_reorder_row(relations[idx], sorted_idxs) for idx in sorted_idxs])
  clusters = [[I] for I in sorted_idxs]

  relations, clusters = combine_identical_clusters(relations, clusters)
  if remove_small:
    relations, clusters = remove_small_clusters(relations, clusters)
    relations, clusters = combine_identical_clusters(relations, clusters)

  return (relations, clusters)

def plot_relations_toposort(relations, clusters, suffix, outf):
  colours = make_colour_matrix(relations, make_colour_from_category)
  labels = ['C%s' % I for I in range(len(clusters))]
  write_table('relations_toposort_%s' % suffix, relations, labels, colours, outf)
  write_cluster_map(clusters, outf)

  return clusters

def extract_B_A_rels(relations):
  ssmidxs = list(range(len(relations)))
  A_B_rels, B_A_rels = set(), set()
  B_A_adjlist = {}

  for sidx1 in ssmidxs:
    B_A_adjlist[sidx1] = set()

    for sidx2 in ssmidxs:
      rel = relations[sidx1,sidx2]
      if rel == Models.A_B:
        A_B_rels.add((sidx1, sidx2))
      elif rel == Models.B_A:
        B_A_rels.add((sidx1, sidx2))
        B_A_adjlist[sidx1].add(sidx2)

  assert len(A_B_rels) == len(B_A_rels)
  reversed_B_A_rels = set([(S2, S1) for (S1, S2) in B_A_rels])
  assert reversed_B_A_rels == A_B_rels

  return B_A_adjlist

def toposort(relations):
  # Algorithm taken from https://en.wikipedia.org/w/index.php?title=Topological_sorting&oldid=779516160#Kahn.27s_algorithm.
  B_A_rels = extract_B_A_rels(relations)
  topo_sorted = []
  indeg_zero = [vert for vert, ancs in B_A_rels.items() if len(ancs) == 0]

  while len(indeg_zero) > 0:
    indeg_zero.sort(reverse=True)
    node = indeg_zero.pop()
    topo_sorted.append(node)
    for child, ancs in B_A_rels.items():
      if node not in ancs:
        continue
      ancs.remove(node)
      if len(ancs) == 0:
        indeg_zero.append(child)

  assert set([len(ancs) for ancs in B_A_rels.values()]) == set([0]), 'Graph has cycle'
  if len(topo_sorted) != len(B_A_rels):
    raise Exception('Graph has cycle')
  return topo_sorted

def euclid_dist(A, B):
  return np.sqrt(np.sum((B - A)**2))

def find_closest(tensor, target):
  # Find matrix in tensor that's closest to `target` in Euclidean distance.
  min_dist = float('inf')
  best_idx = None

  for idx, mat in enumerate(tensor):
    dist = euclid_dist(mat, target)
    if dist < min_dist:
      min_dist = dist
      best_idx = idx

  assert best_idx is not None
  return best_idx

def assign_missing(clusters, model_probs):
  K = len(clusters)
  M = len(model_probs)

  already_assigned = set([midx for cluster in clusters for midx in cluster])
  missing = set(range(M)) - already_assigned

  # Mean relations of members.
  cluster_rels = np.array([np.mean(np.array([model_probs[midx] for midx in cluster]), axis=0) for cluster in clusters])
  assert cluster_rels.shape == (K, M, len(Models._all))

  for midx in missing:
    closest_cluster_idx = find_closest(cluster_rels, model_probs[midx])
    clusters[closest_cluster_idx].append(midx)

  for cluster in clusters:
    cluster.sort()

  assigned = set([midx for cluster in clusters for midx in cluster])
  assert len(assigned) == M
  assert  assigned == already_assigned | missing
  return clusters

def find_root_node(adj):
  K = len(adj)
  assert adj.shape == (K, K)
  assert np.array_equal(np.diag(adj), np.ones(K))
  assert np.sum(adj) == K + (K - 1)
  assert np.array_equal(np.sort(np.sum(adj, axis=0)), np.array([1] + (K - 1)*[2]))
  num_parents = np.sum(adj, axis=0) - 1
  root = np.where(num_parents == 0)[0][0]
  return int(root) # Must convert from NumPy int.

def add_normal_root(adj):
  old_root = find_root_node(adj)
  for axis in (0, 1):
    # Insert additional first row & column of zeroes.
    adj = np.insert(adj, 0, 0, axis=axis)
  # Add root non-cancerous node.
  adj[0,0] = adj[0,old_root+1] = 1
  return adj

def plot(sampid, model_probs, output_type, ssmfn, paramsfn, spreadsheetfn, outfn, treesummfn, mutlistfn):
  should_cluster = not (output_type == 'unclustered')
  model_probs_tensor = create_model_prob_tensor(model_probs)
  relations = calc_relations(model_probs_tensor)
  variants = parse_ssms(ssmfn)

  with open(outfn, 'w') as outf:
    write_header(sampid, output_type, outf)
    write_legend(outf)
    if output_type != 'condensed':
      plot_individual(model_probs, should_cluster, outf)
      plot_relations(relations, should_cluster, outf)
    for remove_small in (False, True):
      clustered_relations, clusters = cluster_relations(relations, remove_small)

      if remove_small:
        assign_missing(clusters, model_probs_tensor)
        sampled_adjm, sampled_llh = tree_sampler.sample_trees(model_probs_tensor, clusters)

        handbuilt_adjm = tree_builder.make_adj(clustered_relations)
        sampled_adjm.insert(0, handbuilt_adjm)
        sampled_llh.insert(0, 0)

        sampled_adjm = [add_normal_root(adj) for adj in sampled_adjm]
        clusters.insert(0, [])

        phi = phi_fitter.fit_phis(sampled_adjm[0], clusters, variants)
        print('phi', phi)

        json_writer.write_json(sampid, variants, clusters, sampled_adjm, sampled_llh, treesummfn, mutlistfn)

      suffix = remove_small and 'small_excluded' or 'small_included'
      plot_relations_toposort(clustered_relations, clusters, suffix, outf)
      plot_vaf_matrix(clusters, variants, paramsfn, spreadsheetfn, outf)

def load_model_probs(model_probs_fn):
  with open(model_probs_fn) as F:
    return json.load(F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--output-type', dest='output_type', choices=('clustered', 'unclustered', 'condensed'), default='clustered')
  parser.add_argument('sampid')
  parser.add_argument('model_probs_fn')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('spreadsheet_fn')
  parser.add_argument('out_fn')
  parser.add_argument('treesumm_fn')
  parser.add_argument('mutlist_fn')
  args = parser.parse_args()

  model_probs = load_model_probs(args.model_probs_fn)
  plot(args.sampid, model_probs, args.output_type, args.ssm_fn, args.params_fn, args.spreadsheet_fn, args.out_fn, args.treesumm_fn, args.mutlist_fn)

main()
