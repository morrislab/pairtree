import csv
import argparse
import json
import numpy as np
import sklearn.cluster
import colorlover as cl
from common import parse_ssms, Models
import vaf_plotter
from collections import defaultdict
import tree_sampler
import tree_builder
import json_writer
import phi_fitter
import handbuilt

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

def make_table_row(entries, visibilities, colours, is_header=False):
  elem = 'th' if is_header else 'td'
  assert len(entries) == len(visibilities)
  entries = [E if V else ('<span>%s</span>' % E) for E, V in zip(entries, visibilities)]
  return '<tr>' + ''.join(['<%s style="background-color: %s">%s</%s>' % (elem, C, E, elem) for E, C in zip(entries, colours)]) + '</tr>'

def write_header(sampid, extra, outf):
  print('<script type="text/javascript" src="https://code.jquery.com/jquery-3.2.1.slim.min.js"></script>', file=outf)
  print('<script type="text/javascript" src="highlight_table_labels.js"></script>', file=outf)
  print('<link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet">', file=outf)
  print('<h1>%s (%s)</h1>' % (sampid, extra), file=outf)
  print('<style>td, th, table { padding: 5px; margin: 0; border-collapse: collapse; font-weight: normal; } span { visibility: hidden; } td:hover > span { visibility: visible; } .highlighted { background-color: black !important; color: white; }</style>', file=outf)

def write_table(model, mat, labels, colours, outf):
  print('<h2>%s</h2>' % model, file=outf)
  print('<table class="matrix"><thead>', file=outf)

  N = len(mat)

  entries        = [''] + labels
  visibility     = len(entries)*[True]
  header_colours = len(entries)*['transparent']
  print(make_table_row(entries, visibility, header_colours, is_header=True), file=outf)
  print('</thead><tbody>', file=outf)

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

def cluster_relations(relations):
  sorted_idxs = toposort(relations)
  # Sort rows by toposorted indexes.
  relations = np.array([_reorder_row(relations[idx], sorted_idxs) for idx in sorted_idxs])
  clusters = [[I] for I in sorted_idxs]

  relations, clusters = combine_identical_clusters(relations, clusters)
  return (relations, clusters)

def plot_relations_toposort(relations, clusters, outf):
  colours = make_colour_matrix(relations, make_colour_from_category)
  labels = ['C%s' % I for I in range(len(clusters))]
  write_table('relations_toposort', relations, labels, colours, outf)
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

def make_vaf_matrix(variants):
  vaf = np.array([V['var_reads'] / V['total_reads'] for V in variants.values()])
  return vaf

def cluster_variants(model_probs_tensor, vaf_matrix):
  M = len(model_probs_tensor)
  features = np.zeros((M, M*len(Models._all)))
  for idx, mat in enumerate(model_probs_tensor):
    features[idx] = np.ravel(mat)

  #features = vaf_matrix
  # Insert clustering algorithm here. Bleh.
  clusterer = hdbscan.HDBSCAN(min_cluster_size=2)
  labels = clusterer.fit_predict(vaf_matrix)

  clusters = [list() for _ in range(np.max(labels) + 1)]
  unclustered = []
  for idx, label in enumerate(labels):
    dest = clusters[label] if label >= 0 else unclustered
    dest.append(idx)

  return clusters

def make_trees(variants, model_probs, clusters):
  model_probs_tensor = create_model_prob_tensor(model_probs)

  #assign_missing(clusters, model_probs_tensor)
  sampled_adjm, sampled_llh = tree_sampler.sample_trees(model_probs_tensor, clusters, 1000)
  #sampled_adjm = [add_normal_root(adj) for adj in sampled_adjm]
  #clusters.insert(0, [])

  nsamples = len(list(variants.values())[0]['total_reads'])
  ntrees = len(sampled_adjm)
  phi = np.ones((ntrees, nsamples, len(clusters)))
  phi[-1,:,:] = phi_fitter.fit_phis(sampled_adjm[-1], clusters, variants)

  return (sampled_adjm, sampled_llh, phi)

def remove_garbage(garbage_ids, model_probs, variants, clusters):
  garbage_variants = {'g%s' % gidx: variants['s%s' % vidx] for gidx, vidx in enumerate(garbage_ids)}
  for gvid, var in garbage_variants.items():
    var['id'] = gvid

  varid_map = {}
  for varid in sorted(variants.keys(), key = lambda S: int(S[1:])):
    if int(varid[1:]) in garbage_ids:
      del variants[varid]
    else:
      varid_map[varid] = 's%s' % len(varid_map)

  remapped_model_probs = {}
  remapped_model_probs['models'] = model_probs['models']
  remapped_model_probs['var_names'] = {}
  for varid in model_probs['var_names'].keys():
    if int(varid[1:]) in garbage_ids:
      continue
    remapped_model_probs['var_names'][varid_map[varid]] = model_probs['var_names'][varid]

  for K in ('model_probs', 'model_evidence'):
    remapped_model_probs[K] = {}
    for model in model_probs[K].keys():
      remapped_model_probs[K][model] = {}
      for pair in list(model_probs[K][model].keys()):
        vidx1, vidx2 = [int(V[1:]) for V in pair.split(',')]
        if vidx1 in garbage_ids or vidx2 in garbage_ids:
          pass
        else:
          vid1new, vid2new = varid_map['s%s' % vidx1], varid_map['s%s' % vidx2]
          remapped_model_probs[K][model]['%s,%s' % (vid1new, vid2new)] = model_probs[K][model][pair]

  renamed_variants = {varid_map[S]: variants[S] for S in variants.keys()}
  for var in renamed_variants.values():
    var['id'] = varid_map[var['id']]

  for cidx, cluster in enumerate(clusters):
    clusters[cidx] = [int(varid_map['s%s' % S][1:]) for S in cluster]

  return (remapped_model_probs, renamed_variants, garbage_variants, clusters)

def plot(sampid, model_probs, output_type, ssmfn, paramsfn, spreadsheetfn, handbuiltfn, outfn, treesummfn, mutlistfn):
  variants = parse_ssms(ssmfn)
  garbage_ids = handbuilt.load_garbage(handbuiltfn)
  clusters = handbuilt.load_clusters(handbuiltfn, variants)
  model_probs, variants, garbage_variants, clusters = remove_garbage(garbage_ids, model_probs, variants, clusters)

  should_cluster = not (output_type == 'unclustered')
  model_probs_tensor = create_model_prob_tensor(model_probs)
  relations = calc_relations(model_probs_tensor)

  handbuilt_adjm = handbuilt.load_tree(handbuiltfn)
  handbuilt_mutrel = tree_sampler.make_mutrel_tensor_from_cluster_adj(handbuilt_adjm, clusters)
  handbuilt_llh = tree_sampler.calc_llh(model_probs_tensor, handbuilt_mutrel)

  #vaf = make_vaf_matrix(variants)
  #clusters = cluster_variants(model_probs_tensor, vaf)

  with open(outfn, 'w') as outf:
    write_header(sampid, output_type, outf)
    write_legend(outf)
    if output_type != 'condensed':
      plot_individual(model_probs, should_cluster, outf)
      plot_relations(relations, should_cluster, outf)

    #clustered_relations, clusters = cluster_relations(relations, remove_small)
    clustered_relations, _ = cluster_relations(relations)

    sampled_adjm, sampled_llh, phi = make_trees(variants, model_probs, clusters)
    sampled_adjm.insert(0, handbuilt_adjm)
    sampled_llh.insert(0, handbuilt_llh)
    phi = np.insert(phi, 0, phi[0,:,:], axis=0)

    # TODO: adjust this to use "super SSMs" generated by summing `a` and `d` for clustered mutations.
    #mle_adjm = tree_builder.make_adj(clustered_relations)
    #mle_adjm = add_normal_root(mle_adjm)
    #mle_mutrel = tree_sampler.make_mutrel_tensor_from_cluster_adj(mle_adjm, clusters)
    #sampled_adjm.insert(0, mle_adjm)
    #sampled_llh.insert(0, 0)

    json_writer.write_json(sampid, variants, clusters, sampled_adjm, sampled_llh, phi, treesummfn, mutlistfn)
    plot_relations_toposort(clustered_relations, clusters, outf)
    vaf_plotter.plot_vaf_matrix(clusters, variants, garbage_variants, paramsfn, spreadsheetfn, outf)

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
  parser.add_argument('handbuilt_fn')
  parser.add_argument('out_fn')
  parser.add_argument('treesumm_fn')
  parser.add_argument('mutlist_fn')
  args = parser.parse_args()

  model_probs = load_model_probs(args.model_probs_fn)
  plot(args.sampid, model_probs, args.output_type, args.ssm_fn, args.params_fn, args.spreadsheet_fn, args.handbuilt_fn, args.out_fn, args.treesumm_fn, args.mutlist_fn)

main()
