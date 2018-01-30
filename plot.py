import csv
import argparse
import json
import numpy as np
import colorlover as cl
import common
from common import Models
import vaf_plotter
from collections import defaultdict
import tree_sampler
import tree_builder
import json_writer
import phi_fitter
import handbuilt
import pairwise
import eta_plotter
import vaf_correcter
import cluster

np.set_printoptions(threshold=np.nan)
np.random.seed(1)

def plot_cluster_mle_relations(supervars, svid2svidx, svidx2svid, outf):
  posterior, evidence = pairwise.calc_posterior(supervars)
  results = pairwise.generate_results(posterior, evidence, supervars)
  model_probs_tensor = create_model_prob_tensor(results, svid2svidx)
  supervar_relations = calc_relations(model_probs_tensor)

  # Sort rows by toposorted indexes.
  sorted_idxs = toposort(supervar_relations)
  sortmap = dict(enumerate(sorted_idxs))
  supervar_relations = np.array([_reorder_row(supervar_relations[idx], sorted_idxs) for idx in sorted_idxs])
  svid2svidx = {K: sortmap[V] for K, V in svid2svidx.items()}
  svidx2svid = {sortmap[K]: V for K, V in svidx2svid.items()}

  plot_relations(supervar_relations, False, svidx2svid, outf)
  plot_individual(results, False, svid2svidx, svidx2svid, outf)

  return supervar_relations

def create_matrix(model, model_probs, variants, vid2vidx):
  N = len(variants)
  mat = np.zeros((N, N))
  for vids, P in model_probs.items():
    assert 0 <= P <= 1
    vidx1, vidx2 = [vid2vidx[V] for V in vids.split(',')]
    mat[(vidx1, vidx2)] = P
  if model in ('cocluster', 'diff_branches'):
    # These should be symmetric.
    assert np.allclose(mat, mat.T)

  return mat

def _reorder_row(R, idxs):
  return np.array([R[idx] for idx in idxs])

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
  print('<script src="https://d3js.org/d3.v4.min.js"></script>', file=outf)
  print('<script src="https://d3js.org/d3-scale-chromatic.v1.min.js"></script>', file=outf)
  print('<script type="text/javascript" src="highlight_table_labels.js"></script>', file=outf)
  print('<script type="text/javascript" src="tree_plotter.js"></script>', file=outf)
  print('<link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet">', file=outf)
  print('<h1>%s (%s)</h1>' % (sampid, extra), file=outf)
  print('<link href="tree.css" rel="stylesheet">', file=outf)
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

def plot_individual(model_probs, should_cluster, vid2vidx, vidx2vid, outf):
  for model in Models._all:
    mat = create_matrix(model, model_probs['model_probs'][model], model_probs['variants'], vid2vidx)
    if should_cluster:
      mat, ssmidxs = common.reorder_square_matrix(mat)
    else:
      ssmidxs = list(range(len(mat)))
    if model in ('cocluster', 'diff_branches'):
      # These should be symmetric.
      assert np.allclose(mat, mat.T)
    colours = make_colour_matrix(mat, make_colour_from_intensity)
    write_table(model, mat, [vidx2vid[I] for I in ssmidxs], colours, outf)

def write_legend(outf):
  print('<br><table class="table"><tbody>', file=outf)
  for midx, model in enumerate(Models._all):
    colour = make_colour_from_category(midx)
    print('<tr><td style="background-color: %s">%s</td></tr>' % (colour, model), file=outf)
  print('</tbody></table>', file=outf)

def create_model_prob_tensor(model_probs, vid2vidx):
  M = len(model_probs['variants'])
  num_models = len(Models._all)
  tensor = np.zeros((M, M, num_models))
  for midx, mdl in enumerate(Models._all):
    tensor[:,:,midx] = create_matrix(mdl, model_probs['model_probs'][mdl], model_probs['variants'], vid2vidx)
  return tensor

def calc_relations(model_probs):
  M = len(model_probs)
  relations = np.argmax(model_probs, axis=2)
  assert relations.shape == (M, M)
  return relations

def plot_relations(relations, should_cluster, vidx2vid, outf):
  if should_cluster:
    relations, ssmidxs = common.reorder_square_matrix(relations)
  else:
    ssmidxs = list(range(len(relations)))

  colours = make_colour_matrix(relations, make_colour_from_category)
  write_table('relations', relations, [vidx2vid[I] for I in ssmidxs], colours, outf)

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

def fit_phis(adjm, variants, clusters, tidxs=None):
  ntrees = len(adjm)
  nsamples = len(list(variants.values())[0]['total_reads'])

  eta = np.ones((ntrees, nsamples, len(clusters)))
  phi = np.ones((ntrees, nsamples, len(clusters)))

  if tidxs is None:
    tidxs = range(ntrees)
  for tidx in tidxs:
    phi[tidx,:,:], eta[tidx,:,:] = phi_fitter.fit_phis(adjm[tidx], clusters, variants)

  return (phi, eta)

def remove_garbage(garbage_ids, model_probs, variants, clusters):
  garbage_variants = {}

  for varid in sorted(variants.keys(), key = lambda S: int(S[1:])):
    if int(varid[1:]) in garbage_ids:
      garbage_variants[varid] = variants[varid]
      del variants[varid]
      del model_probs['variants'][varid]

  for K in ('model_probs', 'model_evidence'):
    for model in model_probs[K].keys():
      for pair in list(model_probs[K][model].keys()):
        vidx1, vidx2 = [int(V[1:]) for V in pair.split(',')]
        if vidx1 in garbage_ids or vidx2 in garbage_ids:
          del model_probs[K][model][pair]

  return garbage_variants

def load_sampnames(paramsfn):
  with open(paramsfn) as P:
    params = json.load(P)
  sampnames = params['samples']
  return sampnames

def cluster_samples(variants, sampnames):
  vaf = make_vaf_matrix(variants)
  _, idxs = cluster_cols(vaf)

  for var in variants.values():
    for K in ('ref_reads', 'var_reads', 'total_reads'):
      var[K] = var[K][idxs]
  sampnames = [sampnames[idx] for idx in idxs]

  return (variants, sampnames)

def write_trees(sampid, colourings, outf):
  print('''<div id="trees"></div>''', file=outf)
  print('''<script type="text/javascript">$(document).ready(function() {''', file=outf)
  for tname, tidx in (('handbuilt', 0), ('sampled', 1000)):
    for colouring in colourings:
      print('''(new TreePlotter()).plot('%s.summ.json', %s, '%s', '%s', '%s', '#trees');''' % (
        sampid,
        tidx,
        tname,
        colouring['left'],
        colouring['right'],
      ), file=outf)
  print('});</script>', file=outf)
  #for tree_type in ('handbuilt', 'sampled'):
    #print('<h2>%s</h2><div id="tree-%s"></div>' % (tree_type, tree_type), file=outf)

def write_phi_matrix(sampid, outf):
  print('''<script type="text/javascript">$(document).ready(function() {
  (new PhiMatrix()).plot('%s', '%s.phi.json', '#phi_matrix');
  });</script>''' % (sampid, sampid), file=outf)
  print('<div id="phi_matrix" style="margin: 30px"></div>', file=outf)

def plot(sampid, model_probs, output_type, tree_type, ssmfn, paramsfn, spreadsheetfn, handbuiltfn, outfn, treesummfn, mutlistfn, phifn):
  sampnames = load_sampnames(paramsfn)
  variants = common.parse_ssms(sampid, ssmfn)
  if tree_type == 'handbuilt.patient':
    variants, sampnames = common.extract_patient_samples(variants, sampnames)
  #cluster.cluster_patient_vars(sampid, variants, sampnames)

  garbage_ids = handbuilt.load_garbage(handbuiltfn, tree_type)
  clusters, handbuilt_adjm, node_colourings = handbuilt.load_clusters_and_tree(handbuiltfn, variants, tree_type, sampnames)
  samporders = handbuilt.load_samporders(handbuiltfn, tree_type)

  garbage_variants = remove_garbage(garbage_ids, model_probs, variants, clusters)
  vidxs = sorted(model_probs['variants'].keys(), key = lambda V: int(V[1:]))
  vidx2vid = dict(enumerate(vidxs))
  vid2vidx = {V: K for K, V in vidx2vid.items()}

  supervars, svid2svidx, svidx2svid = common.make_cluster_supervars(clusters, variants)

  should_cluster = not (output_type == 'unclustered')
  model_probs_tensor = create_model_prob_tensor(model_probs, vid2vidx)
  ssm_relations = calc_relations(model_probs_tensor)
  #supervar_relations = plot_cluster_mle_relations(supervars, svid2svidx, svidx2svid, outf)
  sampled_adjm, sampled_llh = tree_sampler.sample_trees(model_probs_tensor, clusters, vid2vidx, 1000)

  with open(outfn, 'w') as outf:
    write_header(sampid, '%s, %s' % (output_type, tree_type), outf)
    #try:
    #  mle_adjm = tree_builder.make_adj(supervar_relations)
    #  mle_adjm = add_normal_root(mle_adjm)
    #except tree_builder.CannotBuildTreeException:
    #  mle_adjm = None
    mle_adjm = None

    for adjm in (mle_adjm, handbuilt_adjm):
      if adjm is None:
        continue
      mutrel = tree_sampler.make_mutrel_tensor_from_cluster_adj(adjm, clusters, vid2vidx)
      llh = tree_sampler.calc_llh(model_probs_tensor, mutrel)
      sampled_adjm.insert(0, adjm)
      sampled_llh.insert(0, llh)

    phi, eta = fit_phis(sampled_adjm, variants, clusters, tidxs=(0, -1))
    json_writer.write_json(sampid, sampnames, variants, clusters, sampled_adjm, sampled_llh, phi, treesummfn, mutlistfn)

    eta_plotter.plot_eta(eta[0].T, sampnames, samporders, outf)
    eta_plotter.write_phi_json(phi[0].T, sampnames, samporders, phifn)
    write_trees(sampid, node_colourings, outf)
    write_phi_matrix(sampid, outf)
    for correct_vafs in (True, False):
      if correct_vafs is True and not vaf_correcter.has_corrections(sampid):
        continue
      vaf_plotter.plot_vaf_matrix(sampid, clusters, variants, supervars, garbage_variants, phi[0].T, sampnames, spreadsheetfn, correct_vafs, outf)
    vaf_plotter.plot_unclustered_vafs(sampid, variants, None,             sampnames, spreadsheetfn, outf, patient_samples_only=False)
    #vaf_plotter.plot_unclustered_vafs(sampid, variants, garbage_variants, sampnames, spreadsheetfn, outf, patient_samples_only=True)
    #plot_individual(model_probs, should_cluster, vid2vidx, vidx2vid, outf)
    #plot_relations(ssm_relations, should_cluster, vidx2vid, outf)
    write_legend(outf)

def load_model_probs(model_probs_fn):
  with open(model_probs_fn) as F:
    return json.load(F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--output-type', dest='output_type', choices=('clustered', 'unclustered', 'condensed'), default='clustered')
  parser.add_argument('--tree-type', dest='tree_type', required=True)
  parser.add_argument('sampid')
  parser.add_argument('model_probs_fn')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('spreadsheet_fn')
  parser.add_argument('handbuilt_fn')
  parser.add_argument('out_fn')
  parser.add_argument('treesumm_fn')
  parser.add_argument('mutlist_fn')
  parser.add_argument('phi_fn')
  args = parser.parse_args()

  model_probs = load_model_probs(args.model_probs_fn)
  plot(args.sampid, model_probs, args.output_type, args.tree_type, args.ssm_fn, args.params_fn, args.spreadsheet_fn, args.handbuilt_fn, args.out_fn, args.treesumm_fn, args.mutlist_fn, args.phi_fn)

if __name__ == '__main__':
  main()
