import csv
import argparse
import json
import numpy as np
import sklearn.cluster
import colorlover as cl
from common import parse_ssms
from vaf_plotter import plot_vaf_matrix
np.set_printoptions(threshold=np.nan)

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
  scalenum = 4
  assert 0 <= cat < scalenum
  scale = cl.scales[str(scalenum)]['qual']['Dark2']
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
  for model in model_probs['models']:
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

def write_legend(models, outf):
  print('<br><table class="table"><tbody>', file=outf)
  for midx, model in enumerate(models):
    colour = make_colour_from_category(midx)
    print('<tr><td style="background-color: %s">%s</td></tr>' % (colour, model), file=outf)
  print('</tbody></table>', file=outf)

def write_cluster_map(cmap, cidxs, outf):
  assert len(cmap) == len(cidxs)
  print('<br><table class="table"><thead><tr><th>Cluster</th><th>SSMs</th></tr></thead><tbody>', file=outf)
  for cidx, ssmidxs in zip(cidxs, cmap):
    ssmidxs = ', '.join(['s%s' % I for I in sorted(ssmidxs)])
    print('<tr><td style="font-weight: bold">C%s</td><td>%s</td></tr>' % (cidx, ssmidxs), file=outf)
  print('</tbody></table>', file=outf)

def calc_mle(model_probs):
  mats = np.array([create_matrix(M, model_probs['model_probs'][M], model_probs['var_names']) for M in model_probs['models']])
  mle = np.argmax(mats, axis=0)
  return mle

def plot_mle(model_probs, should_cluster, outf):
  mle = calc_mle(model_probs)

  if should_cluster:
    mle, ssmidxs = cluster_rows(mle)
  else:
    ssmidxs = list(range(len(mle)))

  colours = make_colour_matrix(mle, make_colour_from_category)
  write_legend(model_probs['models'], outf)
  write_table('mle', mle, ['s%s' % I for I in ssmidxs], colours, outf)

def collapse_identical(mat):
  '''Collapse identical rows & columns.'''
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

def remove_small_clusters(mat, clusters, cidxs, threshold=1):
  assert len(clusters) == len(mat) == len(cidxs)

  N = len(mat)
  to_remove = set([idx for idx, C in enumerate(clusters) if len(C) <= threshold])
  to_keep = [idx for idx in range(N) if idx not in to_remove]
  assert len(to_remove) + len(to_keep) == N

  filtered_mat = mat[to_keep,:][:,to_keep]
  filtered_clusters = [C for idx, C in enumerate(clusters) if idx not in to_remove]
  filtered_cidxs = [C for idx, C in enumerate(cidxs) if idx not in to_remove]
  assert len(filtered_clusters) == len(to_keep) == len(filtered_cidxs)

  return (filtered_mat, filtered_clusters, filtered_cidxs)

def cluster_mle(model_probs):
  mle = calc_mle(model_probs)
  sidxs_toposort = toposort(mle, model_probs['models'])
  # Sort rows by toposorted indexes.
  mle_toposort = np.array([_reorder_row(mle[idx], sidxs_toposort) for idx in sidxs_toposort])

  collapsed, idxmap = collapse_identical(mle_toposort)
  row_to_sidx_map = dict(enumerate(sidxs_toposort))
  clusters = [[row_to_sidx_map[rowidx] for rowidx in cluster] for cluster in idxmap]
  cidxs = range(len(clusters))

  return (collapsed, clusters, cidxs)

def plot_mle_toposort(model_probs, outf, remove_small=False):
  mle, clusters, cidxs = cluster_mle(model_probs)
  if remove_small:
    mle, clusters, cidxs = remove_small_clusters(mle, clusters, cidxs)

  colours = make_colour_matrix(mle, make_colour_from_category)
  labels = ['C%s' % I for I in cidxs]
  write_cluster_map(clusters, cidxs, outf)
  suffix = remove_small and 'small_excluded' or 'small_included'
  write_table('mle_toposort_%s' % suffix, mle, labels, colours, outf)

def extract_B_A_rels(mle, models):
  ssmidxs = list(range(len(mle)))
  A_B_rels, B_A_rels = set(), set()
  B_A_adjlist = {}

  for sidx1 in ssmidxs:
    B_A_adjlist[sidx1] = set()

    for sidx2 in ssmidxs:
      rel = mle[sidx1,sidx2]
      if models[rel] == 'A_B':
        A_B_rels.add((sidx1, sidx2))
      elif models[rel] == 'B_A':
        B_A_rels.add((sidx1, sidx2))
        B_A_adjlist[sidx1].add(sidx2)

  assert len(A_B_rels) == len(B_A_rels)
  reversed_B_A_rels = set([(S2, S1) for (S1, S2) in B_A_rels])
  assert reversed_B_A_rels == A_B_rels

  return B_A_adjlist

def toposort(mle, models):
  # Algorithm taken from https://en.wikipedia.org/w/index.php?title=Topological_sorting&oldid=779516160#Kahn.27s_algorithm.
  B_A_rels = extract_B_A_rels(mle, models)
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

def plot(sampid, model_probs, should_cluster, ssmfn, paramsfn, spreadsheetfn, outfn):
  with open(outfn, 'w') as outf:
    write_header(sampid, should_cluster and 'clustered' or 'unclustered', outf)
    plot_individual(model_probs, should_cluster, outf)
    plot_mle(model_probs, should_cluster, outf)
    plot_mle_toposort(model_probs, outf)
    plot_mle_toposort(model_probs, outf, remove_small=True)

    _, clusters, cidxs = cluster_mle(model_probs)
    plot_vaf_matrix(clusters, cidxs, ssmfn, paramsfn, spreadsheetfn, outf)

def load_model_probs(model_probs_fn):
  with open(model_probs_fn) as F:
    return json.load(F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--cluster', dest='should_cluster', action='store_true')
  parser.add_argument('sampid')
  parser.add_argument('model_probs_fn')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('spreadsheet_fn')
  parser.add_argument('out_fn')
  args = parser.parse_args()

  model_probs = load_model_probs(args.model_probs_fn)
  plot(args.sampid, model_probs, args.should_cluster, args.ssm_fn, args.params_fn, args.spreadsheet_fn, args.out_fn)

main()
