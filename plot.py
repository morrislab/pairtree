import csv
import argparse
import json
import numpy as np
import sklearn.cluster
import colorlover as cl
np.set_printoptions(threshold=np.nan)

def load_spreadsheet(spreadsheetfn):
  with open(spreadsheetfn) as S:
    reader = csv.DictReader(S)
    return list(reader)

def find_gene_name(chrom, pos, spreadsheet_rows):
  for row in spreadsheet_rows:
    assert row['CHR'].startswith('chr')
    C = row['CHR'][3:].upper()
    P = int(row['WU_HG19_POS'])

    if chrom == C and pos == P:
      return row['GENENAME']
  return 'unknown'
  raise Exception('Could not find gene name for %s_%s' % (chrom, pos))

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
  print('<h1>%s (%s)</h1>' % (sampid, extra), file=outf)
  print('<style>td, table { padding: 5px; margin: 0; border-collapse: collapse; } span { visibility: hidden; } td:hover > span { visibility: visible; }</style>', file=outf)

def write_table(model, mat, ssmidxs, colours, outf):
  print('<h2>%s</h2>' % model, file=outf)
  print('<table><thead></thead><tbody>', file=outf)

  N = len(mat)
  vids = ['s%s' % idx for idx in ssmidxs]

  entries        = [''] + vids
  visibility     = len(entries)*[True]
  header_colours = len(entries)*['transparent']
  print(make_table_row(entries, visibility, header_colours), file=outf)

  for vid, row, row_colours in zip(vids, mat, colours):
    entries     = [vid]           + ['%.2f' % P for P in row]
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
        assert np.all(mat == mat.T)
    else:
      ssmidxs = list(range(len(mat)))
    colours = make_colour_matrix(mat, make_colour_from_intensity)
    write_table(model, mat, ssmidxs, colours, outf)

def write_legend(models, outf):
  print('<table><tbody>', file=outf)
  for midx, model in enumerate(models):
    colour = make_colour_from_category(midx)
    print('<tr><td style="background-color: %s">%s</td></tr>' % (colour, model), file=outf)
  print('</tbody></table>', file=outf)

def plot_joint(model_probs, outf):
  mats = np.array([create_matrix(M, model_probs['model_probs'][M], model_probs['var_names']) for M in model_probs['models']])
  mle = np.argmax(mats, axis=0)
  ssmidxs = list(range(len(mle)))
  colours = make_colour_matrix(mle, make_colour_from_category)
  write_table('joint', mle, ssmidxs, colours, outf)
  write_legend(model_probs['models'], outf)

def plot(sampid, model_probs, should_cluster, outfn):
  with open(outfn, 'w') as outf:
    write_header(sampid, should_cluster and 'clustered' or 'unclustered', outf)
    plot_individual(model_probs, should_cluster, outf)
    plot_joint(model_probs, outf)

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
  parser.add_argument('out_fn')
  args = parser.parse_args()

  model_probs = load_model_probs(args.model_probs_fn)
  plot(args.sampid, model_probs, args.should_cluster, args.out_fn)

main()
