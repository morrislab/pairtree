from common import Models, NUM_MODELS, ALL_MODELS
import colorlover as cl
import numpy as np
import common

def find_ml_relations(mutrel_posterior):
  M = len(mutrel_posterior)
  relations = np.argmax(mutrel_posterior, axis=2)
  assert relations.shape == (M, M)
  return relations

def make_colour_from_category(cat):
  scalenum = NUM_MODELS
  assert 0 <= cat < scalenum
  scale = cl.scales[str(scalenum)]['qual']['Set1']
  return scale[cat]

def make_colour_from_intensity(intensity):
  val = np.int(np.round(255*(1 - intensity)))
  return 'rgb(255, %s, %s)' % (2*(val,))

def make_colour_matrix(vals, colour_maker):
  N, M = vals.shape
  colours = N*[None]

  for i in range(N):
    colours[i] = M*[None]
    for j in range(M):
      colours[i][j] = colour_maker(vals[i,j])

  return colours

def make_table_row(entries, is_label, colours):
  assert len(entries) == len(is_label)
  entries = [E if L else ('<span>%s</span>' % E) for E, L in zip(entries, is_label)]
  elems = ['th' if L else 'td' for L in is_label]
  return '<tr>' + ''.join(['<%s style="background-color: %s">%s</%s>' % (elem, C, E, elem) for E, C, elem in zip(entries, colours, elems)]) + '</tr>'

def write_table(title, mat, labels, colours, outf):
  print('<h2>%s</h2>' % title, file=outf)
  print('<table class="matrix"><thead>', file=outf)

  N = len(mat)

  entries        = [''] + labels
  is_label       = len(entries)*[True]
  header_colours = len(entries)*['transparent']
  print(make_table_row(entries, is_label, header_colours), file=outf)
  print('</thead><tbody>', file=outf)

  for label, row, row_colours in zip(labels, mat, colours):
    entries     = [label]           + ['%.2f' % P for P in row]
    is_label    = [True]          + N*[False]
    row_colours = ['transparent'] + row_colours
    print(make_table_row(entries, is_label, row_colours), file=outf)

  print('</tbody></table>', file=outf)

def write_legend(outf):
  print('<br><table class="matrix"><tbody>', file=outf)
  for midx, model in enumerate(ALL_MODELS):
    colour = make_colour_from_category(midx)
    print('<tr><td style="background-color: %s">%s</td></tr>' % (colour, model), file=outf)
  print('</tbody></table>', file=outf)

def plot_ml_relations(mutrel_posterior, outf):
  ml_relations = find_ml_relations(mutrel_posterior.rels)
  ml_relations, vidxs = common.reorder_square_matrix(ml_relations)
  colours = make_colour_matrix(ml_relations, make_colour_from_category)
  write_table('ML relations', ml_relations, [mutrel_posterior.vids[vidx] for vidx in vidxs], colours, outf)
  write_legend(outf)

def plot_separate_relations(mutrel_posterior, outf):
  for midx, model in enumerate(ALL_MODELS):
    mat = mutrel_posterior.rels[:,:,midx]
    mat, vidxs = common.reorder_square_matrix(mat)
    if model in ('cocluster', 'diff_branches'):
      # These should be symmetric.
      assert np.allclose(mat, mat.T)
    colours = make_colour_matrix(mat, make_colour_from_intensity)
    write_table(model, mat, [mutrel_posterior.vids[vidx] for vidx in vidxs], colours, outf)
