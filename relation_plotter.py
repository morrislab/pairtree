from common import Models
import colorlover as cl
import numpy as np
import common

def find_ml_relations(model_probs_tensor):
  M = len(model_probs_tensor)
  relations = np.argmax(model_probs_tensor, axis=2)
  assert relations.shape == (M, M)
  return relations

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
  for midx, model in enumerate(Models._all):
    colour = make_colour_from_category(midx)
    print('<tr><td style="background-color: %s">%s</td></tr>' % (colour, model), file=outf)
  print('</tbody></table>', file=outf)

def plot_ml_relations(model_probs_tensor, outf):
  ml_relations = find_ml_relations(model_probs_tensor)
  ml_relations, ssmidxs = common.reorder_square_matrix(ml_relations)
  colours = make_colour_matrix(ml_relations, make_colour_from_category)
  write_table('ML relations', ml_relations, ['C%s' % I for I in ssmidxs], colours, outf)
  write_legend(outf)
