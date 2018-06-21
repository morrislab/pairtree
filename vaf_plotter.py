import colorlover as cl
import csv
import json
import numpy as np
from collections import defaultdict
import vaf_correcter
import common

def find_gene_name(chrom, pos, spreadsheet_rows):
  for row in spreadsheet_rows:
    assert row['CHR'].startswith('chr')
    C = row['CHR'][3:].upper()
    P = int(row['WU_HG19_POS'])

    if chrom == C and pos == P:
      return row['GENENAME']
  return 'unknown'
  raise Exception('Could not find gene name for %s_%s' % (chrom, pos))

def load_spreadsheet(spreadsheetfn):
  with open(spreadsheetfn) as S:
    reader = csv.DictReader(S)
    return list(reader)

def munge_samp_names(sampnames):
  return [S.replace('Diagnosis Xeno ', 'DX').replace('Relapse Xeno ', 'RX') for S in sampnames]

def euclid_dist(A, B):
  return np.sqrt(np.sum((B - A)**2))

def find_closest(vec, mat):
  # Find matrix in tensor that's closest to `target` in Euclidean distance.
  min_dist = float('inf')
  best_idx = None

  for idx, candidate in enumerate(mat):
    if candidate is None:
      continue
    dist = euclid_dist(candidate, vec)
    if dist < min_dist:
      min_dist = dist
      best_idx = idx

  assert best_idx is not None
  return best_idx

def partition_garbage_variants(cluster_supervars, garbage_variants):
  supervars = list(cluster_supervars.values())
  supervafs = [C['vaf'] for C in supervars]
  parted = defaultdict(list)

  for gvar in sorted(garbage_variants.values(), key = lambda V: int(V['id'][1:])):
    gvar['cluster'] = None
    closest_idx = find_closest(gvar['vaf'], supervafs)
    cidx = supervars[closest_idx]['cluster']
    parted[cidx].append(gvar)
  return parted

def make_phi_pseudovars(phi):
  V = [{
    'gene': None,
    'id': 'P%s' % cidx,
    'name': 'P%s' % cidx,
    'chrom': None,
    'pos': None,
    'cluster': cidx,
    'vaf': 0.5 * row,
  } for cidx, row in enumerate(phi)]
  return V

def print_vaftable_header(sampnames, outf):
  print('<style>.vafmatrix td, .vafmatrix { padding: 5px; margin: 0; border-collapse: collapse; } .vafmatrix th { transform: rotate(45deg); font-weight: normal !important; } .vafmatrix span { visibility: hidden; } .vafmatrix td:hover > span { visibility: visible; }</style>', file=outf)
  print('<br><br><br><table class="vafmatrix matrix"><thead>', file=outf)
  header = ['Gene', 'ID', 'Chrom', 'Locus', 'Cluster']
  header += munge_samp_names(sampnames)
  print(''.join(['<th>%s</th>' % H for H in header]), file=outf)
  print('</thead><tbody>', file=outf)

def print_vaftable_row(V, bgcolour, outf):
  td = ['<td>%s</td>' % (V[K] if V[K] is not None else '&mdash;') for K in ('gene', 'id', 'chrom', 'pos', 'cluster')]
  td += ['<td style="background-color: %s"><span>%s</span></td>' % (make_colour(v), make_vaf_label(v)) for v in V['vaf']]
  print('<tr style="background-color: %s">%s</tr>' % (
    bgcolour,
    ''.join(td)
  ), file=outf)

def print_vaftable_footer(outf):
  print('</tbody></table>', file=outf)

def print_vafs(clustered_vars, supervars, garbage_variants, phi, sampnames, outf):
  nclusters = len(clustered_vars)
  cluster_colours = assign_colours(nclusters)
  parted_garbage_vars = partition_garbage_variants(supervars, garbage_variants)
  phi_pseudovars = make_phi_pseudovars(phi)

  print_vaftable_header(sampnames, outf)

  for cidx, cluster in enumerate(clustered_vars):
    if len(cluster) == 0:
      continue
    supervar = supervars['C%s' % cidx]
    garbage = parted_garbage_vars[cidx] if cidx in parted_garbage_vars else []
    phi_pseudovar = phi_pseudovars[cidx]
    for V in [phi_pseudovar, supervar] + cluster + garbage:
      print_vaftable_row(V, cluster_colours[V['cluster']], outf)

  print_vaftable_footer(outf)

def reorder_variants(variants, sampnames):
  vids = list(variants.keys())
  vaf = np.array([variants[vid]['vaf'] for vid in vids])

  vaf, vidxs = common.reorder_rows(vaf)
  vids = [vids[I] for I in vidxs]
  vaf, sampidxs = common.reorder_cols(vaf)
  sampnames = [sampnames[I] for I in sampidxs]

  ordered_variants = []
  for vid in vids:
    variant = dict(variants[vid])
    for K in ('total_reads', 'ref_reads', 'var_reads', 'vaf'):
      variant[K] = [variant[K][I] for I in sampidxs]
    ordered_variants.append(variant)

  return (ordered_variants, sampnames)

def print_unclustered_vafs(variants, sampnames, outf, patient_samples_only=False):
  if patient_samples_only:
    variants, sampnames = common.extract_patient_samples(variants, sampnames)
  ordered_variants, sampnames = reorder_variants(variants, sampnames)
  print_vaftable_header(sampnames, outf)
  for variant in ordered_variants:
    print_vaftable_row(variant, '#fff', outf)
  print_vaftable_footer(outf)

def make_colour(vaf):
  val = int(255*(1 - float(vaf)))
  return 'rgb(255, %s, %s)' % (2*(val,))

def make_vaf_label(vaf):
  return '%.2f' % float(vaf)

def assign_colours(num_colours):
  colours = {None: '#fff'}
  for cidx in range(num_colours):
    colours[cidx] = get_next_colour()
  return colours

def get_next_colour():
  scale = cl.scales['11']['qual']['Set3']
  L = len(scale)
  idx = get_next_colour._last_idx + 1
  if idx >= L:
    idx = 0
  get_next_colour._last_idx = idx
  return scale[idx]
get_next_colour._last_idx = -1

def augment_variant(V, spreadsheet, correct_vaf):
  if spreadsheet is None:
    return
  if correct_vaf:
    V['vaf'] *= V['vaf_correction']
  V['gene'] = find_gene_name(V['chrom'], V['pos'], spreadsheet)

def plot_unclustered_vafs(sampid, variants, garbage_variants, sampnames, spreadsheetfn, outf, patient_samples_only=False):
  print('<h2>Unclustered VAFs (corrected, %s samples)</h2>' % ('patient' if patient_samples_only else 'all') , file=outf)
  print('<h3>Corrected variants: %s</h3>' % ', '.join(vaf_correcter.corrected_vars(sampid)), file=outf)
  if spreadsheetfn is not None:
    spreadsheet = load_spreadsheet(spreadsheetfn)
  else:
    spreadsheet = None

  # Copy variant so we don't modify original dict.
  variants = {vid: dict(variants[vid]) for vid in variants.keys()}
  if garbage_variants is not None:
    variants.update({vid: dict(garbage_variants[vid]) for vid in garbage_variants.keys()})
  for V in variants.values():
    augment_variant(V, spreadsheet, True)
    V['cluster'] = None

  print_unclustered_vafs(variants, sampnames, outf, patient_samples_only)

def plot_vaf_matrix(sampid, clusters, variants, supervars, garbage_variants, phi, sampnames, spreadsheetfn, correct_vafs, outf):
  print('<h2>VAFs (%s)</h2>' % ('corrected' if correct_vafs else 'uncorrected',), file=outf)
  if correct_vafs is True:
    print('<h3>Corrected variants: %s</h3>' % ', '.join(vaf_correcter.corrected_vars(sampid)), file=outf)
  spreadsheet = load_spreadsheet(spreadsheetfn)

  # Copy variant so we don't modify original dict.
  variants         = {vid: dict(variants[vid])         for vid in variants.keys()}
  garbage_variants = {vid: dict(garbage_variants[vid]) for vid in garbage_variants.keys()}
  for V in list(variants.values()) + list(garbage_variants.values()):
    augment_variant(V, spreadsheet, correct_vafs)

  clustered_vars = [[variants['s%s' % vid] for vid in C] for C in clusters]
  for cidx, cluster in enumerate(clustered_vars):
    for var in cluster:
      var['cluster'] = cidx
  print_vafs(clustered_vars, supervars, garbage_variants, phi, sampnames, outf)
