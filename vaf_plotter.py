import colorlover as cl
import csv
import json
import numpy as np
from collections import defaultdict

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

def print_vafs(clustered_vars, supervars, garbage_variants, phi, sampnames, outf):
  nclusters = len(clustered_vars)
  cluster_colours = assign_colours(nclusters)

  print('<style>#vafmatrix td, #vafmatrix { padding: 5px; margin: 0; border-collapse: collapse; } #vafmatrix th { transform: rotate(45deg); font-weight: normal !important; } #vafmatrix span { visibility: hidden; } #vafmatrix td:hover > span { visibility: visible; }</style>', file=outf)
  print('<br><br><br><table id="vafmatrix" class="matrix"><thead>', file=outf)
  header = ['Gene', 'ID', 'Chrom', 'Locus', 'Cluster']
  header += munge_samp_names(sampnames)
  print(''.join(['<th>%s</th>' % H for H in header]), file=outf)
  print('</thead><tbody>', file=outf)

  parted_garbage_vars = partition_garbage_variants(supervars, garbage_variants)
  phi_pseudovars = make_phi_pseudovars(phi)

  for cidx, cluster in enumerate(clustered_vars):
    if len(cluster) == 0:
      continue
    supervar = supervars['C%s' % cidx]
    garbage = parted_garbage_vars[cidx] if cidx in parted_garbage_vars else []
    phi_pseudovar = phi_pseudovars[cidx]

    for V in [phi_pseudovar, supervar] + cluster + garbage:
      td = ['<td>%s</td>' % (V[K] if V[K] is not None else '&mdash;') for K in ('gene', 'id', 'chrom', 'pos', 'cluster')]
      td += ['<td style="background-color: %s"><span>%s</span></td>' % (make_colour(v), make_vaf_label(v)) for v in V['vaf']]
      print('<tr style="background-color: %s">%s</tr>' % (
        cluster_colours[V['cluster']],
        ''.join(td)
      ), file=outf)
  print('</tbody></table>', file=outf)

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

def plot_vaf_matrix(clusters, variants, supervars, garbage_variants, phi, paramsfn, spreadsheetfn, outf):
  spreadsheet = load_spreadsheet(spreadsheetfn)
  with open(paramsfn) as P:
    params = json.load(P)
  sampnames = params['samples']

  for V in list(variants.values()) + list(garbage_variants.values()):
    V['chrom'], V['pos'] = V['name'].split('_')
    V['pos'] = int(V['pos'])
    V['gene'] = find_gene_name(V['chrom'], V['pos'], spreadsheet)
    V['vaf'] = V['var_reads'] / V['total_reads']

  # Copy variant so we don't modify original dict.
  clustered_vars = [[dict(variants['s%s' % vid]) for vid in C] for C in clusters]
  for cidx, cluster in enumerate(clustered_vars):
    for var in cluster:
      var['cluster'] = cidx
  print_vafs(clustered_vars, supervars, garbage_variants, phi, sampnames, outf)
