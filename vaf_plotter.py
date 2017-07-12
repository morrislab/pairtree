import colorlover as cl
import csv
import json
import numpy as np

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

def make_cluster_supervars(ordered_variants):
  cluster_supervars = []
  for cidx, cluster in enumerate(ordered_variants):
    if len(cluster) == 0:
      cluster_supervars.append(None)
      continue

    cluster_var_reads = np.array([V['var_reads'] for V in cluster])
    cluster_total_reads = np.array([V['total_reads'] for V in cluster])
    cluster_supervars.append({
      'gene': None,
      'id': None,
      'chrom': None,
      'pos': None,
      'cluster': cidx,
      'vaf': np.sum(cluster_var_reads, axis=0) / np.sum(cluster_total_reads, axis=0)
    })

  return cluster_supervars

def partition_garbage_variants(cluster_supervars, garbage_variants):
  supervafs = [C['vaf'] if C is not None else None for C in cluster_supervars]
  parted = [list() for _ in cluster_supervars]
  for gvar in garbage_variants.values():
    gvar['cluster'] = None
    closest = find_closest(gvar['vaf'], supervafs)
    parted[closest].append(gvar)
  return parted

def print_vafs(ordered_variants, garbage_variants, sampnames, outf):
  nclusters = len(ordered_variants)
  cluster_colours = assign_colours(nclusters)

  print('<style>#vafmatrix td, #vafmatrix { padding: 5px; margin: 0; border-collapse: collapse; } #vafmatrix th { transform: rotate(45deg); font-weight: normal !important; } #vafmatrix span { visibility: hidden; } #vafmatrix td:hover > span { visibility: visible; }</style>', file=outf)
  print('<br><br><br><table id="vafmatrix" class="matrix"><thead>', file=outf)
  header = ['Gene', 'ID', 'Chrom', 'Locus', 'Cluster']
  header += munge_samp_names(sampnames)
  print(''.join(['<th>%s</th>' % H for H in header]), file=outf)
  print('</thead><tbody>', file=outf)

  supervars = make_cluster_supervars(ordered_variants)
  parted_garbage_vars = partition_garbage_variants(supervars, garbage_variants)

  for supervar, cluster, cluster_garbage_vars in zip(supervars, ordered_variants, parted_garbage_vars):
    if len(cluster) == 0:
      continue
    for V in [supervar] + cluster + cluster_garbage_vars:
      td = ['<td>%s</td>' % V[K] if V[K] is not None else '&mdash;' for K in ('gene', 'id', 'chrom', 'pos', 'cluster')]
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

def plot_vaf_matrix(clusters, variants, garbage_variants, paramsfn, spreadsheetfn, outf):
  spreadsheet = load_spreadsheet(spreadsheetfn)
  with open(paramsfn) as P:
    params = json.load(P)
  sampnames = params['samples']

  for V in list(variants.values()) + list(garbage_variants.values()):
    V['chrom'], V['pos'] = V['name'].split('_')
    V['pos'] = int(V['pos'])
    V['gene'] = find_gene_name(V['chrom'], V['pos'], spreadsheet)
    V['vaf'] = V['var_reads'] / V['total_reads']

  ordered_variants = []
  for cidx, C in enumerate(clusters):
    ordered_variants.append([])
    for V in C:
      # Copy variant so we don't modify original dict.
      variant = dict(variants['s%s' % V])
      variant['cluster'] = cidx
      ordered_variants[cidx].append(variant)
  print_vafs(ordered_variants, garbage_variants, sampnames, outf)
