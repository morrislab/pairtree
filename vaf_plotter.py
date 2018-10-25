import colorlover as cl
import csv
import json
import numpy as np
from collections import defaultdict
import common
import itertools

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
  if spreadsheetfn is None:
    return None
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
  omega_v = 0.5
  V = [{
    'gene': None,
    'id': 'P%s' % cidx,
    'name': 'P%s' % cidx,
    'chrom': None,
    'pos': None,
    'cluster': cidx,
    'vaf': omega_v * row,
    'omega_v': omega_v,
  } for cidx, row in enumerate(phi)]
  return V

def print_vaftable_header(sampnames, outf):
  print('<style>.vafmatrix td, .vafmatrix { padding: 5px; margin: 0; border-collapse: collapse; } .vafmatrix th { transform: rotate(45deg); font-weight: normal !important; } .vafmatrix span { visibility: hidden; } .vafmatrix td:hover > span { visibility: visible; }</style>', file=outf)
  print('''<script type="text/javascript">$(document).ready(function() {''', file=outf)
  print('new VafMatrix("#vafmatrix_toggles");', file=outf)
  print('});</script>', file=outf);
  print('<div id="vafmatrix_toggles" class="btn-group" data-toggle="buttons">', file=outf)
  print('<label class="btn btn-primary active toggle_phi"><input type="checkbox" autocomplete="off" checked> &phi;</label>', file=outf)
  print('<label class="btn btn-primary active toggle_cluster_means"><input type="checkbox" autocomplete="off" checked> Cluster means</label>', file=outf)
  print('<label class="btn btn-primary active toggle_cluster_members"><input type="checkbox" autocomplete="off" checked> Cluster members</label>', file=outf)
  print('<label class="btn btn-primary active toggle_garbage"><input type="checkbox" autocomplete="off" checked> Garbage</label>', file=outf)
  print('</div>', file=outf)
  print('<br><br><br>', file=outf)
  print('<table class="vafmatrix matrix"><thead><tr>', file=outf)
  header = ['Gene', 'ID', 'Chrom', 'Locus', 'Cluster']
  header += munge_samp_names(sampnames)
  print(''.join(['<th>%s</th>' % H for H in header]), file=outf)
  print('</tr></thead><tbody>', file=outf)

def print_vaftable_row(V, bgcolour, should_correct_vaf, outf):
  if should_correct_vaf:
    vaf = V['vaf'] / V['omega_v']
  else:
    vaf = V['vaf']
  td = ['<td class="%s">%s</td>' % (K, V[K] if K in V and V[K] is not None else '&mdash;') for K in ('gene', 'id', 'chrom', 'pos', 'cluster')]
  td += ['<td style="background-color: %s"><span>%s</span></td>' % (make_colour(v), make_vaf_label(v)) for v in vaf]
  print('<tr style="background-color: %s">%s</tr>' % (
    bgcolour,
    ''.join(td)
  ), file=outf)

def print_vaftable_footer(outf):
  print('</tbody></table>', file=outf)

def print_vafs(clustered_vars, supervars, garbage_variants, phi, sampnames, should_correct_vaf, outf):
  nclusters = len(clustered_vars)
  cluster_colours = assign_colours(nclusters)
  parted_garbage_vars = partition_garbage_variants(supervars, garbage_variants)
  phi_pseudovars = make_phi_pseudovars(phi)

  print_vaftable_header(sampnames, outf)

  for cidx, cluster in enumerate(clustered_vars):
    if len(cluster) == 0:
      continue
    supervar = supervars['C%s' % (cidx - 1)]
    garbage = parted_garbage_vars[cidx] if cidx in parted_garbage_vars else []
    phi_pseudovar = phi_pseudovars[cidx]

    cluster_rows = [phi_pseudovar, supervar]
    cluster_rows += cluster
    cluster_rows += garbage
    for V in cluster_rows:
      print_vaftable_row(V, cluster_colours[V['cluster']], should_correct_vaf, outf)

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

def augment_variant(V, spreadsheet):
  if spreadsheet is None:
    return
  V['gene'] = find_gene_name(V['chrom'], V['pos'], spreadsheet)

def print_distances(sampid, supervars, phi):
  def _compute_dist(V1, V2):
    dist = np.sum(np.abs(V1['vaf'] - V2['vaf']))
    normdist = dist / len(V1['vaf'])
    return (dist, normdist)

  phi_pseudovars = {P['id']: P for P in make_phi_pseudovars(phi)}
  sv_ids = sorted(supervars.keys(), key = lambda S: int(S[1:]))
  phi_ids = sorted(phi_pseudovars.keys(), key = lambda P: int(P[1:]))

  # Add 1 to supervariant IDs to match cluster numbering in the results --
  # supervariant C_i corresponds to cluter K_(i + 1).  Actually, maybe this is
  # the wrong thing to do, since I'm already removing empty clusters, meaning
  # that my numbering doesn't match Steph's. Meh.
  for V1, V2 in itertools.combinations(sv_ids, 2):
    dist, normdist = _compute_dist(supervars[V1], supervars[V2])
    print('cluster_dist', int(V1[1:]) + 1, int(V2[1:]) + 1, dist, normdist, sep=',')
  for sv in sv_ids:
    phivid = 'P%s' % (int(sv[1:]) + 1)
    dist, normdist = _compute_dist(supervars[sv], phi_pseudovars[phivid])
    print('phi_dist', int(phivid[1:]), int(sv[1:]) + 1, dist, normdist, sep=',')

def plot_vaf_matrix(sampid, clusters, variants, supervars, garbage_variants, phi, sampnames, spreadsheetfn, should_correct_vaf, outf):
  print('<h2>VAFs (%s)</h2>' % ('corrected' if should_correct_vaf else 'uncorrected'), file=outf)
  spreadsheet = load_spreadsheet(spreadsheetfn)

  # Copy variant so we don't modify original dict.
  variants         = {vid: dict(variants[vid])         for vid in variants.keys()}
  garbage_variants = {vid: dict(garbage_variants[vid]) for vid in garbage_variants.keys()}
  for V in list(variants.values()) + list(garbage_variants.values()):
    augment_variant(V, spreadsheet)

  clustered_vars = [[variants['s%s' % vid] for vid in C] for C in clusters]
  for cidx, cluster in enumerate(clustered_vars):
    for var in cluster:
      var['cluster'] = cidx
  print_vafs(clustered_vars, supervars, garbage_variants, phi, sampnames, should_correct_vaf, outf)
  print_distances(sampid, supervars, phi)
