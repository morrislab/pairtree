import colorlover as cl
import csv
import json
import numpy as np
from collections import defaultdict
import common
import itertools
import random
import string

def make_random_id(K=8):
  return ''.join(random.choices(string.ascii_letters, k=K))

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

def partition_garbage_variants(supervars, garbage_variants):
  for gvar in garbage_variants.values():
    gvar['cluster'] = None
  if supervars is None:
    return {0: list(garbage_variants.values())}

  supervafs = [C['vaf'] for C in supervars]
  parted = defaultdict(list)

  for gvar in sorted(garbage_variants.values(), key = lambda V: int(V['id'][1:])):
    closest_idx = find_closest(gvar['vaf'], supervafs)
    cidx = supervars[closest_idx]['cluster']
    parted[cidx].append(gvar)
  return parted

def make_phi_pseudovars(phi):
  M, S = phi.shape
  omega_v = 0.5 * np.ones(S)
  V = [{
    'gene': None,
    'id': 'P%s' % cidx,
    'name': 'P%s' % cidx,
    'chrom': None,
    'pos': None,
    'cluster': cidx,
    'vaf': omega_v * row,
    'omega_v': omega_v,
  } for cidx, row in enumerate(phi[1:])]
  return V

def print_vaftable_header(sampnames, outf):
  container_id = 'vafmatrix_' + make_random_id()
  container = f'#{container_id}'

  print('<style>{container} .matrix td, {container} .matrix {{ padding: 5px; margin: 0; border-collapse: collapse; }} {container} .matrix th {{ transform: rotate(45deg); font-weight: normal !important; }} {container} .matrix span {{ visibility: hidden; }} {container} .matrix td:hover > span {{ visibility: visible; }}</style>'.format(container=container), file=outf)
  print('''<script type="text/javascript">$(document).ready(function() {''', file=outf)
  print(f'new VafMatrix("{container}");', file=outf)
  print('});</script>', file=outf)

  print(f'<div id="{container_id}">', file=outf)
  print('<p><input type="text" class="filter" placeholder="s0,s1,..."></p>', file=outf)
  print('<div class="vafmatrix_toggles btn-group" data-toggle="buttons">', file=outf)
  print('<label class="btn btn-primary active toggle_phi"><input type="checkbox" autocomplete="off" checked> &phi;</label>', file=outf)
  print('<label class="btn btn-primary active toggle_supervars"><input type="checkbox" autocomplete="off" checked> Supervariants</label>', file=outf)
  print('<label class="btn btn-primary active toggle_cluster_members"><input type="checkbox" autocomplete="off" checked> Cluster members</label>', file=outf)
  print('<label class="btn btn-primary active toggle_garbage"><input type="checkbox" autocomplete="off" checked> Garbage</label>', file=outf)
  print('</div>', file=outf)
  print('<br><br><br>', file=outf)
  print('<table class="matrix"><thead><tr>', file=outf)
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
  print('</div></tbody></table>', file=outf)

def print_vafs(clustered_vars, supervars, garbage_variants, phi, sampnames, should_correct_vaf, outf):
  nclusters = len(clustered_vars)
  cluster_colours = assign_colours(nclusters)
  parted_garbage_vars = partition_garbage_variants(supervars, garbage_variants)
  if phi is not None:
    phi_pseudovars = make_phi_pseudovars(phi)

  print_vaftable_header(sampnames, outf)

  for cidx, cluster in enumerate(clustered_vars):
    assert len(cluster) > 0
    garbage = parted_garbage_vars[cidx] if cidx in parted_garbage_vars else []

    cluster_rows = []
    if phi is not None:
      cluster_rows.append(phi_pseudovars[cidx])
    if supervars is not None:
      cluster_rows.append(supervars[cidx])
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

def plot_vaf_matrix(clusters, variants, supervars, garbage_vids, phi, sampnames, should_correct_vaf, outf):
  assert len(supervars) == len(clusters)
  print('<h2>VAFs (%s)</h2>' % ('corrected' if should_correct_vaf else 'uncorrected'), file=outf)

  # Copy variant so we don't modify original dict.
  variants = {vid: dict(variants[vid]) for vid in variants.keys()}
  garbage_variants = {vid: dict(variants[vid]) for vid in garbage_vids}

  clustered_vars = [[variants[vid] for vid in C] for C in clusters]
  for cidx, cluster in enumerate(clustered_vars):
    if supervars is not None:
      supervars[cidx]['cluster'] = cidx
    for var in cluster:
      var['cluster'] = cidx

  print_vafs(clustered_vars, supervars, garbage_variants, phi, sampnames, should_correct_vaf, outf)
  #print_distances(sampid, supervars, phi)
