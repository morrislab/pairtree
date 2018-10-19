import json
import numpy as np

def convert_adj_matrix_to_adj_list(adjm):
  adjm = np.copy(adjm)
  np.fill_diagonal(adjm, 0)
  adjl = {}

  for parent, child in zip(*np.nonzero(adjm)):
    # JSON keys must be strings.
    parent = str(parent)
    # Must convert from NumPy ints to Python ints. Ugh.
    child = int(child)
    if parent not in adjl:
      adjl[parent] = []
    adjl[parent].append(child)

  return adjl

def generate_treesumm(sampnames, clusters, adjmats, llh, phi):
  result = {
    'params': {'samples': sampnames},
    'trees': {},
  }

  assert len(adjmats) == len(llh)
  for tidx, adjmat, L in zip(range(len(llh)), adjmats, llh):
    pops = { str(pidx): {
      'num_ssms': len(ssms),
      'num_cnvs': 0,
      'cellular_prevalence': list(phi[tidx,pidx,:])
      } for pidx, ssms in enumerate(clusters)
    }
    result['trees'][str(tidx)] = {
      'llh': L,
      'populations': pops,
      'structure': convert_adj_matrix_to_adj_list(adjmat),
      'root': 0,
    }

  return result

def generate_mutlist(variants):
  mutlist = {'cnvs': {}, 'ssms': {}}

  varids = sorted(variants.keys(), key = lambda S: int(S[1:]))
  for varid in varids:
    var = variants[varid]
    mutlist['ssms'][varid] = {
      'name': var['name'],
      'ref_reads': list(var['ref_reads']),
      'total_reads': list(var['total_reads']),
      'expected_ref_in_ref': var['mu_r'],
      'expected_ref_in_variant': 1 - var['omega_v'],
    }

  return mutlist

def write_json(sampid, sampnames, variants, clusters, adjmats, llh, phi, treesummfn, mutlistfn):
  treesumm = generate_treesumm(sampnames, clusters, adjmats, llh, phi)
  mutlist = generate_mutlist(variants)

  for results in (treesumm, mutlist):
    results['dataset_name'] = sampid

  with open(treesummfn, 'w') as outf:
    json.dump(treesumm, outf)
  with open(mutlistfn, 'w') as outf:
    json.dump(mutlist, outf)
