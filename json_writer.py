import json
import numpy as np
import common

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
  clusters = [[]] + clusters.tolist()
  assert len(adjmats) == len(llh) == len(phi)

  result = {
    'params': {'samples': sampnames},
    'trees': {},
  }

  for tidx, adjmat, L in zip(range(len(llh)), adjmats, llh):
    C, S = phi[tidx].shape
    assert C == len(adjmat) == len(clusters)
    assert S == len(sampnames)

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

  for varid in common.extract_vids(variants):
    var = variants[varid]
    mutlist['ssms'][varid] = {
      'name': var['name'],
      'ref_reads': var['ref_reads'].tolist(),
      'total_reads': var['total_reads'].tolist(),
      'var_read_prob': var['omega_v'].tolist(),
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
