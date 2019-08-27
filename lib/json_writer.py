import json
import numpy as np
import common

def generate_treesumm(sampnames, clusters, adjmats, llh, phi):
  clusters = list(clusters)
  clusters = [[]] + clusters
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
      'structure': common.convert_adj_matrix_to_json_adjlist(adjmat),
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
