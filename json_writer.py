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

def generate_treesumm(clusters, nsamples, adjmats, llh, handbuilt_adjm):
  result = {
    'trees': {},
  }
  pops = { str(pidx): {
    'num_ssms': len(ssms),
    'num_cnvs': 0,
    'cellular_prevalence': nsamples*[1] }
    for pidx, ssms in enumerate(clusters)
  }

  assert len(adjmats) == len(llh)
  for tidx, adjmat, llh in zip(range(len(llh)), adjmats, llh):
    result['trees'][str(tidx)] = {
      'llh': llh,
      'populations': pops,
      'structure': convert_adj_matrix_to_adj_list(adjmat),
      'root': 0,
    }

  if handbuilt_adjm is not None:
    result['trees']['-1'] = {
      'llh': 0,
      'populations': pops,
      'structure': convert_adj_matrix_to_adj_list(handbuilt_adjm),
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
      'expected_ref_in_ref': 0.999,     # Dummy value
      'expected_ref_in_variant': 0.499, # Dummy value
    }

  return mutlist

def write_json(sampid, clusters, adjmats, llh, handbuilt_adjm, variants, treesummfn, mutlistfn):
  nsamples = len(list(variants.values())[0]['total_reads'])
  treesumm = generate_treesumm(clusters, nsamples, adjmats, llh, handbuilt_adjm)
  mutlist = generate_mutlist(variants)

  for results in (treesumm, mutlist):
    results['dataset_name'] = sampid

  with open(treesummfn, 'w') as outf:
    json.dump(treesumm, outf)
  with open(mutlistfn, 'w') as outf:
    json.dump(mutlist, outf)
