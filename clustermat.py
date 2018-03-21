import handbuilt
import common
import numpy as np
import plot
import json

def write_cluster_matrix(sampid, outf):
  print('''<script type="text/javascript">$(document).ready(function() {
  (new ClusterMatrix()).plot('%s.handbuilt.json', 'xeno', 'patient', '#cluster_matrix');
  });</script>''' % sampid, file=outf)
  print('<div id="cluster_matrix" style="margin: 30px"></div>', file=outf)

def load_xeno_and_patient_clusters(sampid, handbuiltfn, paramsfn, ssmfn):
  xeno = {
    'variants':  common.parse_ssms(sampid, ssmfn),
    'sampnames': plot.load_sampnames(paramsfn),
    'tree_type': 'handbuilt.xeno',
  }
  patient = {
    'tree_type': 'handbuilt.patient'
  }
  patient['variants'], patient['sampnames'] = common.extract_patient_samples(xeno['variants'], xeno['sampnames'])

  xeno['garbage']    = handbuilt.load_garbage(handbuiltfn, xeno['tree_type'])
  patient['garbage'] = handbuilt.load_garbage(handbuiltfn, patient['tree_type'])

  xeno['clusters'], _, _    = handbuilt.load_clusters_and_tree(handbuiltfn, xeno['variants'], xeno['tree_type'], xeno['sampnames'])
  patient['clusters'], _, _ = handbuilt.load_clusters_and_tree(handbuiltfn, patient['variants'], patient['tree_type'], patient['sampnames'])

  assert len(xeno['clusters'][0]) == len(patient['clusters'][0]) == 0
  xeno['clusters'][0] = xeno['garbage']
  patient['clusters'][0] = patient['garbage']
  assert sum([len(C) for C in xeno['clusters']]) == sum([len(C) for C in patient['clusters']])
  return (xeno['clusters'], patient['clusters'])

def make_clustermat(xeno_clusters, patient_clusters):
  M, N = len(xeno_clusters), len(patient_clusters)
  clustermat = np.zeros((M, N))
  T = sum([len(C) for C in patient_clusters])
  for i in range(M):
    for j in range(N):
      clustermat[i,j] = len(set(xeno_clusters[i]) & set(patient_clusters[j]))
  assert np.isclose(T, np.sum(clustermat))

  clustermat /= T
  # Make each row sum to 1, showing distribution of its assignments. Rather
  # than vectorizing operation, do row-wise so that we don't divide by zero on
  # empty rows. (But the only empty row should be the garbage row.)
  for idx, row in enumerate(clustermat):
    if np.sum(row) == 0:
      # This should be the garbage row -- no actual cluster should be empty.
      assert idx == 0
      continue
    clustermat[idx] = row / np.sum(row)
  return clustermat

def write_clustermat_json(sampid, handbuiltfn, paramsfn, ssmfn, clustermatfn):
  xeno_clusters, patient_clusters = load_xeno_and_patient_clusters(sampid, handbuiltfn, paramsfn, ssmfn)
  clustermat = make_clustermat(xeno_clusters, patient_clusters)

  with open(clustermatfn, 'w') as outf:
    json.dump({
      'clustermat': clustermat.tolist()
    }, outf)
