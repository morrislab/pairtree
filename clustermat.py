import handbuilt
import common
import numpy as np
import plot
import json

def make_clustermat(clusters1, garbage1, clusters2, garbage2):
  assert len(clusters1[0]) == len(clusters2[0]) == 0
  clusters1[0] = garbage1
  clusters2[0] = garbage2
  assert sum([len(C) for C in clusters1]) == sum([len(C) for C in clusters2])

  M, N = len(clusters1), len(clusters2)
  clustermat = np.zeros((M, N))
  T = sum([len(C) for C in clusters2])
  for i in range(M):
    for j in range(N):
      clustermat[i,j] = len(set(clusters1[i]) & set(clusters2[j]))
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

def write_clustermat_json(clustermat, clustermatfn):
  with open(clustermatfn, 'w') as outf:
    json.dump({
      'clustermat': clustermat.tolist()
    }, outf)

def make_clustermat_html(clustermatfn):
  return '''<script type="text/javascript">$(document).ready(function() {
  (new ClusterMatrix()).plot('%s', '#cluster_matrix');
  });</script>
  <div id="cluster_matrix" style="margin: 30px"></div>''' % clustermatfn
