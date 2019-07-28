import handbuilt
import common
import numpy as np
import plot
import json
import argparse
import os
import random

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
  rand = random.randint(0, 2**32)
  elemid = 'cluster_matrix_%s' % rand
  return '''<script type="text/javascript">$(document).ready(function() {
  (new ClusterMatrix()).plot('%s', '#%s');
  });</script>
  <div id="%s" style="margin: 30px"></div>''' % (clustermatfn, elemid, elemid)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--sampid', dest='sampid', required=True)
  parser.add_argument('--ssms', dest='ssm_fn', required=True)
  parser.add_argument('--params', dest='params_fn', required=True)
  parser.add_argument('--handbuilt1', dest='handbuilt_fn1', required=True)
  parser.add_argument('--handbuilt2', dest='handbuilt_fn2', required=True)
  parser.add_argument('--treetype1', dest='treetype1', required=True)
  parser.add_argument('--treetype2', dest='treetype2', required=True)
  parser.add_argument('clustermat_fn')
  args = parser.parse_args()

  variants = common.parse_ssms(args.sampid, args.ssm_fn)
  sampnames = common.load_sampnames(args.params_fn)
  clusters1 = handbuilt.load_clusters(args.handbuilt_fn1, variants, args.treetype1, sampnames)
  clusters2 = handbuilt.load_clusters(args.handbuilt_fn2, variants, args.treetype2, sampnames)
  garbage1 = handbuilt.load_garbage(args.handbuilt_fn1, args.treetype1)
  garbage2 = handbuilt.load_garbage(args.handbuilt_fn2, args.treetype2)

  clustermat = make_clustermat(clusters1, garbage1, clusters2, garbage2)
  write_clustermat_json(clustermat, args.clustermat_fn)
  html = make_clustermat_html(os.path.basename(args.clustermat_fn))
  print(html)

if __name__ == '__main__':
  main()
