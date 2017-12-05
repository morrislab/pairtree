import argparse
import json
from collections import OrderedDict
import common
import handbuilt

def calc_jaccard(S1, S2):
  S1, S2 = set(S1), set(S2)
  if len(S1 | S2) == 0:
    return 0
  return float(len(S1 & S2)) / len(S1 | S2)

def compare(clusters1, clusters2):
  closest = OrderedDict()
  for J, cluster_j in enumerate(clusters1):
    max_jaccard = float('-inf')
    for K, cluster_k in enumerate(clusters2):
      jac = calc_jaccard(cluster_j, cluster_k)
      if jac > max_jaccard:
        max_jaccard = jac
        closest[J] = (K, cluster_j, cluster_k, jac)
    assert max_jaccard >= 0
  return dict(closest)

def load_sampnames(paramsfn):
  with open(paramsfn) as P:
    params = json.load(P)
  sampnames = params['samples']
  return sampnames

def load_clusters(sampid, tree_type, handbuiltfn, ssmfn, paramsfn):
  sampnames = load_sampnames(paramsfn)
  variants = common.parse_ssms(sampid, ssmfn)
  if tree_type == 'handbuilt.patient':
    variants, sampnames = common.extract_patient_samples(variants, sampnames)
  clusters, handbuilt_adjm, node_colourings = handbuilt.load_clusters_and_tree(handbuiltfn, variants, tree_type, sampnames)
  return clusters

def print_results(closest):
  for J, (K, cluster_j, cluster_k, jac) in closest.items():
    print(J, K, jac, len(cluster_j), len(cluster_k), sorted(cluster_j), sorted(cluster_k), sep='\t')

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sampid')
  parser.add_argument('handbuiltfn')
  parser.add_argument('ssmfn')
  parser.add_argument('paramsfn')
  parser.add_argument('tree_type_1')
  parser.add_argument('tree_type_2')
  args = parser.parse_args()

  clusters1 = load_clusters(args.sampid, args.tree_type_1, args.handbuiltfn, args.ssmfn, args.paramsfn)
  clusters2 = load_clusters(args.sampid, args.tree_type_2, args.handbuiltfn, args.ssmfn, args.paramsfn)
  closest = compare(clusters1, clusters2)
  print_results(closest)

if __name__ == '__main__':
  main()
