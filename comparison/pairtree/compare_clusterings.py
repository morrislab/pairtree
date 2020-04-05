import argparse
import sklearn.metrics

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import common
import inputparser

def convert_clustering_to_assignment(clusters):
  mapping = {vid: cidx for cidx, cluster in enumerate(clusters) for vid in cluster}
  vids = common.sort_vids(mapping.keys())
  assign = [cidx for vid, cidx in mapping.items()]
  return (vids, assign)

def extract_assignment(paramsfn):
  params = inputparser.load_params(paramsfn)
  C = len(params['clusters'])
  vids, assign = convert_clustering_to_assignment(params['clusters'])
  return (C, vids, assign)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('params1_fn')
  parser.add_argument('params2_fn')
  args = parser.parse_args()

  C1, vids1, assign1 = extract_assignment(args.params1_fn)
  C2, vids2, assign2 = extract_assignment(args.params2_fn)
  assert vids1 == vids2

  homo, comp, vm = sklearn.metrics.homogeneity_completeness_v_measure(assign1, assign2)
  print(C1, C2, homo, comp, vm)

if __name__ == '__main__':
  main()
