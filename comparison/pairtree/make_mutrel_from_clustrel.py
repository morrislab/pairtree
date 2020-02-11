import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import resultserializer
import evalutil
from mutrel import Mutrel
from common import Models, NUM_MODELS

def perturb_clustrel(clustrel):
  K = len(clustrel.rels)
  new_rels = np.zeros((K+1, K+1, NUM_MODELS))
  new_rels[0,0,Models.cocluster] = 1
  new_rels[0,1:,Models.A_B] = 1
  new_rels[1:,0,Models.B_A] = 1
  new_rels[1:,1:,:] = clustrel.rels
  new_vids = list(clustrel.vids) + ['S%s' % len(clustrel.vids)]
  return Mutrel(vids=new_vids, rels=new_rels)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('pairtree_results_fn')
  parser.add_argument('clustrel_mutrel_fn')
  args = parser.parse_args()

  results = resultserializer.load(args.pairtree_results_fn)
  clusters = [[]] + list(results['clusters'])
  garbage = list(results['garbage'])
  all_vids = set([V for C in results['clusters'] for V in C] + garbage)

  clustrel = perturb_clustrel(results['clustrel_posterior'])
  clustrel_mutrel = evalutil.make_mutrel_from_clustrel(clustrel, clusters)
  clustrel_mutrel = evalutil.add_garbage(clustrel_mutrel, garbage)
  assert set(clustrel_mutrel.vids) == all_vids
  evalutil.save_sorted_mutrel(clustrel_mutrel, args.clustrel_mutrel_fn)

if __name__ == '__main__':
  main()
