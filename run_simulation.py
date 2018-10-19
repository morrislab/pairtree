import argparse
import numpy as np
import pickle

import common
import simulator
import pairwise
import vaf_plotter
import relation_plotter

from pairtree import write_header

def convert_to_array(pairs):
  N = int(np.sqrt(len(pairs)))
  assert len(pairs) == N**2
  arr = np.nan * np.ones((N, N, len(common.Models._all)))

  for pair in pairs.keys():
    V1, V2 = [int(V[1:]) for V in pair]
    arr[V1,V2,:] = pairs[pair]

  return arr

def write_results(data, posterior, evidence, datafn, pairwisefn, resultsfn):
  with open(datafn, 'wb') as outf:
    pickle.dump(data, outf)

  posterior = convert_to_array(posterior)
  evidence = convert_to_array(evidence)
  np.savez_compressed(pairwisefn, posterior=posterior, evidence=evidence)

  with open(resultsfn, 'w') as outf:
    sampid = 'A sample ID'
    write_header(sampid, outf)
    relation_plotter.plot_ml_relations(posterior, outf)
    relation_plotter.plot_relation_probs(posterior, outf)

    supervars = common.make_cluster_supervars(data['clusters'], data['variants_good'])
    garbage = {}
    vaf_plotter.plot_vaf_matrix(
      sampid,
      data['clusters'],
      data['variants_good'],
      supervars,
      data['variants_garbage'],
      data['phi'],
      data['sampnames'],
      None,
      False,
      outf)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--seed', dest='seed', type=int)
  parser.add_argument('--parallel', dest='parallel', type=int, default=1)
  parser.add_argument('-K', dest='K', type=int, default=4, help='Number of clusters')
  parser.add_argument('-S', dest='S', type=int, default=3, help='Number of samples')
  parser.add_argument('-T', dest='T', type=int, default=4000, help='Total reads per mutation')
  parser.add_argument('-M', dest='M', type=int, default=10, help='Number of non-garbage mutations')
  parser.add_argument('-G', dest='G', type=int, default=4, help='Number of garbage mutations')
  parser.add_argument('datafn')
  parser.add_argument('pairwisefn')
  parser.add_argument('resultsfn')
  args = parser.parse_args()

  if args.seed is None:
    seed = np.random.randint(2**31)
  else:
    seed = args.seed
  np.random.seed(args.seed)

  data = simulator.generate_data(args.K, args.S, args.T, args.M, args.G)
  data['seed'] = seed

  posterior, evidence = pairwise.calc_posterior(data['variants_all'], parallel=args.parallel, include_garbage_in_posterior=True, include_cocluster_in_posterior=True)
  write_results(data, posterior, evidence, args.datafn, args.pairwisefn, args.resultsfn)

if __name__ == '__main__':
  main()
