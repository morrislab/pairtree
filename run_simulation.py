import argparse
import numpy as np
import pickle
from collections import defaultdict

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

def make_variants(data):
  variants = {}
  for midx in range(len(data['omega_v'])):
    variant = {
      'id': 's%s' % midx,
      'name': 'Happy variant %s' % midx,
      'var_reads': data['V'][midx],
      'total_reads': data['T'][midx],
      'omega_v': data['omega_v'][midx],
      'vaf_correction': 1.,
    }
    variant['ref_reads'] = variant['total_reads'] - variant['var_reads']
    variant['vaf'] = variant['var_reads'] / variant['total_reads']
    variants[variant['id']] = variant
  return variants

def make_clusters(mutass):
  clusters = defaultdict(list)
  for midx, cidx in enumerate(mutass):
    clusters[cidx].append(midx)
  assert set(clusters.keys()) == set(range(len(clusters)))

  clusters = [clusters[cidx] for cidx in sorted(clusters.keys())]
  return clusters

def write_results(variants, data, posterior, datafn, resultsfn):
  with open(datafn, 'wb') as outf:
    pickle.dump(data, outf)

  with open(resultsfn, 'w') as outf:
    sampid = 'A sample ID'
    write_header(sampid, outf)
    posterior = convert_to_array(posterior)
    relation_plotter.plot_ml_relations(posterior, outf)
    relation_plotter.plot_relation_probs(posterior, outf)

    # TODO: remove empty first cluster.
    clusters = [[]] + make_clusters(data['mutass'])
    data['phi'] = np.insert(data['phi'], 0, 1, axis=0)

    supervars = common.make_cluster_supervars(clusters, variants)
    garbage = {}
    vaf_plotter.plot_vaf_matrix(
      sampid,
      clusters,
      variants,
      supervars,
      garbage,
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
  parser.add_argument('datafn')
  parser.add_argument('resultsfn')
  args = parser.parse_args()
  if args.seed is not None:
    np.random.seed(args.seed)

  K = 4
  S = 3
  T = 4000
  M = 10
  data = simulator.generate_data(K, S, T, M)
  variants = make_variants(data)

  posterior, evidence = pairwise.calc_posterior(variants, parallel=args.parallel, include_garbage_in_posterior=True, include_cocluster_in_posterior=True)
  write_results(variants, data, posterior, args.datafn, args.resultsfn)

if __name__ == '__main__':
  main()
