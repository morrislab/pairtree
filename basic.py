import argparse
import os
import numpy as np
import scipy.integrate
import warnings
import random
import multiprocessing
import sys

import common
import pairwise
import inputparser
import tree_sampler
import clustermaker
import resultserializer
import plotter

def main():
  np.set_printoptions(linewidth=400, precision=3, threshold=sys.maxsize, suppress=True)
  np.seterr(divide='raise', invalid='raise')
  warnings.simplefilter('ignore', category=scipy.integrate.IntegrationWarning)

  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--verbose', action='store_true')
  parser.add_argument('--seed', dest='seed', type=int)
  parser.add_argument('--parallel', dest='parallel', type=int, default=None)
  parser.add_argument('--params', dest='params_fn')
  parser.add_argument('--burnin-per-chain', dest='burnin_per_chain', type=int, default=1000)
  parser.add_argument('--trees-per-chain', dest='trees_per_chain', type=int, default=2000)
  parser.add_argument('--tree-chains', dest='tree_chains', type=int, default=None)
  parser.add_argument('--phi-iterations', dest='phi_iterations', type=int, default=10000)
  parser.add_argument('--phi-fitter', dest='phi_fitter', choices=('graddesc', 'projection', 'rprop'), default='graddesc')
  parser.add_argument('--tree_perturbations', dest='tree_perturbations', type=int, default=100)
  parser.add_argument('ssm_fn')
  parser.add_argument('results_fn')
  args = parser.parse_args()

  common.debug.DEBUG = args.verbose

  # Note that multiprocessing.cpu_count() returns number of logical cores, so
  # if you're using a hyperthreaded CPU, this will be more than the number of
  # physical cores you have.
  parallel = args.parallel if args.parallel is not None else multiprocessing.cpu_count()
  tree_chains = args.tree_chains if args.tree_chains is not None else parallel
  prior = {'garbage': 0.001}

  if args.seed is not None:
    seed = args.seed
  else:
    # Maximum seed is 2**32 - 1.
    seed = np.random.randint(2**32)
  np.random.seed(seed)
  random.seed(seed)

  variants = inputparser.load_ssms(args.ssm_fn)
  params = inputparser.load_params(args.params_fn)

  if os.path.exists(args.results_fn):
    results = resultserializer.load(args.results_fn)
  else:
    results = {'seed': seed}

  if 'clustrel_posterior' not in results:
    if 'clusters' in params and 'garbage' in params:
      supervars, results['clustrel_posterior'], results['clustrel_evidence'], results['clusters'], results['garbage'] = clustermaker.use_pre_existing(
        variants,
        prior,
        parallel,
        params['clusters'],
        params['garbage'],
      )
    else:
      if 'mutrel_posterior' not in results:
        results['mutrel_posterior'], results['mutrel_evidence'] = pairwise.calc_posterior(variants, prior=prior, rel_type='variant', parallel=parallel)
        resultserializer.save(results, args.results_fn)
      clustermaker._plot.prefix = os.path.basename(args.ssm_fn).split('.')[0]
      supervars, results['clustrel_posterior'], results['clustrel_evidence'], results['clusters'], results['garbage'] = clustermaker.cluster_and_discard_garbage(
        variants,
        results['mutrel_posterior'],
        results['mutrel_evidence'],
        prior,
        parallel,
      )
    resultserializer.save(results, args.results_fn)
  else:
    supervars = clustermaker.make_cluster_supervars(results['clusters'], variants)

  superclusters = clustermaker.make_superclusters(supervars)
  # Add empty initial cluster, which serves as tree root.
  superclusters.insert(0, [])

  if 'adjm' not in results:
    if 'structure' not in params:
      results['adjm'], results['phi'], results['llh'] = tree_sampler.sample_trees(
        results['clustrel_posterior'],
        supervars,
        superclusters,
        args.trees_per_chain,
        args.burnin_per_chain,
        tree_chains,
        args.phi_fitter,
        args.phi_iterations,
        args.tree_perturbations,
        seed,
        parallel,
      )
      resultserializer.save(results, args.results_fn)
    else:
      adjlist = inputparser.load_structure(params['structure'])
      adjm = common.convert_adjlist_to_adjmatrix(adjlist)
      results['adjm'], results['phi'], results['llh'] = tree_sampler.use_existing_structure(
        adjm,
        supervars,
        superclusters,
        args.phi_fitter,
        args.phi_iterations,
        parallel
      )
      resultserializer.save(results, args.results_fn)

if __name__ == '__main__':
  main()
