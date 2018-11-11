import argparse
import os
import numpy as np
import scipy.integrate
import warnings
import random
import multiprocessing

import common
import pairwise
import vaf_plotter
import relation_plotter
import inputparser
import json_writer
import tree_sampler
import phi_fitter
import clustermaker
import resultserializer

def read_file(fn):
  basedir = os.path.abspath(os.path.dirname(__file__))
  with open(os.path.join(basedir, fn)) as F:
    return F.read()

def write_header(sampid, outf):
  print('<script type="text/javascript" src="https://code.jquery.com/jquery-3.2.1.slim.min.js"></script>', file=outf)
  print('<script src="https://d3js.org/d3.v4.min.js"></script>', file=outf)
  print('<script src="https://d3js.org/d3-scale-chromatic.v1.min.js"></script>', file=outf)
  print('<script type="text/javascript">%s</script>' % read_file('highlight_table_labels.js'), file=outf)
  print('<script type="text/javascript">%s</script>' % read_file('tree_plotter.js'), file=outf)
  print('<link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet">', file=outf)
  print('<h1>%s</h1>' % sampid, file=outf)
  print('<style type="text/css">%s</style>' % read_file('tree.css'), file=outf)
  print('<style type="text/css">%s</style>' % read_file('matrix.css'), file=outf)
  print('<style type="text/css">td, th, table { padding: 5px; margin: 0; border-collapse: collapse; font-weight: normal; } span { visibility: hidden; } td:hover > span { visibility: visible; } .highlighted { background-color: black !important; color: white; }</style>', file=outf)

def lolsomething(variants, prior, args):
  #task = 'embed'
  task = 'calc'
  #task = 'plot'
  pairwisefn = '%s.npz' % args.sampid

  if task == 'embed':
    from IPython import embed
    embed()
  elif task == 'calc':
    posterior, evidence = pairwise.calc_posterior(variants, prior=prior, parallel=parallel)
    posterior = convert_to_array(posterior)
    evidence = convert_to_array(evidence)
    np.savez_compressed(pairwisefn, posterior=posterior, evidence=evidence)
  elif task == 'plot':
    if not os.path.exists(pairwisefn):
      import sys
      sys.exit()
    results = resultserializer.load(pairwisefn)
    with open('%s.results.html' % args.sampid, 'w') as outf:
      write_header(args.sampid, outf)
      relation_plotter.plot_ml_relations(results['posterior'], outf)
      relation_plotter.plot_relation_probs(results['posterior'], outf)
  import sys
  sys.exit()

def fit_phis(adjm, supervars, superclusters, tidxs, iterations, parallel):
  ntrees = len(adjm)
  nsamples = len(next(iter(supervars.values()))['total_reads'])

  N, K, S = ntrees, len(superclusters), nsamples
  eta = np.ones((N, K, S))
  phi = np.ones((N, K, S))

  for tidx in tidxs:
    phi[tidx,:,:], eta[tidx,:,:] = phi_fitter.fit_phis(
      adjm[tidx],
      superclusters,
      supervars,
      iterations,
      parallel
    )

  return (phi, eta)

def load_sampnames(params):
  if 'samples' in params:
    return params['samples']
  else:
    return None

def main():
  np.set_printoptions(linewidth=400, precision=3, threshold=np.nan, suppress=True)
  np.seterr(divide='raise', invalid='raise')
  warnings.simplefilter('ignore', category=scipy.integrate.IntegrationWarning)

  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--seed', dest='seed', type=int)
  parser.add_argument('--parallel', dest='parallel', type=int, default=None)
  parser.add_argument('--params', dest='params_fn')
  parser.add_argument('--trees-per-chain', dest='trees_per_chain', type=int, default=1000)
  parser.add_argument('--tree-chains', dest='tree_chains', type=int, default=None)
  parser.add_argument('--phi-iterations', dest='phi_iterations', type=int, default=1000)
  parser.add_argument('sampid')
  parser.add_argument('ssm_fn')
  args = parser.parse_args()

  # Note that multiprocessing.cpu_count() returns number of logical cores, so
  # if you're using a hyperthreaded CPU, this will be more than the number of
  # physical cores you have.
  parallel = args.parallel if args.parallel is not None else multiprocessing.cpu_count()
  tree_chains = args.tree_chains if args.tree_chains is not None else parallel
  prior = {'garbage': 0.001}

  if args.seed is not None:
    np.random.seed(args.seed)
    random.seed(args.seed)

  variants = inputparser.load_ssms(args.ssm_fn)
  params = inputparser.load_params(args.params_fn)
  sampnames = load_sampnames(params)

  resultsfn = args.sampid + '.results.npz'
  import os
  if os.path.exists(resultsfn):
    results = resultserializer.load(resultsfn)
  else:
    results = {}

  if 'mutrel_posterior' not in results:
    results['mutrel_posterior'], results['mutrel_evidence'] = pairwise.calc_posterior(variants, prior=prior, parallel=parallel)

  if 'clustrel_posterior' not in results:
    if 'clusters' in params and 'garbage' in params:
      supervars, results['clustrel_posterior'], results['clustrel_evidence'], results['clusters'], results['garbage'] = clustermaker.use_pre_existing(
        variants,
        results['mutrel_posterior'],
        prior,
        parallel,
        params['clusters'],
        params['garbage'],
      )
    else:
      supervars, results['clustrel_posterior'], results['clustrel_evidence'], results['clusters'], results['garbage'] = clustermaker.cluster_and_discard_garbage(
        variants,
        results['mutrel_posterior'],
        prior,
        parallel,
      )
  else:
    supervars = clustermaker.make_cluster_supervars(results['clusters'], variants)

  superclusters = clustermaker.make_superclusters(supervars)
  # Add empty initial cluster, which serves as tree root.
  superclusters.insert(0, [])

  if 'adjm' not in results:
    results['adjm'], results['llh'] = tree_sampler.sample_trees(
      results['clustrel_posterior'],
      supervars,
      superclusters,
      trees_per_chain = args.trees_per_chain,
      nchains = tree_chains,
      parallel = parallel,
    )

  if 'phi' not in results:
    results['phi'], results['eta'] = fit_phis(
      results['adjm'],
      supervars,
      superclusters,
      tidxs = (-1,),
      iterations = args.phi_iterations,
      parallel = parallel,
    )

  #json_writer.write_json(
  #  args.sampid,
  #  sampnames,
  #  variants,
  #  clusters,
  #  results['adjm'],
  #  results['llh'],
  #  results['phi'],
  #  '%s.summ.json' % args.sampid,
  #  '%s.muts.json' % args.sampid,
  #)

  resultserializer.save(results, resultsfn)
  #lolsomething(variants, prior, args)

if __name__ == '__main__':
  main()
