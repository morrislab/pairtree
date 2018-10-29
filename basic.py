import argparse
import os
import numpy as np
import scipy.integrate
import warnings

import common
import pairwise
import vaf_plotter
import relation_plotter
import inputparser

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

def convert_to_array(pairs):
  N = int(np.sqrt(len(pairs)))
  assert len(pairs) == N**2
  arr = np.nan * np.ones((N, N, len(common.Models._all)))

  for pair in pairs.keys():
    V1, V2 = [int(V[1:]) for V in pair]
    arr[V1,V2,:] = pairs[pair]

  return arr

def lolpairwise(variants, sampid, parallel, pairwisefn):
  prior = {'garbage': 0.001}
  posterior, evidence = pairwise.calc_posterior(variants, prior=prior, parallel=parallel)
  posterior = convert_to_array(posterior)
  evidence = convert_to_array(evidence)
  np.savez_compressed(pairwisefn, posterior=posterior, evidence=evidence)

def main():
  np.set_printoptions(linewidth=400, precision=3, threshold=np.nan, suppress=True)
  np.seterr(divide='raise', invalid='raise')
  warnings.simplefilter('ignore', category=scipy.integrate.IntegrationWarning)

  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--seed', dest='seed', type=int)
  parser.add_argument('--parallel', dest='parallel', type=int, default=1)
  parser.add_argument('--params', dest='params_fn')
  parser.add_argument('sampid')
  parser.add_argument('ssm_fn')
  args = parser.parse_args()

  if args.seed is not None:
    np.random.seed(args.seed)

  variants = inputparser.load_ssms(args.ssm_fn)
  params = inputparser.load_params(args.params_fn)
  if 'samples' in params:
    sampnames = params['samples']
  else:
    sampnames = None

  #task = 'embed'
  task = 'calc'
  #task = 'plot'
  pairwisefn = '%s.npz' % args.sampid

  if task == 'embed':
    from IPython import embed
    embed()
  elif task == 'calc':
    lolpairwise(variants, args.sampid, args.parallel, pairwisefn)
  elif task == 'plot':
    if not os.path.exists(pairwisefn):
      import sys
      sys.exit()
    results = np.load(pairwisefn)
    with open('%s.results.html' % args.sampid, 'w') as outf:
      write_header(args.sampid, outf)
      relation_plotter.plot_ml_relations(results['posterior'], outf)
      relation_plotter.plot_relation_probs(results['posterior'], outf)

if __name__ == '__main__':
  main()
