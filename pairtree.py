import argparse
import numpy as np

import common
import handbuilt
import pairwise
import json_writer
import json
import tree_sampler
import phi_fitter

def create_matrix(model, model_probs, variants):
  N = len(variants)
  mat = np.zeros((N, N))
  for vids, P in model_probs.items():
    assert 0 <= P <= 1
    vidx1, vidx2 = [int(V[1:]) for V in vids.split(',')]
    mat[(vidx1, vidx2)] = P
  if model in ('cocluster', 'diff_branches'):
    # These should be symmetric.
    assert np.allclose(mat, mat.T)

  return mat

def create_model_prob_tensor(model_probs):
  M = len(model_probs['variants'])
  num_models = len(common.Models._all)
  tensor = np.zeros((M, M, num_models))
  for midx, mdl in enumerate(common.Models._all):
    tensor[:,:,midx] = create_matrix(mdl, model_probs['model_probs'][mdl], model_probs['variants'])
  return tensor

def fit_phis(adjm, variants, clusters, tidxs=None, parallel=1):
  ntrees = len(adjm)
  nsamples = len(list(variants.values())[0]['total_reads'])

  eta = np.ones((ntrees, nsamples, len(clusters)))
  phi = np.ones((ntrees, nsamples, len(clusters)))

  if tidxs is None:
    tidxs = range(ntrees)
  for tidx in tidxs:
    phi[tidx,:,:], eta[tidx,:,:] = phi_fitter.fit_phis(adjm[tidx], clusters, variants, parallel)

  return (phi, eta)

def write_header(sampid, outf):
  print('<script type="text/javascript" src="https://code.jquery.com/jquery-3.2.1.slim.min.js"></script>', file=outf)
  print('<script src="https://d3js.org/d3.v4.min.js"></script>', file=outf)
  print('<script src="https://d3js.org/d3-scale-chromatic.v1.min.js"></script>', file=outf)
  print('<script type="text/javascript" src="highlight_table_labels.js"></script>', file=outf)
  print('<script type="text/javascript" src="tree_plotter.js"></script>', file=outf)
  print('<link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet">', file=outf)
  print('<h1>%s</h1>' % sampid, file=outf)
  print('<link href="tree.css" rel="stylesheet">', file=outf)
  print('<style>td, th, table { padding: 5px; margin: 0; border-collapse: collapse; font-weight: normal; } span { visibility: hidden; } td:hover > span { visibility: visible; } .highlighted { background-color: black !important; color: white; }</style>', file=outf)

def write_trees(sampid, tidx, outf):
  colourings = {
    'SJBALL022609': [{'left': 'D', 'right': 'R1'}],
    'SJETV010steph': [{'left': 'D', 'right': 'R2'}],
    'SJBALL022610steph': [{'left': 'D', 'right': 'R1'}],
  }

  print('''<div id="trees"></div>''', file=outf)
  print('''<script type="text/javascript">$(document).ready(function() {''', file=outf)
  for colouring in colourings[sampid]:
    print('''(new TreePlotter()).plot('%s.summ.json', %s, '%s', '%s', '%s', '#trees');''' % (
      sampid,
      tidx,
      'sampled',
      colouring['left'],
      colouring['right'],
    ), file=outf)
  print('});</script>', file=outf)

def write_phi_matrix(sampid, outf):
  print('''<script type="text/javascript">$(document).ready(function() {
  (new PhiMatrix()).plot('%s', '%s.phi.json', '#phi_matrix');
  });</script>''' % (sampid, sampid), file=outf)
  print('<div id="phi_matrix" style="margin: 30px"></div>', file=outf)

def write_phi_json(phi, sampnames, phifn):
  with open(phifn, 'w') as outf:
    json.dump({
      'samples': sampnames,
      'phi': phi.tolist()
    }, outf)

def main():
  np.set_printoptions(threshold=np.nan, linewidth=120)
  np.seterr(divide='raise', invalid='raise')

  np.random.seed(1)

  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--struct-iterations', dest='struct_iterations', type=int, default=1000)
  parser.add_argument('--parallel', dest='parallel', type=int, default=1)
  parser.add_argument('sampid')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('clusters_fn')
  parser.add_argument('tree_type')
  args = parser.parse_args()

  variants = common.parse_ssms(args.sampid, args.ssm_fn)
  sampnames = common.load_sampnames(args.params_fn)
  clusters = handbuilt.load_clusters(args.clusters_fn, variants, args.tree_type, sampnames)

  supervars = common.make_cluster_supervars(clusters, variants)
  superclusters = common.make_superclusters(supervars)

  posterior, evidence = pairwise.calc_posterior(supervars, parallel=args.parallel)
  results = pairwise.generate_results(posterior, evidence, supervars)

  model_probs_tensor = create_model_prob_tensor(results)
  sampled_adjm, sampled_llh = tree_sampler.sample_trees(model_probs_tensor, superclusters, nsamples=args.struct_iterations, nchains=args.parallel)
  phi, eta = fit_phis(sampled_adjm, supervars, superclusters, tidxs=(-1,), parallel=args.parallel)
  json_writer.write_json(
    args.sampid,
    sampnames,
    variants,
    clusters,
    sampled_adjm,
    sampled_llh,
    phi,
    '%s.summ.json' % args.sampid,
    '%s.muts.json' % args.sampid,
  )
  write_phi_json(phi[-1].T, sampnames, '%s.phi.json' % args.sampid)

  with open('%s.results.html' % args.sampid, 'w') as outf:
    write_header(args.sampid, outf)
    write_trees(args.sampid, len(sampled_adjm) - 1, outf)
    write_phi_matrix(args.sampid, outf)

main()
