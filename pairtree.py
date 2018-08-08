import argparse
import os
import numpy as np

import common
import handbuilt
import pairwise
import json_writer
import json
import tree_sampler
import phi_fitter
import vaf_plotter

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

def create_model_prob_tensor(model_probs, K='model_probs'):
  M = len(model_probs['variants'])
  num_models = len(common.Models._all)
  tensor = np.zeros((M, M, num_models))
  for midx, mdl in enumerate(common.Models._all):
    tensor[:,:,midx] = create_matrix(mdl, model_probs[K][mdl], model_probs['variants'])
  return tensor

def fit_phis(adjm, variants, clusters, tidxs, iterations, parallel=1):
  ntrees = len(adjm)
  nsamples = len(list(variants.values())[0]['total_reads'])

  N, K, S = ntrees, len(clusters), nsamples
  eta = np.ones((N, K, S))
  phi = np.ones((N, K, S))

  for tidx in tidxs:
    phi[tidx,:,:], eta[tidx,:,:] = phi_fitter.fit_phis(adjm[tidx], clusters, variants, iterations, parallel)

  return (phi, eta)

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
  print('<style type="text/css">td, th, table { padding: 5px; margin: 0; border-collapse: collapse; font-weight: normal; } span { visibility: hidden; } td:hover > span { visibility: visible; } .highlighted { background-color: black !important; color: white; }</style>', file=outf)

def write_trees(sampid, tidx, outf):
  colourings = {
    'SJBALL022609': [{'left': 'D', 'right': 'R1'}],
    'SJETV010steph': [{'left': 'D', 'right': 'R2'}],
    'SJBALL022610steph': [{'left': 'D', 'right': 'R1'}],
  }
  if sampid not in colourings:
    colourings[sampid] = [{'left': 'D', 'right': 'R1'}]

  print('''<div id="trees"></div>''', file=outf)
  print('''<script type="text/javascript">$(document).ready(function() {''', file=outf)
  print('new VafMatrix("#vafmatrix_toggles");', file=outf)
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

def print_error(phi, supervars, llh):
  svkeys = sorted(supervars.keys(), key = lambda S: int(S[1:]))
  supervar_vaf = np.array([supervars[S]['vaf'] for S in svkeys])
  assert np.allclose(phi[0], 1)
  phi = phi[1:] # Discard phi for normal node, which should always be 1.
  assert supervar_vaf.shape == phi.shape

  diff = np.abs(2*supervar_vaf - phi)
  M, S = phi.shape
  phi_error = np.sum(diff) / (M*S)

  llh_error = -llh / (M*S*np.log(2))
  print('phi_error = %.3f, llh_error = %.3f' % (phi_error, llh_error))

def remove_garbage(garbage_ids, variants):
  garbage_variants = {}
  for varid in sorted(variants.keys(), key = lambda S: int(S[1:])):
    if int(varid[1:]) in garbage_ids:
      garbage_variants[varid] = variants[varid]
      del variants[varid]
  return garbage_variants

def main():
  np.set_printoptions(threshold=np.nan, linewidth=120)
  np.seterr(divide='raise', invalid='raise')
  np.random.seed(1)

  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--phi-iterations', dest='phi_iterations', type=int, default=1000)
  parser.add_argument('--tree-samples', dest='tree_samples', type=int, default=1000)
  parser.add_argument('--tree-chains', dest='tree_chains', type=int, default=1)
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
  garbage_vids = handbuilt.load_garbage(args.clusters_fn, args.tree_type)
  garbage_variants = remove_garbage(garbage_vids, variants)

  supervars = common.make_cluster_supervars(clusters, variants)
  superclusters = common.make_superclusters(supervars)

  adjm = handbuilt.load_tree(args.clusters_fn, args.tree_type)
  if adjm is not None:
    sampled_adjm = [adjm]
    sampled_llh = [0]
  else:
    posterior, evidence = pairwise.calc_posterior(supervars, parallel=args.parallel, include_garbage_in_posterior=False, include_cocluster_in_posterior=False)
    results = pairwise.generate_results(posterior, evidence, supervars)

    model_probs_tensor = create_model_prob_tensor(results)
    sampled_adjm, sampled_llh = tree_sampler.sample_trees(model_probs_tensor, supervars, superclusters, nsamples=args.tree_samples, nchains=args.tree_chains, parallel=args.parallel)

  phi, eta = fit_phis(sampled_adjm, supervars, superclusters, tidxs=(-1,), iterations=args.phi_iterations, parallel=args.parallel)
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
  write_phi_json(phi[-1], sampnames, '%s.phi.json' % args.sampid)

  with open('%s.results.html' % args.sampid, 'w') as outf:
    write_header(args.sampid, outf)
    write_trees(args.sampid, len(sampled_adjm) - 1, outf)
    write_phi_matrix(args.sampid, outf)
    vaf_plotter.plot_vaf_matrix(args.sampid, clusters, variants, supervars, garbage_variants, phi[-1], sampnames, None, False, outf)

  print_error(phi[-1], supervars, sampled_llh[-1])

main()
