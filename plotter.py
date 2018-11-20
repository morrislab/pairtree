import argparse
import os

import common
import resultserializer
import clustermaker
import inputparser

import vaf_plotter
import relation_plotter
import json_writer

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

# TODO: eliminate the phi.json file, and instead make the JS extract it from
# the tree summ.
def write_phi_json(phi, sampnames, phifn):
  import json
  with open(phifn, 'w') as outf:
    json.dump({
      'samples': sampnames,
      'phi': phi.tolist()
    }, outf)

def write_phi_matrix(sampid, outf):
  print('''<script type="text/javascript">$(document).ready(function() {
  (new PhiMatrix()).plot('%s', '%s.phi.json', '#phi_matrix');
  });</script>''' % (sampid, sampid), file=outf)
  print('<div id="phi_matrix" style="margin: 30px"></div>', file=outf)

def write_trees(sampid, tidx, outf):
  colourings = [{'left': 'D', 'right': 'R1'}]

  print('''<div id="trees"></div>''', file=outf)
  print('''<script type="text/javascript">$(document).ready(function() {''', file=outf)
  for colouring in colourings:
    print('''(new TreePlotter()).plot('%s.summ.json', %s, '%s', '%s', '%s', '#trees');''' % (
      sampid,
      tidx,
      'sampled',
      colouring['left'],
      colouring['right'],
    ), file=outf)
  print('});</script>', file=outf)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sampid')
  parser.add_argument('ssm_fn')
  parser.add_argument('params_fn')
  parser.add_argument('results_fn')
  parser.add_argument('out_fn')
  args = parser.parse_args()

  results = resultserializer.load(args.results_fn)
  variants = inputparser.load_ssms(args.ssm_fn)
  params = inputparser.load_params(args.params_fn)
  supervars = clustermaker.make_cluster_supervars(results['clusters'], variants)
  supervars = [supervars[vid] for vid in common.sort_vids(supervars.keys())]

  if 'phi' in results:
    json_writer.write_json(
      args.sampid,
      params['samples'],
      variants,
      results['clusters'],
      results['adjm'],
      results['llh'],
      results['phi'],
      '%s.summ.json' % args.sampid,
      '%s.muts.json' % args.sampid,
    )
    write_phi_json(results['phi'][-1], params['samples'], '%s.phi.json' % args.sampid)

  with open(args.out_fn, 'w') as outf:
    write_header(args.sampid, outf)
    if 'adjm' in results:
      write_trees(args.sampid, len(results['adjm']) - 1, outf)
    if 'phi' in results:
      write_phi_matrix(args.sampid, outf)
    relation_plotter.plot_ml_relations(results['mutrel_posterior'], outf)
    relation_plotter.plot_separate_relations(results['mutrel_posterior'], outf)
    vaf_plotter.plot_vaf_matrix(
      results['clusters'],
      variants,
      supervars,
      results['garbage'],
      results['phi'][-1] if 'phi' in results else None,
      params['samples'],
      should_correct_vaf=True,
      outf=outf,
    )

if __name__ == '__main__':
  main()
