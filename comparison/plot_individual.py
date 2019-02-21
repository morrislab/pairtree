import pandas as pd
import argparse
import plotly
import plotly.figure_factory as ff
from collections import defaultdict
import re
import numpy as np
import sys

SIM_PARAMS = 'GKMST'
SIM_PARAM_LABELS = {
  'G': 'Garbage mutations',
  'K': 'Number of clusters',
  'M': 'Non-garbage mutations',
  'S': 'Samples',
  'T': 'Read depth',
}
HAPPY_METHOD_NAMES = {
  'mle_unconstrained': 'MLE lineage frequencies',
  'pairtree_trees_llh': 'Pairtree (automated tree search)',
  'pairtree_handbuilt': 'Pairtree (manually constructed trees)',
  'pwgs_allvars_single_llh': 'PhyloWGS (no clustering enforced)',
  'pwgs_allvars_multi_llh': 'Multi-chain PhyloWGS (no clustering enforced)',
  'pwgs_supervars_single_llh': 'PhyloWGS (enforced clustering)',
  'pwgs_supervars_multi_llh': 'Multi-chain PhyloWGS (enforced clustering)',
}

def augment(results):
  if len(SIM_PARAMS) < 1:
    return results
  simparams = defaultdict(list)
  sim_param_names = ''.join(SIM_PARAMS)
  for rid in results.runid:
    for token in re.findall(rf'([{sim_param_names}]\d+)', rid):
      K, V = token[0], token[1:]
      simparams[K].append(int(V))

  lengths = np.array([len(V) for V in simparams.values()])
  if not np.all(lengths == len(results)):
    raise Exception('Some SIM_PARAMS missing from some points')

  for K in simparams.keys():
    results[K] = simparams[K]
  return results

def load_results(resultsfn):
  results = pd.read_csv(resultsfn)
  #results = results.drop(['pairtree_clustrel'], axis=1)
  results = augment(results)
  methods = get_method_names(results)

  for M in methods:
    inf_idxs = np.isinf(results[M])
    if np.any(inf_idxs):
      # Sometimes the LLH may be -inf. (Thanks, PASTRI.) Consider this to be a
      # failed run.
      print('%s has %s infs' % (M, np.sum(inf_idxs)), file=sys.stderr)
      results.loc[inf_idxs,M] = -1
    #missing_idxs = results[M] == -1
    #results.loc[missing_idxs,M] = 0

  return (results, methods)

def get_method_names(results):
  methnames = [K for K in results.keys() if K != 'runid' and K not in SIM_PARAMS]
  return methnames

def plot_distplot(results, methods):
  hist = [results[M] for M in methods]
  fig = ff.create_distplot(
    hist,
    methods,
    bin_size = 0.05,
  )
  return fig

def plot_violin(results, methods):
  points = [(M, V) for M in methods for V in results[M]]
  X, Y = zip(*points)
  traces = [{
    'type': 'violin',
    'x': X,
    'y': Y,
    'line': { 'color': 'black' },
    'meanline': { 'visible': True },
    'fillcolor': '#8dd3c7',
    'opacity': 0.6,
    'x0': 'lol'
  }]

  fig = {
    'data': traces,
    'layout': {
      'title': 'blah',
      'yaxis': { 'zeroline': False },
    }
  }
  return fig

def plot_bar(results, methods, key='S'):
  methods = sorted(methods)

  traces = []
  for V in sorted(pd.unique(results[key])):
    Y = []
    valid_methods = []
    hovertext = []
    for M in methods:
      if M.endswith('uniform'):
        continue
      R = results[results[key] == V][M]
      failed = R == -1
      R_succ = R[np.logical_not(failed)]
      R_succ -= results[results[key] == V]['truth'][np.logical_not(failed)]
      if len(R_succ) == 0:
        continue
      Y.append(np.mean(R_succ))
      #Y.append(len(R_succ) / len(R))
      valid_methods.append(M)
      hovertext.append('complete=%.3f' % (len(R_succ) / len(R)))

    traces.append({
      'type': 'bar',
      'x': valid_methods,
      'y': Y,
      'name': 'S=%s' % V,
      'text': hovertext
    })

  fig = {
    'data': traces,
    'layout': {
      'title': 'blah',
      'barmode': 'group',
    }
  }
  return fig

def write_fig(fig, outfn):
  plot = plotly.offline.plot(
    fig,
    output_type = 'div',
    include_plotlyjs = 'cdn',
  )
  with open(outfn, 'w') as outf:
    print(plot, file=outf)

def print_stats(results, methods):
  #for key in ('S', 'T'):
  for key in ('S',):
    uniq_vals = pd.unique(results[key])
    for M in methods:
      for V in sorted(uniq_vals):
        if M.endswith('uniform'):
          continue

        R = results[results[key] == V][M]
        failed = R == -1
        R_succ = R[np.logical_not(failed)]
        if len(R_succ) == 0:
          continue

        print(
          key,
          V,
          M,
          len(R),
          '%.3f' % (np.sum(failed) / len(R)),
          '%.3f' % np.mean(R_succ),
          '%.3f' % np.median(R_succ),
        )
        print('')
        #print(results[key == V])

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('results_fn')
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  results, methods = load_results(args.results_fn)
  #fig = plot_distplot(results, methods)
  #fig = plot_violin(results, methods)
  fig = plot_bar(results, methods)
  write_fig(fig, args.plot_fn)

  #fig = print_stats(results, methods)

if __name__ == '__main__':
  main()
