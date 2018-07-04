from __future__ import print_function
from pwgsresults.result_loader import ResultLoader as PwgsResults
import re
import argparse
import numpy as np
import itertools

import plotly
import plotly.graph_objs as go

# Interpreting: cd ~/work/steph/data/pairwise.xeno.nocns && for foo in *tsv; do cat $foo | grep -v absent | tail -n+2 | tr ' ' , | while read blah; do echo $(echo $foo | cut -d. -f1),$blah; done; done  | awk -F',' '{print $0","(($6 >= 0.3 && $5 >= 0) ? "True" : "False") }' | grep True | sort -nk5 -t, > /tmp/blah.txt

DEBUG = False

def lpdist(A, B, p=2):
  A, B = np.array(A), np.array(B)
  dist = np.sum(np.abs(A - B)**p)**(1./p)
  assert dist >= 0
  if DEBUG:
    from IPython import embed
    embed()
  return dist

def compare(s1idx, s2idx, tidx, treesumm):
  pops = treesumm[tidx]['populations']
  s1cp = []
  s2cp = []

  for pop in pops.values():
    s1cp.append(pop['cellular_prevalence'][s1idx])
    s2cp.append(pop['cellular_prevalence'][s2idx])

  dist = lpdist(s1cp, s2cp, p=1)
  norm_dist = dist / len(pops)
  return (dist, norm_dist)

def plot_concordance_cdf(dists, samples, num_samps, outfn):
  pairs = [(D['A'], D['B']) for D in dists]
  pairdists = [D['dist'] for D in dists]
  num_pairs = len(pairs)

  X, Y, L = cdf(pairdists, pairs)
  traces = []
  traces.append(go.Scatter(
    mode='lines',
    x = X,
    y = Y,
    text = ['%s, %s' % (samples[A], samples[B]) for (A, B) in L],
    name = 'Distance (all DX sample pairs)',
  ))

  is_interesting = [samples[B].startswith(samples[A] + ' ') or samples[A].startswith(samples[B] + ' ') for A, B in L]
  if np.sum(is_interesting) == 0:
    return
  X, Y, L = zip(*[(x, y, l) for x, y, l, i in zip(X, Y, L, is_interesting) if i])
  traces.append(go.Scatter(
    mode='markers',
    x = X,
    y = Y,
    text = ['%s, %s' % (samples[A], samples[B]) for (A, B) in L],
    name = 'Distance (same mouse pairs)',
    marker = {'size': 12}
  ))

  scatter(traces, 'Concordance (%s samples, %s pairs)' % (num_samps, num_pairs), 'Distance', 'ECDF(x)', outfn)

def cdf(arr, labels=None):
  arr = np.array(arr)
  sorted_idxs = np.argsort(arr)
  ret = [
    arr[sorted_idxs],
    np.linspace(0, 1, len(arr), endpoint=False),
  ]
  if labels is not None:
    ret.append([labels[idx] for idx in sorted_idxs])
  return tuple(ret)

def scatter(traces, title, xtitle, ytitle, outfn, logx = False, xmin = None, xmax = None):
  xaxis = {
    'title': xtitle,
    'type': logx and 'log' or 'linear',
    'range': [xmin, xmax],
  }

  layout = go.Layout(
    title = title,
    hovermode = 'closest',
    xaxis = xaxis,
    yaxis = {
      'title': ytitle,
    },
  )
  fig = go.Figure(data=traces, layout=layout)
  plotly.offline.plot(fig, filename=outfn)

def calc_concordance(samples, treesumm, tablefn, plotfn):
  #Use handbuilt tree.
  tidx = 0

  allsampidxs = [idx for idx, S in enumerate(samples) if re.search(r'^Diagnosis Xeno', S)]
  distlist = []
  for A, B in itertools.combinations(allsampidxs, 2):
    prefix = 'Diagnosis Xeno 7'
    DEBUG = samples[A].startswith(prefix) and samples[B].startswith(prefix)

    dist, norm_dist = compare(A, B, tidx, treesumm)
    distlist.append({'A': A, 'B': B, 'dist': dist, 'norm_dist': norm_dist})
  distlist.sort(key = lambda E: E['dist'])

  for pos, dist in enumerate(distlist):
    dist['rank'] = (pos + 1.) / len(distlist)
  plot_concordance_cdf(distlist, samples, len(allsampidxs), plotfn)

  dists_indexed = {}
  for dist in distlist:
    A, B = dist['A'], dist['B']
    dists_indexed[(A, B)] = dists_indexed[(B, A)] = dist

  with open(tablefn, 'w') as tablef:
    print('sample1', 'sample2', 'distance', 'norm_dist', 'rank', sep='\t', file=tablef)
    for sampidx, sampname in enumerate(samples):
      if not re.search(r'^Diagnosis Xeno \d+$', sampname):
        continue
      for suffix in ('CNS', 'Spleen'):
        other = '%s %s' % (sampname, suffix)
        if other not in samples:
          print(sampname, other, 'absent', 'absent', 'absent', 'absent', sep='\t', file=tablef)
          continue
        otheridx = samples.index(other)

        D = dists_indexed[(sampidx, otheridx)]
        print(sampname, other, D['dist'], D['norm_dist'], D['rank'], sep='\t', file=tablef)

def scatter(traces, title, xtitle, ytitle, outfn, extra_shapes=None):
  if extra_shapes is None:
    extra_shapes = []

  layout = go.Layout(
    title = title,
    hovermode = 'closest',
    xaxis = {
      'title': xtitle,
    },
    yaxis = {
      'title': ytitle,
    },
    shapes = extra_shapes,
  )
  fig = go.Figure(data=traces, layout=layout)
  plotly.offline.plot(fig, filename=outfn)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('treesumm_fn')
  parser.add_argument('mutlist_fn')
  parser.add_argument('table_fn')
  parser.add_argument('plot_fn')
  args = parser.parse_args()

  loader = PwgsResults(args.treesumm_fn, args.mutlist_fn, None)
  samples = loader.params['samples']
  treesumm = loader.tree_summary

  calc_concordance(samples, treesumm, args.table_fn, args.plot_fn)

main()
