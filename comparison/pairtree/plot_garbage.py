import plotly.io as pio
import plotly.figure_factory as ff
import plotly.graph_objects as go
import argparse
import os
import re
import json
import numpy as np

GARB_TYPES = ('acquiredtwice', 'wildtypebackmut', 'missedcna', 'uniform')
GARB_NAMES = {
  'wildtypebackmut': 'Wildtype back mutation',
  'acquiredtwice': 'Homoplasy',
  'uniform': 'Uniform',
  'missedcna': 'Missed CNA',
}
MET_NAMES = {
  'prec': 'Precision',
  'recall': 'Recall',
  'spec': 'Specificity',
  'f1': 'F1',
}

def to_html(fig):
  return pio.to_html(
    fig,
    full_html = False,
    include_plotlyjs = 'cdn',
    include_mathjax = 'cdn',
    config = {
      'showLink': True,
      'toImageButtonOptions': {
        'format': 'svg',
        'width': 750,
        'height': 450,
      },
    },
  )

def _load_results(garbresults_fns):
  results = []

  for fn in garbresults_fns:
    tokens = os.path.basename(fn).split('.')[0].split('_')
    assert tokens[0] == 'sim'
    params = {'garb_type': tokens[1]}

    for T in tokens[2:]:
      M = re.match(r'^([A-Za-z]+)([,0-9]+)$', T)
      assert len(M.groups()) == 2
      K, V = M.groups()
      params[K] = V.replace(',', '.')

    with open(fn) as F:
      J = json.load(F)
    results.append((params, J))

  return results

def _combine(results):
  combined = {}
  conf_types = ('tp', 'tn', 'fp', 'fn')

  param_names = set(results[0][0].keys())
  param_names.discard('run')

  for params, res in results:
    key = frozenset([(K, params[K]) for K in param_names])
    if key not in combined:
      combined[key] = {K: 0 for K in conf_types}
    for K in conf_types:
      combined[key][K] += res[K]

  for key in combined:
    combined[key] |= _calc_metrics(combined[key])

  return combined

def _calc_metrics(confusion):
  # Allow us to proceed with division by zero by converting to NumPy types.
  tp = np.int64(confusion['tp'])
  fp = np.int64(confusion['fp'])
  tn = np.int64(confusion['tn'])
  fn = np.int64(confusion['fn'])

  # Allow division by zero and for 0/0.
  np.seterr(divide='ignore', invalid='ignore')
  mets = {
    'acc': (tp + tn) / (tp + tn + fp + fn),
    # prec = p(y=1 | y_hat=1)
    'prec': tp / (tp + fp),
    # `recall`, `sensitivity`, and `tpr` are the same
    # recall = p(y_hat=1 | y=1)
    'recall': tp / (tp + fn),
    # `specificity` and `tnr` are the same
    # spec = 1 - fpr
    # spec = p(y_hat=0 | y=0)
    'spec': tn / (tn + fp),
    'true_garb': tp + fn,
    'result_garb': tp + fp,
    'total': tp + fp + tn + fn,
  }

  for met in ('acc', 'prec', 'recall'):
    if np.isnan(mets[met]):
      # NaN arises when we get 0/0.
      mets[met] = 1.
  mets['f1'] = 2/(1/mets['recall'] + 1/mets['prec'])
  if np.isinf(mets['f1']):
    mets['f1'] = 0.

  for K in mets.keys():
    # Convert to native Python type from NumPy to permit JSON serialization.
    # These will exist as a mix of native and NumPy types, so I need to allow
    # for either.
    mets[K] = getattr(mets[K], 'tolist', lambda: mets[K])()
  return mets

def _make_heatmap(X, Y, Z, title, xtitle, ytitle):
  fig = {
    'data': go.Heatmap(
      x = X,
      y = Y,
      z = Z,
      type = 'heatmap',
      colorscale = 'Viridis',
      zmin = 0.,
      zmax = 1.,
    ),
    'layout': {
      'xaxis': {'title': xtitle, 'type': 'category'},
      'yaxis': {'title': ytitle, 'type': 'category'},
      'title': title,
    },
  }

  return to_html(fig)

def _write_header():
  html = '<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css">'
  html += '<h1><marquee>I LOVE TREES</marquee></h1>'
  return html

def _partition(results, partition_on):
  param_sets = [dict(ps) for ps in results.keys()]
  parted_vals = set([tuple([ps[K] for K in partition_on]) for ps in param_sets])
  parted_vals = sorted(parted_vals, key = lambda A: [float(e) for e in A])
  return parted_vals

def _plot_2dhist(results, mindelta=None):
  mets = ('prec', 'recall', 'f1', 'spec')
  partition_on = ('K', 'mindelta')
  parted_vals = _partition(results, partition_on)

  html = ''
  for gt in GARB_TYPES:
    html += f'<h1>{GARB_NAMES[gt]}</h1>'
    for pv in parted_vals:
      result_key = dict(zip(partition_on, pv))
      if mindelta is not None and result_key['mindelta'] != mindelta:
        continue
      html += '<h2>' + ', '.join([f'{K}={V}' for K,V in result_key.items()]) + '</h2>'
      html += '<table class="table"><thead><tr>' + ''.join(['<th style="text-align: center" class="w-25">%s</th>' % MET_NAMES[met] for met in mets]) + '</tr></thead><tbody><tr>'

      result_key['garb_type'] = gt
      subresults = {K: V for K,V in results.items() if frozenset(result_key.items()).issubset(K)}
      priors   = sorted([dict(K)['prior'] for K in subresults.keys()],   key = lambda A: float(A))
      maxgarbs = sorted([dict(K)['maxgarb'] for K in subresults.keys()], key = lambda A: float(A))

      for met in mets:
        Z = [ [results[frozenset([
          ('prior', prior),
          ('maxgarb', mg),
        ]) | frozenset(result_key.items())][met] for prior in priors] for mg in maxgarbs ]
        html += '<td>%s</td>' % _make_heatmap(
          priors,
          maxgarbs,
          Z,
          MET_NAMES[met],
          'Garbage prior',
          'Max pairwise garbage',
        )
      html += '</tr></tbody></table>'

  return html

def _plot_bar(results, prior=None, max_garbage=None):
  mets = ('prec', 'recall')
  partition_on = ('K', 'mindelta')
  parted_vals = _partition(results, partition_on)

  html = ''
  for pv in parted_vals:
    result_key = dict(zip(partition_on, pv))
    title = ', '.join([f'{K}={V}' for K,V in result_key.items()])
    html += '<h1>%s</h1>' % title

    if prior is not None:
      result_key['prior'] = prior
    if max_garbage is not None:
      result_key['maxgarb'] = max_garbage

    bars = []
    for met in mets:
      X = []
      Y = []
      labels = []
      for gt in GARB_TYPES:
        key = frozenset((result_key | {'garb_type': gt}).items())
        if key not in results:
          continue
        X.append(GARB_NAMES[gt])
        Y.append(results[key][met])
        labels.append('%.0f%%' % (100*Y[-1]))
      bars.append(go.Bar(
        name = MET_NAMES[met],
        x = X,
        y = Y,
        text = labels,
        textposition = 'inside',
        textfont = {'size': 15},
      ))
    fig = go.Figure(data = bars, layout = go.Layout(
      barmode = 'group',
      title = title,
      yaxis = {
        'tickformat': ',.0%',
        'range': (0,1),
      },
      template = 'seaborn',
    ))
    html += to_html(fig)
  return html

def main():
  parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--plot-fn', required=True)
  parser.add_argument('--plot-bars', action='store_true')
  parser.add_argument('--plot-2dhist', action='store_true')
  parser.add_argument('--filter-prior')
  parser.add_argument('--filter-maxgarb')
  parser.add_argument('--filter-mindelta')
  parser.add_argument('garbresults_fns', nargs='+')
  args = parser.parse_args()

  results = _load_results(args.garbresults_fns)
  combined = _combine(results)
  print(json.dumps({str(K): V for K,V in combined.items()}))

  html = _write_header()
  if args.plot_bars:
    html += _plot_bar(combined, args.filter_prior, args.filter_maxgarb)
  if args.plot_2dhist:
    html += _plot_2dhist(combined, args.filter_mindelta)
  with open(args.plot_fn, 'w') as F:
    print(html, file=F)

if __name__ == '__main__':
  main()
