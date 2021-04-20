import plotly.io as pio
import plotly.figure_factory as ff
import plotly.graph_objects as go
import argparse
import os
import re
import json
import numpy as np

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
      M = re.match(r'^([a-z]+)([,0-9]+)$', T)
      assert len(M.groups()) == 2
      K, V = M.groups()
      params[K] = V.replace(',', '.')

    with open(fn) as F:
      J = json.load(F)
    results.append(J | params)

  return results

def _combine(results):
  combined = {}
  conf_types = ('tp', 'tn', 'fp', 'fn')

  for res in results:
    key = tuple([res[K] for K in ('garb_type', 'alpha', 'prior', 'maxgarb')])
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

def _make_ann_heatmap(results, met):
  priors   = set([K[0] for K in results.keys()])
  maxgarbs = set([K[1] for K in results.keys()])
  X = sorted(priors,   key = lambda V: float(V))
  Y = sorted(maxgarbs, key = lambda V: float(V))
  Z = [ [results[(prior,mg)][met] for prior in X] for mg in Y ]
  labels = [ [f'{100*V:.0f}%' for V in row] for row in Z ]

  fig = ff.create_annotated_heatmap(
    Z,
    annotation_text=labels,
    x = ['%s' % x for x in X],
    y = ['%s' % y for y in Y],
    colorscale = 'Viridis',
    zmin = 0.,
    zmax = 1.,
    showscale = True,
    colorbar = {'tickformat': ',.0%'},
  )
  fig.update_xaxes(title = 'Garbage prior', type='category')
  fig.update_yaxes(title = 'Max pairwise garbage', type='category')
  fig.update_layout(title = met)
  return to_html(fig)

def _make_heatmap(results, met):
  priors   = set([K[0] for K in results.keys()])
  maxgarbs = set([K[1] for K in results.keys()])
  X = sorted(priors,   key = lambda V: float(V))
  Y = sorted(maxgarbs, key = lambda V: float(V))
  Z = [ [results[(prior,mg)][met] for prior in X] for mg in Y ]

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
      'xaxis': {'title': 'Garbage prior', 'type': 'category'},
      'yaxis': {'title': 'Max pairwise garbage', 'type': 'category'},
      'title': met,
    },
  }

  return to_html(fig)

def _plot(results):
  mets = ('prec', 'recall', 'f1', 'spec')
  garb_types = set([K[0] for K in results.keys()])
  alphas     = set([K[1] for K in results.keys()])

  html = ''
  html += '<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css">'
  for gt in garb_types:
    html += f'<h1>{gt}</h1>'
    for alpha in alphas:
      html += f'<h2>&alpha;={alpha}</h2>'
      html += '<table class="table"><thead><tr>' + ''.join(['<th style="text-align: center" class="w-25">%s</th>' % met for met in mets]) + '</tr></thead><tbody><tr>'
      for met in mets:
        # Index on 
        subresults = {K[2:]: V for K, V in results.items() if K[:2] == (gt, alpha)}
        html += '<td>%s</td>' % _make_heatmap(subresults, met)
      html += '</tr></tbody></table>'
  return html

def main():
  parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--plot', dest='plot_fn', required=True)
  parser.add_argument('garbresults_fns', nargs='+')
  args = parser.parse_args()

  results = _load_results(args.garbresults_fns)
  combined = _combine(results)
  print(json.dumps({str(K): V for K,V in combined.items()}))

  html = '<h1><marquee>I LOVE TREES</marquee></h1>'
  html += _plot(combined)
  with open(args.plot_fn, 'w') as F:
    print(html, file=F)

if __name__ == '__main__':
  main()
