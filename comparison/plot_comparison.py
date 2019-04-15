import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import re
import pandas as pd
from collections import defaultdict
import numpy as np
import scipy.stats
import os
import sys

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

external_css = [{
  'href': 'https://stackpath.bootstrapcdn.com/bootstrap/4.2.1/css/bootstrap.min.css',
  'rel': 'stylesheet',
  'integrity': 'sha384-GJzZqFGwb1QTTN6wy59ffF1BuGJpLSa9DkKMp0DgiMDm4iYMj70gZWKYbI706tWS',
  'crossorigin': 'anonymous'
}]
app = dash.Dash(__name__, external_stylesheets=external_css)

def augment(results):
  if len(SIM_PARAMS) < 1:
    return results
  simparams = defaultdict(list)
  sim_param_names = ''.join(SIM_PARAMS)
  for rid in results.index:
    for token in re.findall(rf'([{sim_param_names}]\d+)', rid):
      K, V = token[0], token[1:]
      simparams[K].append(int(V))

  lengths = np.array([len(V) for V in simparams.values()])
  if not np.all(lengths == len(results)):
    raise Exception('Some SIM_PARAMS missing from some points')

  for K in simparams.keys():
    results[K] = simparams[K]
  return results

def get_method_names(results):
  methnames = [K for K in results.columns if K not in SIM_PARAMS and 'uniform' not in K]
  return methnames

def init_content(result_types):
  controls = []

  controls.append(html.Label(children='Result type'))
  controls.append(dcc.Dropdown(
    id = 'result_type',
    options = [{'label': rt, 'value': rt} for rt in result_types],
    value = list(result_types)[0],
  ))

  for axis in ('x', 'y'):
    controls.append(html.Label(children=f'Method {axis.upper()}'))
    controls.append(dcc.Dropdown(id = 'method_%s' % axis))

  for param_name in SIM_PARAMS:
    controls.append(html.Label(children=SIM_PARAM_LABELS[param_name]))
    controls.append(dcc.Dropdown(id = f'{param_name}_chooser', multi=True))

  controls.append(html.Label(children='Jitter'))
  # Note that `values` must be specified, or I get a weird exception concerning
  # `TypeError: Object of type method is not JSON serializable`.
  controls.append(dcc.Slider(
    id = 'jitter',
    min = 0,
    max = 0.1,
    step = 0.02,
    value = 0.05
  ))

  controls.append(dcc.Checklist(
    id = 'toggleable_options',
    options = [
      {'label': 'Show missing data', 'value': 'show_missing'},
      {'label': 'Make Quaid happy', 'value': 'make_quaid_happy'},
    ],
    values = [],
  ))

  plot = dcc.Graph(
    id = f'results_plot',
    config = {
      'showLink': True,
      'toImageButtonOptions': {
        'format': 'svg',
        'width': 750,
        'height': 450,
      },
    }
  )
  children = [
    html.Div(className='row', children=html.H1(children='Mutation relations')),
    html.Div(className='row', children=[
      html.Div(className='col-sm-8', children=plot),
      html.Div(className='col-sm-4', children=controls),
    ])
  ]

  return children

def jitter_points(P, magnitude=0.01):
  if len(P) == 0:
    return P
  # If all points in `P` have the same value, avoid making the jitter 0.
  scale = np.max((0.01, magnitude * (np.max(P) - np.min(P))))
  noise = scale * np.random.uniform(low=-1, high=1, size=len(P))
  return P + noise
      
def update_plot(plot_type, methx, methy, results, jitter, toggleable_options, *simparams):
  if methx is None or methy is None or None in simparams:
    return
  show_missing = 'show_missing' in toggleable_options
  make_quaid_happy = 'make_quaid_happy' in toggleable_options

  filters = [np.array(len(results)*[True])]
  for param_name, param_vals in zip(SIM_PARAMS, simparams):
    filters.append(np.array(results[param_name].isin(param_vals)))
  visible = results[np.logical_and.reduce(filters)]
  visible = visible[np.logical_not(np.logical_and(np.isnan(visible[methx]), np.isnan(visible[methy])))]

  total = len(visible)
  present_X = total - np.sum(np.isnan(visible[methx]))
  present_Y = total - np.sum(np.isnan(visible[methy]))
  if not show_missing:
    visible = visible[np.logical_not(np.logical_or(np.isnan(visible[methx]), np.isnan(visible[methy])))]
  X = np.copy(visible[methx])
  Y = np.copy(visible[methy])
  L = np.copy(visible.index)
  assert np.logical_not(np.any(np.logical_or(np.isinf(X), np.isinf(Y))))

  if len(X) > 0 and len(Y) > 0:
    threshold = 0.01
    Y_gt_X = np.sum((Y - X) > threshold) / float(len(X))
    X_gt_Y = np.sum((X - Y) > threshold) / float(len(X))
    Y_approx_X = np.sum(np.abs(Y - X) <= threshold) / float(len(X))
  else:
    Y_gt_X = 0
    X_gt_Y = 0
    Y_approx_X = 0

  X_orig = X
  Y_orig = Y
  if jitter > 0:
    X = jitter_points(X, jitter)
    Y = jitter_points(Y, jitter)

  if len(X) + len(Y) > 0:
    max_val = np.max(np.concatenate((X, Y)))
  else:
    max_val = 1

  marker = { 'size': 14 }
  if len(X) == 0:
    diag_topright = 1.
    diag_bottomleft = 0.
  else:
    diag_topright = np.min([np.max(X), np.max(Y)])
    diag_bottomleft = np.max((np.min(X), np.min(Y)))

  if total > 0:
    present_X_ratio = present_X / total
    present_Y_ratio = present_Y / total
  else:
    present_X_ratio = 0
    present_Y_ratio = 0

  if make_quaid_happy:
    axis_label = {
      'mutrel': 'Average pairwise score for %s',
      'mutphi': 'Error of %s (bits)',
      'merged': 'Something',
    }[plot_type]
    xtitle = axis_label % HAPPY_METHOD_NAMES.get(methx, methx)
    ytitle = axis_label % HAPPY_METHOD_NAMES.get(methy, methy)
    title = None
    template = 'plotly_white'
  else:
    xtitle = methx
    ytitle = methy
    title = '%s %s vs. %s (total=%s, X=%.3f, Y=%.3f, <br>Y_approx_X=%.3f, Y_gt_X=%.3f, X_gt_Y=%.3f)' % (
      plot_type,
      methy,
      methx,
      len(visible),
      present_X / total,
      present_Y / total,
      Y_approx_X,
      Y_gt_X,
      X_gt_Y,
    )
    template = 'seaborn'

    # SciPy wants input in shape (# of dims, # of data).
    XY = np.vstack((X, Y))
    try:
      density = scipy.stats.gaussian_kde(XY)(XY)
    except (ValueError, np.linalg.LinAlgError):
      density = np.zeros(X.shape)
    marker['color'] = density
    marker['colorscale'] = 'Viridis'

  assert len(L) == len(X_orig) == len(Y_orig)
  labels = [f'{runid}<br>orig=({xorig:.2f}, {yorig:.2f})' for runid, xorig, yorig in zip(L, X_orig, Y_orig)]

  axis_range = None
  plot = {
    'data': [go.Scatter(
      x = X,
      y = Y,
      text = labels,
      mode = 'markers',
      marker = marker,
    )],
    'layout': go.Layout(
      xaxis = {'title': xtitle, 'range': axis_range},
      yaxis = {'title': ytitle, 'range': axis_range},
      template = template,
      hovermode = 'closest',
      height = 800,
      title = title,
      font = {'size': 12},
      shapes = [{
        'type': 'line',
        'xref': 'x',
        'yref': 'y',
        'x0': diag_bottomleft,
        'y0': diag_bottomleft,
        'x1': diag_topright,
        'y1': diag_topright,
        'opacity': 0.3,
        'line': {
          'width': 3,
          'dash': 'dot',
        },
      }]
    )
  }
  return plot

def create_callbacks(results):
  for target in ('method_x', 'method_y'):
    app.callback(
      dash.dependencies.Output(target, 'options'),
      [dash.dependencies.Input('result_type', 'value')]
    )(lambda rt: [{'label': M, 'value': M} for M in get_method_names(results[rt])])
    app.callback(
      dash.dependencies.Output(target, 'value'),
      [dash.dependencies.Input(target, 'options')]
    )(lambda opts: opts[0]['value'])

  def _make_param_val_cb(P):
    return lambda rt: [{'label': V, 'value': V} for V in sorted(results[rt][P].unique())]
  for P in SIM_PARAMS:
    target = f'{P}_chooser'
    app.callback(
      dash.dependencies.Output(target, 'options'),
      [dash.dependencies.Input('result_type', 'value')]
    )(_make_param_val_cb(P))
    app.callback(
      dash.dependencies.Output(target, 'value'),
      [dash.dependencies.Input(target, 'options')],
    )(lambda opts: [opt['value'] for opt in opts])

  fig_inputs = [
    dash.dependencies.Input('result_type', 'value'),
    dash.dependencies.Input('method_x', 'value'),
    dash.dependencies.Input('method_y', 'value'),
    dash.dependencies.Input('jitter', 'value'),
    dash.dependencies.Input('toggleable_options', 'values'),
  ]
  fig_inputs += [dash.dependencies.Input(f'{P}_chooser', 'value') for P in SIM_PARAMS]

  def _update_plot(result_type, methx, methy, jitter, show_missing, *simparams):
    return update_plot(result_type, methx, methy, results[result_type], jitter, show_missing, *simparams)
  app.callback(
    dash.dependencies.Output('results_plot', 'figure'),
    fig_inputs,
  )(_update_plot)

def load_results(*result_types):
  results = {}
  for K in result_types:
    env_name = '%s_RESULTS' % K.upper()
    if env_name in os.environ:
      results[K] = pd.read_csv(os.environ[env_name])
      results[K] = results[K].set_index('runid')
    results[K][results[K] == -1] = np.nan

  if 'mutphi' in results and 'truth' in results['mutphi']:
    for M in get_method_names(results['mutphi']):
      if M == 'truth':
        continue
      results['mutphi'][M] -= results['mutphi']['truth']

  if len(results) > 1:
    to_merge = []
    for rt, R in results.items():
      namemap = {C: '%s_%s' % (C, rt) for C in R.columns}
      R = R.rename(axis='columns', mapper=namemap)
      to_merge.append(R)

    while len(to_merge) > 1:
      other = to_merge.pop()
      to_merge[0] = to_merge[0].join(other=other, how='outer')
    results['merged'] = to_merge[0]

  return results

def run():
  np.set_printoptions(linewidth=400, precision=3, threshold=sys.maxsize, suppress=True)
  np.seterr(all='raise')#divide='raise', invalid='raise')

  results = load_results('mutrel', 'mutphi')
  assert len(results) > 0

  for K in results:
    results[K] = augment(results[K])
    for M in get_method_names(results[K]):
      inf_idxs = np.isinf(results[K][M])
      if np.any(inf_idxs):
        # Sometimes the LLH may be -inf. (Thanks, PASTRI.) Consider this to be a
        # failed run.
        print('%s %s has %s infs' % (K, M, np.sum(inf_idxs)))
        results[K].loc[inf_idxs,M] = np.nan

  content = init_content(results.keys())
  app.title = 'Pairtree simulation results'
  app.layout = html.Div(children=content, className='container-fluid')
  create_callbacks(results)

if 'SIM_PARAMS' in os.environ:
  SIM_PARAMS = tuple(os.environ['SIM_PARAMS'])
else:
  SIM_PARAMS = tuple()

run()

if __name__ == '__main__':
  app.run_server(debug=True, host='0.0.0.0', port=8000)
else:
  # Expose this for Gunicorn
  server = app.server
