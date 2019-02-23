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

def get_method_names(results):
  methnames = [K for K in results.keys() if K != 'runid' and K not in SIM_PARAMS]
  return methnames

def init_content(method_names, sim_param_vals, result_types):
  method_names = sorted(method_names)
  controls = []

  controls.append(html.Label(children='Result type'))
  controls.append(dcc.Dropdown(
    id = 'result_type',
    options = [{'label': rt, 'value': rt} for rt in result_types],
    value = result_types[0],
  ))

  for axis, default_val in (('x', method_names[0]), ('y', method_names[1])):
    controls.append(html.Label(children=f'Method {axis.upper()}'))
    controls.append(dcc.Dropdown(
      id = 'method_%s' % axis,
      options = [{'label': N, 'value': N} for N in method_names],
      value = default_val,
    ))

  for param_name in SIM_PARAMS:
    controls.append(html.Label(children=SIM_PARAM_LABELS[param_name]))
    controls.append(dcc.Dropdown(
      id = f'{param_name}_chooser',
      options = [{'label': V, 'value': V} for V in sim_param_vals[param_name]],
      value = sim_param_vals[param_name],
      multi = True,
    ))

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
    values = ['show_missing'],
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
  show_missing = 'show_missing' in toggleable_options
  make_quaid_happy = 'make_quaid_happy' in toggleable_options
  plot_padding = 1.05

  filters = [len(results)*[True]]
  for param_name, param_vals in zip(SIM_PARAMS, simparams):
    filters.append(results[param_name].isin(param_vals))
  visible = results[np.logical_and.reduce(filters)]
  visible = visible[np.logical_not(np.logical_and(visible[methx] == -1, visible[methy] == -1))]

  X = np.copy(visible[methx])
  Y = np.copy(visible[methy])
  assert np.logical_not(np.any(np.logical_or(np.isinf(X), np.isinf(Y))))

  if show_missing:
    X[X == -1] = 0
    Y[Y == -1] = 0
  else:
    missing = np.logical_or(X == -1, Y == -1)
    present = np.logical_not(missing)
    X = X[present]
    Y = Y[present]

  if len(X) > 0 and len(Y) > 0:
    threshold = 0.01
    Y_gt_X = np.sum((Y - X) > threshold) / float(len(X))
    X_gt_Y = np.sum((X - Y) > threshold) / float(len(X))
    Y_approx_X = np.sum(np.abs(Y - X) <= threshold) / float(len(X))
    missing_X_prop = np.sum(visible[methx] == -1) / float(len(visible))
    missing_Y_prop = np.sum(visible[methy] == -1) / float(len(visible))
  else:
    Y_gt_X = 0
    X_gt_Y = 0
    Y_approx_X = 0
    missing_X_prop = 0
    missing_Y_prop = 0


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
  if plot_type == 'mutrel':
    axis_range = (-0.1, 1.1)
  elif plot_type == 'mutphi':
    axis_range = (-0.5, plot_padding * max_val)
  else:
    axis_range = None
  if len(X) == 0:
    diag_topright = 1

  if make_quaid_happy:
    axis_label = {
      'mutrel': 'Average pairwise score for %s',
      'mutphi': 'Error of %s (bits',
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
      len(X),
      1 - missing_X_prop,
      1 - missing_Y_prop,
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

  diag_topright = plot_padding * max_val
  labels = [f'{runid}<br>orig=({xorig:.2f}, {yorig:.2f})' for runid, xorig, yorig in zip(visible['runid'], X_orig, Y_orig)]

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
        'x0': 0,
        'y0': 0,
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
  fig_inputs = [
    dash.dependencies.Input('result_type', 'value'),
    dash.dependencies.Input('method_x', 'value'),
    dash.dependencies.Input('method_y', 'value'),
    dash.dependencies.Input('jitter', 'value'),
    dash.dependencies.Input('toggleable_options', 'values'),
  ]
  fig_inputs += [dash.dependencies.Input(f'{P}_chooser', 'value') for P in SIM_PARAMS]

  app.callback(
    dash.dependencies.Output('results_plot', 'figure'),
    fig_inputs,
  )(lambda result_type, methx, methy, jitter, show_missing, *simparams: update_plot(result_type, methx, methy, results[result_type], jitter, show_missing, *simparams))

def load_results(*result_types):
  results = {}
  for K in result_types:
    env_name = '%s_RESULTS' % K.upper()
    if env_name in os.environ:
      results[K] = pd.read_csv(os.environ[env_name])
  return results

def run():
  np.set_printoptions(linewidth=400, precision=3, threshold=np.nan, suppress=True)
  np.seterr(all='raise')#divide='raise', invalid='raise')

  results = load_results('mutrel', 'mutphi')
  assert len(results) > 0
  result_types = list(results.keys())
  first_result_type = result_types[0]

  all_methods = set([M for K in results for M in get_method_names(results[K])])
  for K in results:
    methods = set(get_method_names(results[K]))
    assert methods.issubset(all_methods)
    for missing in all_methods - methods:
      results[K][missing] = -1

    assert np.all(results[K].runid == results[first_result_type].runid)
    results[K] = augment(results[K])

    for M in all_methods:
      inf_idxs = np.isinf(results[K][M])
      if np.any(inf_idxs):
        # Sometimes the LLH may be -inf. (Thanks, PASTRI.) Consider this to be a
        # failed run.
        print('%s %s has %s infs' % (K, M, np.sum(inf_idxs)))
        results[K].loc[inf_idxs,M] = -1

  app.title = 'Pairtree simulation results'
  sim_param_vals = {K: sorted(results[first_result_type][K].unique()) for K in SIM_PARAMS}
  content = init_content(all_methods, sim_param_vals, result_types)
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