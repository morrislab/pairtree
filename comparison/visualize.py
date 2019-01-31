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

SIM_PARAMS = ('G', 'K', 'M', 'S', 'T')
SIM_PARAM_LABELS = {
  'G': 'Garbage mutations',
  'K': 'Number of clusters',
  'M': 'Non-garbage mutations',
  'S': 'Samples',
  'T': 'Read depth',
}

external_css = [{
  'href': 'https://stackpath.bootstrapcdn.com/bootstrap/4.2.1/css/bootstrap.min.css',
  'rel': 'stylesheet',
  'integrity': 'sha384-GJzZqFGwb1QTTN6wy59ffF1BuGJpLSa9DkKMp0DgiMDm4iYMj70gZWKYbI706tWS',
  'crossorigin': 'anonymous'
}]
app = dash.Dash(__name__, external_stylesheets=external_css)

def augment(results):
  simparams = defaultdict(list)
  sim_param_names = ''.join(SIM_PARAMS)
  for rid in results.runid:
    for token in re.findall(rf'([{sim_param_names}]\d+)', rid):
      K, V = token[0], token[1:]
      simparams[K].append(int(V))

  lengths = np.array([len(V) for V in simparams.values()])
  assert np.all(lengths == len(results))

  for K in simparams.keys():
    results[K] = simparams[K]
  return results

def get_method_names(results):
  methnames = [K for K in results.keys() if K != 'runid']
  return methnames

def init_content(method_names, sim_param_vals):
  controls = []
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

  plot = [dcc.Graph(id=f'{plot_type}_results') for plot_type in ('mutrel', 'mutphi')]
  children = [
    html.Div(className='row', children=html.H1(children='Mutation relations')),
    html.Div(className='row', children=[
      html.Div(className='col-sm-8', children=plot),
      html.Div(className='col-sm-4', children=controls),
    ])
  ]

  return children

def jitter_points(P, magnitude=0.01):
  scale = magnitude * (np.max(P) - np.min(P))
  noise = scale * np.random.uniform(low=-1, high=1, size=len(P))
  return P + noise
      
def update_plot(plot_type, methx, methy, results, jitter, *simparams):
  filters = []
  for param_name, param_vals in zip(SIM_PARAMS, simparams):
    filters.append(results[param_name].isin(param_vals))
  visible = results[np.logical_and.reduce(filters)]

  visible_X = visible[methx] != -1
  visible_Y = visible[methy] != -1
  visible = visible[np.logical_and(visible_X, visible_Y)]
  X = visible[methx]
  Y = visible[methy]
  if len(X) > 0:
    threshold = 0.01
    Y_gt_X = np.sum((Y - X) > threshold) / float(len(X))
    X_gt_Y = np.sum((X - Y) > threshold) / float(len(X))
    Y_approx_X = np.sum(np.abs(Y - X) <= threshold) / float(len(X))
  else:
    Y_gt_X = 0
    X_gt_Y = 0
    Y_approx_X = 0

  if len(X) > 0:
    diag_topright = min(max(X), max(Y))
  else:
    diag_topright = 1

  XY = np.vstack((X, Y))
  # SciPy wants input in shape (# of dims, # of data).
  try:
    density = scipy.stats.gaussian_kde(XY)(XY)
  except (ValueError, np.linalg.LinAlgError):
    density = np.zeros(X.shape)

  if jitter > 0:
    X = jitter_points(X, jitter)
    Y = jitter_points(Y, jitter)

  plot = {
    'data': [go.Scatter(
      x = X,
      y = Y,
      text = visible['runid'],
      mode = 'markers',
      marker = {
        'color': density,
        'colorscale': 'Viridis',
      },
    )],
    'layout': go.Layout(
      xaxis = {'title': methx, 'range': (-0.1, 1.1)},
      yaxis = {'title': methy, 'range': (-0.1, 1.1)},
      hovermode = 'closest',
      height = 800,
      title = '%s %s vs. %s (X=%s, Y=%s, total=%s,<br>Y_approx_X=%.3f, Y_gt_X=%.3f, X_gt_Y=%.3f)' % (
        plot_type,
        methx,
        methy,
        np.sum(visible_X),
        np.sum(visible_Y),
        len(X),
        Y_approx_X,
        Y_gt_X,
        X_gt_Y,
      ),
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
  mutrel_fig_inputs = [
    dash.dependencies.Input('method_x', 'value'),
    dash.dependencies.Input('method_y', 'value'),
    dash.dependencies.Input('jitter', 'value'),
  ]
  mutrel_fig_inputs += [dash.dependencies.Input(f'{P}_chooser', 'value') for P in SIM_PARAMS]

  app.callback(
    dash.dependencies.Output('mutrel_results', 'figure'),
    mutrel_fig_inputs,
  )(lambda methx, methy, jitter, *simparams: update_plot('mutrel', methx, methy, results['mutrel'], jitter, *simparams))
  app.callback(
    dash.dependencies.Output('mutphi_results', 'figure'),
    mutrel_fig_inputs,
  )(lambda methx, methy, jitter, *simparams: update_plot('mutphi', methx, methy, results['mutphi'], jitter, *simparams))

def run(mutrel_results_fn, mutphi_results_fn):
  results = {
    'mutrel': pd.read_csv(mutrel_results_fn),
    'mutphi': pd.read_csv(mutphi_results_fn),
  }

  assert np.all(results['mutrel'].runid == results['mutphi'].runid)
  assert get_method_names(results['mutrel']) == get_method_names(results['mutphi'])
  method_names = get_method_names(results['mutrel'])

  for K in results.keys():
    results[K] = augment(results[K])
  sim_param_vals = {K: sorted(results['mutrel'][K].unique()) for K in SIM_PARAMS}

  app.title = 'Pairtree simulation results'
  content = init_content(method_names, sim_param_vals)
  app.layout = html.Div(children=content, className='container-fluid')
  create_callbacks(results)

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('mutrel_results_fn')
  parser.add_argument('mutphi_results_fn')
  args = parser.parse_args()

  run(args.mutrel_results_fn, args.mutphi_results_fn)
  app.run_server(debug=True, host='0.0.0.0')
else:
  run(os.environ['MUTREL_RESULTS'], os.environ['MUTPHI_RESULTS'])
  # To run via Gunicorn:
  #   MUTREL_RESULTS=../scratch/results/sims.mutrel.txt MUTPHI_RESULTS=../scratch/results/sims.mutphi.txt gunicorn -w 4 -b 0.0.0.0:4000 visualize:server
  # Expose this for Gunicorn
  server = app.server
