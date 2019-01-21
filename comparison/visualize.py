import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import re
import argparse
import pandas as pd
from collections import defaultdict
import numpy as np
import scipy.stats

SIM_PARAMS = ('G', 'K', 'M', 'S', 'T')

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

def init_content(method_names, results):
  controls = []
  for axis, default_val in (('x', method_names[0]), ('y', method_names[1])):
    controls.append(html.Label(children=f'Method {axis.upper()}'))
    controls.append(dcc.Dropdown(
      id = 'method_%s' % axis,
      options = [{'label': N, 'value': N} for N in method_names],
      value = default_val,
    ))

  for param_name in SIM_PARAMS:
    controls.append(html.Label(children=param_name))
    param_vals = sorted(results[param_name].unique())
    controls.append(dcc.Dropdown(
      id = f'{param_name}_chooser',
      options = [{'label': V, 'value': V} for V in param_vals],
      value = param_vals,
      multi = True,
    ))

  plot = dcc.Graph(id='mutrel_results')
  children = [
    html.Div(className='row', children=html.H1(children='Mutation relations')),
    html.Div(className='row', children=[
      html.Div(className='col-sm-8', children=plot),
      html.Div(className='col-sm-4', children=controls),
    ])
  ]

  return children
      
def update_plot(methx, methy, results, *simparams):
  visible = results
  for param_name, param_vals in zip(SIM_PARAMS, simparams):
    filter = visible[param_name].isin(param_vals)
    visible = visible[filter]

  X = visible[methx]
  Y = visible[methy]
  if len(X) > 0:
    diag_topright = min(max(X), max(Y))
  else:
    diag_topright = 1

  XY = np.vstack((X, Y))
  # SciPy wants input in shape (# of dims, # of data).
  try:
    density = scipy.stats.gaussian_kde(XY)(XY)
  except ValueError:
    density = np.zeros(X.shape)

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
      title = '%s vs. %s (%s runs shown, %s total)' % (methy, methx, len(visible), len(results)),
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
  ]
  mutrel_fig_inputs += [dash.dependencies.Input(f'{P}_chooser', 'value') for P in SIM_PARAMS]
  app.callback(
    dash.dependencies.Output('mutrel_results', 'figure'),
    mutrel_fig_inputs,
  )(lambda methx, methy, *simparams: update_plot(methx, methy, results, *simparams))

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('resultsfn')
  args = parser.parse_args()

  results = pd.read_csv(args.resultsfn)
  method_names = get_method_names(results)
  results = augment(results)

  app.title = 'Pairtree simulation results'
  content = init_content(method_names, results)
  app.layout = html.Div(children=content, className='container-fluid')
  create_callbacks(results)
  app.run_server(debug=True, host='0.0.0.0')

if __name__ == '__main__':
  main()
