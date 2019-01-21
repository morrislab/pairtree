import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import re
import argparse
import pandas as pd
from collections import defaultdict
import numpy as np

app = dash.Dash(__name__)

def augment(results):
  extra = defaultdict(list)
  for rid in results.runid:
    for token in re.findall(r'([A-Z]\d+)', rid):
      K, V = token[0], token[1:]
      extra[K].append(int(V))

  lengths = np.array([len(V) for V in extra.values()])
  assert np.all(lengths == len(results))

  for K in extra.keys():
    results[K] = extra[K]
  return results

def get_method_names(results):
  methnames = [K for K in results.keys() if K != 'runid']
  return methnames

def init_content(method_names):
  children = [
    html.H1(children='Mutation relations'),
    dcc.Graph(id='mutrel_results'),
  ]
  for axis, default_val in (('x', method_names[0]), ('y', method_names[1])):
    children.append(html.Label(children=f'Method {axis.upper()}'))
    children.append(dcc.Dropdown(
      id = 'method_%s' % axis,
      options = [{'label': N, 'value': N} for N in method_names],
      value = default_val,
    ))
  return children
      
def update_plot(methx, methy, results):
  plot = {
    'data': [go.Scatter(
      x = results[methx],
      y = results[methy],
      text = results['runid'],
      mode = 'markers',
    )],
    'layout': go.Layout(
      xaxis = {'title': methx, 'range': (0, 1)},
      yaxis = {'title': methy, 'range': (0, 1)},
      hovermode = 'closest',
      height = 800,
      title = '%s vs. %s (%s runs)' % (methy, methx, len(results)),
      shapes = [{
        'type': 'line',
        'xref': 'x',
        'yref': 'y',
        'x0': 0,
        'y0': 0,
        'x1': min(max(results[methx]), max(results[methy])),
        'y1': min(max(results[methx]), max(results[methy])),
        'opacity': 0.3,
        'line': {
          'width': 3,
          'dash': 'dot',
        },
      }]
    )
  }
  return plot

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

  content = init_content(method_names)
  app.title = 'Pairtree simulation results'
  app.layout = html.Div(children=content)
  app.callback(
    dash.dependencies.Output('mutrel_results', 'figure'),
    [
      dash.dependencies.Input('method_x', 'value'),
      dash.dependencies.Input('method_y', 'value'),
    ]
  )(lambda methx, methy: update_plot(methx, methy, results))
  app.run_server(debug=True, host='0.0.0.0')

if __name__ == '__main__':
  main()
