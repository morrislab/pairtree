import plotly
import plotly.graph_objs as go
import numpy as np

def plot_eta(eta, sampnames):
  # eta: KxS, K = # of clusters, S = # of samples
  etacum = np.cumsum(eta, axis=0)
  traces = [go.Scatter(
    name = 'Population %s' % cidx,
    x = sampnames,
    y = cumrow,
    text = origrow,
    hoverinfo = 'x+text',
    mode = 'lines',
    fill='tonexty'
  ) for cidx, cumrow, origrow in zip(range(len(etacum)), etacum, eta)]
  traces = list(reversed(traces))
  layout = go.Layout(
    xaxis = {'showgrid': False}
  )
  fig = go.Figure(data=traces, layout=layout)

  plotly.offline.plot(
    fig,
    filename = 'eta.html',
    #image_width = 1368,
    #image_height = 768,
  )
