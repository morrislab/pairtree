import argparse
import numpy as np
import pickle
import pandas as pd
import plotly.express as px
import plotly

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import util
import common

def _calc_depth(adj):
  parents = util.find_parents(adj)
  K = len(adj)
  depth = np.nan + np.zeros(K, dtype=np.int)
  depth[0] = 0
  stack = [0]

  while len(stack) > 0:
    node = stack.pop()
    children = np.flatnonzero(parents == node) + 1
    depth[children] = depth[node] + 1
    stack += children.tolist()

  assert not np.any(np.isnan(depth))
  return depth[1:].astype(np.int)

def _calc_num_desc(adj):
  K = len(adj)
  assert adj.shape == (K,K)
  anc = common.make_ancestral_from_adj(adj)
  C = np.sum(anc, axis=1) - 1
  assert C[0] == K - 1
  return C[1:].astype(np.int)

def _calc_phi_std(truth_fn):
  with open(truth_fn, 'rb') as F:
    truth = pickle.load(F)

  phi_std = np.std(truth['phi'], axis=1)
  depth = _calc_depth(truth['adjm'])
  num_desc = _calc_num_desc(truth['adjm'])
  df = pd.DataFrame({
    'phi_std': phi_std[1:],
    'depth': depth,
    'num_desc': num_desc,
  })
  return df

def _process(truth_fns):
  conjoined = pd.DataFrame()
  for T in truth_fns:
    conjoined = conjoined.append(_calc_phi_std(T))
  return conjoined

def plot(df, plot_fn):
  fig = px.box(df, x='num_desc', y='phi_std')

  K = np.max(df.num_desc)
  N = np.linspace(start=np.min(df.num_desc), stop=np.max(df.num_desc), num=1000)
  var_sum = np.sqrt( (N/K) * (1 - N/K)/(K + 1) )
  lol = pd.DataFrame({'N': N, 'var_sum': var_sum})
  fig.add_scatter(x=lol.N, y=lol.var_sum)

  html = plotly.offline.plot(
    fig,
    output_type = 'div',
    #include_plotlyjs = False,
    config = {
      'showLink': True,
      'toImageButtonOptions': {
        'format': 'png',
        'width': 750,
        'height': 450,
        'filename': os.path.basename(plot_fn).rsplit('.', 1)[0],
      },
    },
  )
  with open(plot_fn, 'w') as F:
    print(html, file=F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('plot_fn')
  parser.add_argument('truth_fns', nargs='+')
  args = parser.parse_args()

  conjoined = _process(args.truth_fns)
  plot(conjoined, args.plot_fn)


if __name__ == '__main__':
  main()
