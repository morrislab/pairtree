import argparse
import numpy as np
import pickle
import json

import inputparser
import simulator

def write_data(data, datafn, paramsfn, ssmfn, should_write_clusters):
  with open(datafn, 'wb') as outf:
    pickle.dump(data, outf)

  with open(paramsfn, 'w') as outf:
    params = {
      'samples': data['sampnames'],
    }
    if should_write_clusters:
      params['clusters'] = data['clusters']
      params['garbage'] = data['vids_garbage']
    json.dump(params, outf)

  inputparser.write_ssms(data['variants'], ssmfn)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--seed', dest='seed', type=int)
  parser.add_argument('--write-clusters', action='store_true')
  parser.add_argument('--tree-type', choices=('monoprimary', 'polyprimary'))
  parser.add_argument('-K', dest='K', type=int, default=4, help='Number of clusters')
  parser.add_argument('-S', dest='S', type=int, default=3, help='Number of samples')
  parser.add_argument('-T', dest='T', type=int, default=4000, help='Total reads per mutation')
  parser.add_argument('-M', dest='M', type=int, default=10, help='Number of non-garbage mutations')
  parser.add_argument('-G', dest='G', type=int, default=4, help='Number of garbage mutations')
  parser.add_argument('datafn')
  parser.add_argument('paramsfn')
  parser.add_argument('ssmfn')
  args = parser.parse_args()

  if args.seed is None:
    seed = np.random.randint(2**31)
  else:
    seed = args.seed
  np.random.seed(args.seed)

  data = simulator.generate_data(args.K, args.S, args.T, args.M, args.G, args.tree_type)
  data['seed'] = seed
  write_data(data, args.datafn, args.paramsfn, args.ssmfn, args.write_clusters)

if __name__ == '__main__':
  main()
