import argparse
import numpy as np
import pickle
import inputparser

import simulator

def write_data(data, datafn, ssmfn):
  with open(datafn, 'wb') as outf:
    pickle.dump(data, outf)
  inputparser.write_ssms(data['variants_all'], ssmfn)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--seed', dest='seed', type=int)
  parser.add_argument('-K', dest='K', type=int, default=4, help='Number of clusters')
  parser.add_argument('-S', dest='S', type=int, default=3, help='Number of samples')
  parser.add_argument('-T', dest='T', type=int, default=4000, help='Total reads per mutation')
  parser.add_argument('-M', dest='M', type=int, default=10, help='Number of non-garbage mutations')
  parser.add_argument('-G', dest='G', type=int, default=4, help='Number of garbage mutations')
  parser.add_argument('datafn')
  parser.add_argument('ssmfn')
  args = parser.parse_args()

  if args.seed is None:
    seed = np.random.randint(2**31)
  else:
    seed = args.seed
  np.random.seed(args.seed)

  data = simulator.generate_data(args.K, args.S, args.T, args.M, args.G)
  data['seed'] = seed
  write_data(data, args.datafn, args.ssmfn)

if __name__ == '__main__':
  main()
