import pickle
import numpy as np
import argparse

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('pickle_fn')
  parser.add_argument('results_fn')
  args = parser.parse_args()

  with open(args.pickle_fn, 'rb') as F:
    truth = pickle.load(F)
  np.savez_compressed(
    args.results_fn, 
    clusters = truth['clusters'],
    garbage = truth['vids_garbage'],
    phi = [truth['phi']],
    adjm = [truth['adjm']],
    llh = [0],
  )

if __name__ == '__main__':
  main()
