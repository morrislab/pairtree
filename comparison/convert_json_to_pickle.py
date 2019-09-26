import argparse
import pickle

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import inputparser
import util

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('params_fn')
  parser.add_argument('pickle_fn')
  args = parser.parse_args()

  params = inputparser.load_params(args.params_fn)
  adjm = util.convert_parents_to_adjmatrix(params['structure'])
  with open(args.pickle_fn, 'wb') as outf:
    pickle.dump({
      'adjm': adjm,
      'clusters': params['clusters'],
      'vids_good': [V for C in params['clusters'] for V in C],
      'vids_garbage': params['garbage'],
    }, outf)

if __name__ == '__main__':
  main()
