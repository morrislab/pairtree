import argparse
import pickle

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import common
import inputparser

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('params_fn')
  parser.add_argument('pickle_fn')
  args = parser.parse_args()

  params = inputparser.load_params(args.params_fn)
  adjlist = inputparser.load_structure(params['structure'])
  adjm = common.convert_adjlist_to_adjmatrix(adjlist)
  with open(args.pickle_fn, 'wb') as outf:
    pickle.dump({
      'adjm': adjm,
      'clusters': params['clusters'],
      'vids_good': [V for C in params['clusters'] for V in C],
      'vids_garbage': params['garbage'],
    }, outf)

if __name__ == '__main__':
  main()
