# This program is useful for inserting true tree structures into Pairtree
# params, allowing you to produce results using those true structures. This
# lets you, e.g., determine how good the muthpi scores are when given the true
# structure.
#
# This program must be updated to use params['structures'] rather than
# params['structure'].
import pickle
import argparse
import json

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import common

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('truth_fn')
  parser.add_argument('params_fn')
  args = parser.parse_args()

  with open(args.truth_fn, 'rb') as F:
    truth = pickle.load(F)
  with open(args.params_fn, 'r') as F:
    params = json.load(F)
  assert 'structure' not in params

  adjl = common.convert_adj_matrix_to_json_adjlist(truth['adjm'])
  params['structure'] = adjl

  with open(args.params_fn, 'w') as F:
    json.dump(params, F)

if __name__ == '__main__':
  main()
