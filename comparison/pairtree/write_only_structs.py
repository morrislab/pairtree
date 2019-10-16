import argparse
import json
import numpy as np

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('orig_params_fn')
  parser.add_argument('pairtree_results_fn')
  parser.add_argument('out_params_fn')
  args = parser.parse_args()

  with open(args.orig_params_fn) as F:
    params = json.load(F)

  results = np.load(args.pairtree_results_fn)
  T = np.argmax(results['prob'])
  #assert T == 0
  assert 'structures' not in params
  params['structures'] = [results['struct'][T].tolist()]


  with open(args.out_params_fn, 'w') as F:
    json.dump(params, F)

if __name__ == '__main__':
  main()
