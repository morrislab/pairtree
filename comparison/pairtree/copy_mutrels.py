import argparse
import numpy as np

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('src_results_fn')
  parser.add_argument('dest_results_fn')
  args = parser.parse_args()

  desired = ['seed', 'clusters', 'garbage']
  for A in ('posterior', 'evidence'):
    for B in ('rels', 'vids'):
      desired.append('clustrel_%s_%s' % (A, B))

  with open(args.src_results_fn, 'rb') as F:
    src_results = np.load(F, allow_pickle=True)
    dest_results = {K: src_results[K] for K in desired}
  np.savez_compressed(args.dest_results_fn, **dest_results)

if __name__ == '__main__':
  main()
