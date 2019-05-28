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

  with open(args.src_results_fn, 'rb') as F:
    src_results = np.load(F)
    # Must load into memory so the arrays persist after the file is closed.
    src_results = {K: V for (K, V) in src_results.items()}

  desired = ['seed', 'clusters', 'garbage']
  for A in ('posterior', 'evidence'):
    for B in ('rels', 'vids'):
      desired.append('clustrel_%s_%s' % (A, B))

  dest_results = {K: src_results[K] for K in desired}
  np.savez_compressed(args.dest_results_fn, **dest_results)

if __name__ == '__main__':
  main()
