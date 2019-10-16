import argparse
import numpy as np

def load_llh(resultsfn):
  F = np.load(resultsfn)
  if 'llh' not in F:
    return None
  llh = F['llh']
  llh /= -1 * F['phi'][0].size * np.log(2)
  return llh

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('results_fn')
  args = parser.parse_args()

  llh = load_llh(args.results_fn)
  if llh is None:
    return
  print(
    'min=%.3f' % np.min(llh),
    'median=%.3f' % np.median(llh),
    'max=%.3f' % np.max(llh),
  )

if __name__ == '__main__':
  main()
