import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import resultserializer
import mutdist

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('neutree_fn')
  parser.add_argument('baseline_fn')
  parser.add_argument('mutdist_fn')
  args = parser.parse_args()

  neutree = np.load(args.neutree_fn, allow_pickle=True)
  baseline = np.load(args.baseline_fn)

  mdist = mutdist.calc_mutdist(neutree['phis'], neutree['logscores'], neutree['clusterings'], baseline, neutree['counts'])
  mutdist.write_mutdist(mdist, args.mutdist_fn)

if __name__ == '__main__':
  main()
