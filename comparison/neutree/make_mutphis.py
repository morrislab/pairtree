import argparse
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
import resultserializer
import mutphi

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('neutree_fn')
  parser.add_argument('ssm_fn')
  parser.add_argument('mutphi_fn')
  args = parser.parse_args()

  neutree = np.load(args.neutree_fn, allow_pickle=True)
  mphi = mutphi.calc_mutphi(neutree['phis'], neutree['logscores'], neutree['clusterings'], args.ssm_fn, neutree['counts'])
  mutphi.write_mutphi(mphi, args.mutphi_fn)

if __name__ == '__main__':
  main()
