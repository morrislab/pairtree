import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
from common import Models
import resultserializer
import mutphi

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('pairtree_results_fn')
  parser.add_argument('ssms_fn')
  parser.add_argument('mutphi_fn')
  args = parser.parse_args()

  results = resultserializer.load(args.pairtree_results_fn)
  clusters = [[]] + list(results['clusters'])

  clusterings = [clusters for idx in range(len(results['phi']))]
  mphi = mutphi.calc_mutphi(results['phi'], results['llh'], clusterings, args.ssms_fn)
  mutphi.write_mutphi(mphi, args.mutphi_fn)

if __name__ == '__main__':
  main()
