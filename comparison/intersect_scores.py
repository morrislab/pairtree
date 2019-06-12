import argparse
import pandas as pd
import numpy as np

# Compare based on strings, so that we don't have str->float and then
# float->str conversion that results in loss of precision.
MISSING = '-1'

def parse_file(fn):
  results = pd.read_csv(fn, dtype=str)
  for M in ('mle_unconstrained', 'truth'):
    if M in results:
      results = results.drop([M], axis=1)
  return results

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('score_fn_A')
  parser.add_argument('score_fn_B')
  args = parser.parse_args()

  results_A = parse_file(args.score_fn_A)
  results_B = parse_file(args.score_fn_B)

  # pairtree_tensor isn't included in the mutphi results, but we do want it in
  # the timing results. Assign a dummy value so that the "are all relevant
  # methods included?" checks below pass.
  results_B = results_B.assign(pairtree_tensor = 1.1)

  methods_A = set(results_A) - set(['runid'])
  methods_B = set(results_B) - set(['runid'])
  print(args.score_fn_A, args.score_fn_B, methods_A - methods_B)
  print(args.score_fn_A, args.score_fn_B, methods_B - methods_A)
  assert methods_A == methods_B
  assert list(results_A['runid']) == list(results_B['runid'])
  assert (results_A['runid'] == results_B['runid']).all()
  methods = methods_A

  for M in methods:
    complete_A = np.asarray(results_A[M], dtype=np.float) == float(MISSING)
    complete_B = np.asarray(results_B[M], dtype=np.float) == float(MISSING)
    mask = np.logical_and(complete_A, np.logical_not(complete_B))
    removing = np.sum(mask)
    if removing > 0:
      print('Removing %s values for %s from %s' % (removing, M, args.score_fn_A))
      results_A.loc[mask,M] = MISSING
  results_A.to_csv(args.score_fn_A, index=False)

  #print('method', 'total', 'A', 'B', 'A&B', 'A&¬B', '¬A&B')
  #for M in methods:
  #  complete_A = np.asarray(results_A[M], dtype=np.float) == float(MISSING)
  #  complete_B = np.asarray(results_B[M], dtype=np.float) == float(MISSING)
  #  print(
  #    M,
  #    len(results_A),
  #    np.sum(complete_A),
  #    np.sum(complete_B),
  #    np.sum(np.logical_and(complete_A, complete_B)),
  #    np.sum(np.logical_and(complete_A, np.logical_not(complete_B))),
  #    np.sum(np.logical_and(np.logical_not(complete_A), complete_B)),
  #  )

if __name__ == '__main__':
  main()
