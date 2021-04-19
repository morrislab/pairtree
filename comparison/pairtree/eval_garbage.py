import argparse
import json
import numpy as np

def sort_vids(vids):
  return sorted(vids, key = lambda V: int(V[1:]))

def _load_params(fns):
  params = {}
  for name, fn in fns.items():
    with open(fn) as F:
      params[name] = json.load(F)
  return params

def _parse_params(params):
  truth_nongarb = set([vid for clust in params['truth']['clusters'] for vid in clust])
  truth_garb = set(params['truth']['garbage'])
  result_garb = set(params['result']['garbage'])
  assert len(truth_nongarb & truth_garb) == 0

  all_vids = sort_vids(truth_nongarb | truth_garb)
  truth = np.array([vid in truth_garb for vid in all_vids])
  result = np.array([vid in result_garb for vid in all_vids])
  return (truth, result)

def _calc_metrics(truth, result):
  notresult = np.logical_not(result)
  nottruth = np.logical_not(truth)

  tp = sum(   truth &    result)
  fp = sum(nottruth &    result)
  tn = sum(nottruth & notresult)
  fn = sum(   truth & notresult)
  # Allow division by zero and for 0/0.
  np.seterr(divide='ignore', invalid='ignore')
  mets = {
    'acc': (tp + tn) / (tp + tn + fp + fn),
    # prec = p(y=1 | y_hat=1)
    'prec': tp / (tp + fp),
    # `recall`, `sensitivity`, and `tpr` are the same
    # recall = p(y_hat=1 | y=1)
    'recall': tp / (tp + fn),
    # `specificity` and `tnr` are the same
    # spec = 1 - fpr
    # spec = p(y_hat=0 | y=0)
    'spec': tn / (tn + fp),
    'true_garb': tp + fn,
    'result_garb': tp + fp,
    'muts': tp + fp + tn + fn,
  }

  for met in ('acc', 'prec', 'recall'):
    if np.isnan(mets[met]):
      # NaN arises when we get 0/0.
      mets[met] = 1.
  mets['f1'] = 2/(1/mets['recall'] + 1/mets['prec'])
  if np.isinf(mets['f1']):
    mets['f1'] = 0.

  for K in mets.keys():
    # Convert to native Python type from NumPy to permit JSON serialization.
    # These will exist as a mix of native and NumPy types, so I need to allow
    # for either.
    mets[K] = getattr(mets[K], 'tolist', lambda: mets[K])()

  return mets

def main():
  parser = argparse.ArgumentParser(
    description='HELLO',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--runid')
  parser.add_argument('params_true_fn')
  parser.add_argument('params_garbremoved_fn')
  args = parser.parse_args()

  params = _load_params({'truth': args.params_true_fn, 'result': args.params_garbremoved_fn,})
  truth, result = _parse_params(params)
  mets = _calc_metrics(truth, result)
  print(json.dumps(mets))

if __name__ == '__main__':
  main()
