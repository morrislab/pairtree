import argparse
import json
import numpy as np

def collapse(clusters, to_collapse):
  parsed = []
  for group in to_collapse:
    clusts = [int(C) - 1 for C in group.split(',')]
    assert len(clusts) > 1
    assert np.all(1 <= np.array(clusts)) and np.all(np.array(clusts) < len(clusters))
    parsed.append(set(clusts))

  for I in range(len(parsed)):
    for J in range(I+1, len(parsed)):
      assert len(parsed[I] & parsed[J]) == 0

  mapped = {idx: clust for idx, clust in enumerate(clusters)}
  for group in parsed:
    group = sorted(group)
    target = group[0]
    for other in group[1:]:
      mapped[target] += mapped[other]
      del mapped[other]

  collapsed = [mapped[idx] for idx in sorted(mapped.keys())]
  return collapsed

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('in_params_fn')
  parser.add_argument('out_params_fn')
  parser.add_argument('to_collapse', nargs='+')
  args = parser.parse_args()

  with open(args.in_params_fn) as F:
    params = json.load(F)

  collapsed = collapse(params['clusters'], args.to_collapse)
  params['clusters'] = collapsed

  with open(args.out_params_fn, 'w') as F:
    json.dump(params, F)

if __name__ == '__main__':
  main()
