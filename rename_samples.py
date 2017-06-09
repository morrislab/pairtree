# Rename/hide samples in the params.json input file.
from __future__ import print_function
from collections import defaultdict
import argparse
import sys
import json

def parse_hidden(hidden_fn):
  hidden = defaultdict(list)
  with open(hidden_fn) as F:
    for line in F:
      dataset, sample = line.strip().split(',')
      hidden[dataset].append(sample)
  return hidden

def parse_renamed(renamed_fn):
  renamed = defaultdict(dict)
  with open(renamed_fn) as F:
    for line in F:
      dataset, old_name, new_name = line.strip().split(',')
      renamed[dataset][old_name] = new_name
  return renamed

def munge(dataset, json_fn, hidden, renamed):
  with open(json_fn) as F:
    J = json.load(F)
  samples = J['samples']
  hidden = set(hidden[dataset])
  renamed = renamed[dataset]

  # Steph wants a sample hidden, but she has also renamed it.
  # Copy set so we're not updating it as we're iterating over it.
  for oldname in set(hidden) & set(renamed.keys()):
    if oldname not in samples:
      continue
    newname = renamed[oldname]
    print(json_fn, oldname, newname)
    hidden.remove(oldname)
    hidden.add(newname)

  J['hidden_samples'] = list(hidden)

  for oldname, newname in renamed.items():
    if oldname not in samples:
      continue
    idx = samples.index(oldname)
    samples[idx] = newname

  with open(json_fn, 'w') as F:
    json.dump(J, F)
  
def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('dataset')
  parser.add_argument('hidden_fn')
  parser.add_argument('renamed_fn')
  parser.add_argument('json_fn')
  args = parser.parse_args()

  hidden = parse_hidden(args.hidden_fn)
  renamed = parse_renamed(args.renamed_fn)
  munge(args.dataset, args.json_fn, hidden, renamed)

main()
