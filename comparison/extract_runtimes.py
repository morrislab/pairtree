import argparse
import csv
import os
import re

MISSING = -1

def extract_runid(cmd):
  match = re.search(r'(sim_[^\.]+)', cmd)
  assert match is not None
  return match.groups()[0]

def extract_runtimes(logfn):
  runtimes = {}

  with open(logfn) as F:
    reader = csv.DictReader(F, delimiter='\t')
    for line in reader:
      runid = extract_runid(line['Command'])
      runtime = float(line['JobRuntime'])
      runtimes[runid] = runtime

  return runtimes

def load_batches(batches):
  results = {}
  for logfn in batches:
    batch = os.path.basename(logfn).split('.')[1]
    results[batch] = extract_runtimes(logfn)
  return results

def combine_batches(results):
  runids = set([runid for batch in results.values() for runid in batch.keys()])
  combined = {}
  for K in runids:
    combined[K] = {}
    for batch_name, batch_results in results.items():
      if K in batch_results:
        combined[K][batch_name] = batch_results[K]
      else:
        combined[K][batch_name] = MISSING
  return combined

def print_runtimes(methods, results):
  print('runid', *methods, sep=',')
  for runid in sorted(results.keys()):
    times = [str(results[runid][M]) for M in methods]
    print(runid, *times, sep=',')

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('batches', nargs='+')
  args = parser.parse_args()

  batches = load_batches(args.batches)
  methods = sorted(batches.keys())
  results = combine_batches(batches)
  print_runtimes(methods, results)

if __name__ == '__main__':
  main()
