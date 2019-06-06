import argparse
import os
import numpy as np

MISSING = -1

def parse_runtime(runtime_path):
  with open(runtime_path) as F:
    times = F.readline().split()
  if len(times) == 3:
    times = dict(zip(('real', 'user', 'sys'), [float(T) for T in times]))
  else:
    times = {'real': 24*3600, 'user': 24*3600, 'sys': 0}
  times['total'] = times['user'] + times['sys']
  return times

def load_runtimes(runtime_args):
  runtimes = {}

  for runtime_arg in runtime_args:
    name, runtime_path = runtime_arg.split('=', 1)
    assert name not in runtimes
    if not os.path.exists(runtime_path):
      runtimes[name] = None
      continue
    runtimes[name] = parse_runtime(runtime_path)

  return runtimes

def print_runtimes(runtimes, time_type):
  methods = sorted(runtimes.keys())
  key = {'wall': 'real', 'cpu': 'total'}[time_type]
  times = [runtimes[M][key] if runtimes[M] is not None else MISSING for M in methods]
  print(*methods, sep=',')
  print(*times, sep=',')

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--time-type', choices=('wall', 'cpu'), required=True)
  parser.add_argument('runtimes', nargs='+')
  args = parser.parse_args()

  runtimes = load_runtimes(args.runtimes)
  print_runtimes(runtimes, args.time_type)

if __name__ == '__main__':
  main()
