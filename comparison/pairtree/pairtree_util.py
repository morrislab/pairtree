import numpy as np

def choose_subset(results, subset_size):
  keys = ('adjm', 'llh', 'phi')
  lengths = set([len(results[K]) for K in keys])
  assert len(lengths) == 1
  total = lengths.pop()
  assert 0 < subset_size <= total
  assert total % subset_size == 0

  partitions = int(total / subset_size)
  part_idx = np.random.randint(partitions)
  start  = part_idx * subset_size
  end = start + subset_size
  for K in keys:
    results[K] = results[K][start:end]
