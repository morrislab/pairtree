import numpy as np
from common import Mutrel

def _write_mutrel(results, key):
  if key not in results:
    return
  results['%s_rels' % key] = results[key].rels
  results['%s_vids' % key] = results[key].vids
  del results[key]

def _load_mutrel(results, key):
  if not (('%s_vids' % key) in results and ('%s_rels' % key) in results):
    return
  results[key] = Mutrel(vids=results['%s_vids' % key], rels=results['%s_rels' % key])
  del results['%s_vids' % key]
  del results['%s_rels' % key]

def save(results, resultsfn):
  serialized = dict(results)
  for reltype in ('mutrel', 'clustrel'):
    for arrtype in ('posterior', 'evidence'):
      _write_mutrel(serialized, '%s_%s' % (reltype, arrtype))
  np.savez_compressed(resultsfn, **serialized)

def load(resultsfn):
  results = np.load(resultsfn)
  # Convert from NpzFile type.
  results = dict(results)

  for reltype in ('mutrel', 'clustrel'):
    for arrtype in ('posterior', 'evidence'):
      _load_mutrel(results, '%s_%s' % (reltype, arrtype))
  return results
