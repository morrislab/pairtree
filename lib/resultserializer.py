import zipfile
import io
import os
import numpy as np
import json
import time

import mutrel

class Results:
  # Currently, this is backed by an LZMA-compressed zip archive, where each
  # file within is either JSON-encoded or NumPy-encoded.
  # However, switching the implementation so it's backed by HDF5 instead should
  # be easy if I decide to do so -- the API as I've designed it here should
  # allow this.
  #
  # Initially, I used a tar rather than zip archive. However, tar is in an
  # inferior format insofar as it does not allow seeking to specific files
  # being read -- it alternates file headers with file names, so (potentially)
  # the entire file must be decompressed to locate a specific file the user
  # wants. Zip does not suffer this limitation.
  #
  # Pickle is a strictly unacceptable choice since it allows arbitrary code
  # execution. If I used Pickle (as I initially did), users could not safely
  # exchange results with one another.

  def __init__(self, fn):
    self._fn = fn
    self._to_add = {}
    self._compress_type = zipfile.ZIP_LZMA

    if self._file_exists():
      with self._open() as F:
        self._names = set([self._resolve_name(fullname) for fullname in F.namelist()])
    else:
      self._names = set()

  def _resolve_name(self, fullname):
    return fullname.rsplit('.', 1)[0]

  def _file_exists(self):
    return os.path.exists(self._fn)

  def _open(self, mode='r'):
    assert mode in ('r', 'w')
    return zipfile.ZipFile(self._fn, mode, compression=self._compress_type,)

  def has(self, name):
    return name in self._names

  def save(self):
    if self._file_exists():
      with self._open() as F:
        for zi in F.infolist():
          fullname = zi.filename
          name = self._resolve_name(fullname)
          if name in self._to_add:
            continue
          with F.open(zi) as G:
            self._to_add[name] = {
              'full_name': fullname,
              'bytes': G.read(),
              'timestamp': zi.date_time,
            }

    with self._open('w') as F:
      for name, data in self._to_add.items():
        zi = zipfile.ZipInfo(
          filename=data['full_name'],
          date_time=data['timestamp'],
        )
        F.writestr(zi, data['bytes'], compress_type=self._compress_type)

    self._to_add = {}

  def add(self, name, data):
    if isinstance(data, np.ndarray):
      output = io.BytesIO()
      np.save(output, data)
      output = output.getvalue()
      data_type = 'npy'
    else:
      # I like having a newline on my JSON files, since it makes `cat`ting them on the terminal more pleasant.
      output = (json.dumps(data) + '\n').encode('utf-8')
      data_type = 'json'

    self._to_add[name] = {
      'full_name': '%s.%s' % (name, data_type),
      'timestamp': time.localtime(time.time()),
      'bytes': output,
    }
    self._names.add(name)

  def add_mutrel(self, name, mutrel):
    self.add('%s_vids' % name, mutrel.vids)
    self.add('%s_rels' % name, mutrel.rels)

  def get_mutrel(self, name):
    data = self.get_many(['%s_%s' % (name, T) for T in ('vids', 'rels')])
    return mutrel.Mutrel(vids=data['%s_vids' % name], rels=data['%s_rels' % name])

  def has_mutrel(self, name):
    return self.has('%s_rels' % name) and self.has('%s_vids' % name)

  def _load(self, full_name, data_type, F):
    data = F.read(full_name)
    if data_type == 'npy':
      bio = io.BytesIO(data)
      # `allow_pickle` was only made False by default in Numpy 1.16.3, so force
      # it here for older Numpy versions.
      return np.load(bio, allow_pickle=False)
    elif data_type == 'json':
      return json.loads(data.decode('utf-8'))
    else:
      raise Exception('Unknown data type: %s' % data_type)

  def get(self, name):
    return self.get_many((name,))[name]

  def get_many(self, names):
    results = {}

    with self._open() as F:
      present = set(F.namelist())
      for name in names:
        for data_type in ('npy', 'json'):
          full_name = '%s.%s' % (name, data_type)
          if full_name in present:
            results[name] = self._load(full_name, data_type, F)
            break
        else:
          results[name] = None

    return results
