import tarfile
import io
import os
import numpy as np
import json
import time

import mutrel

class Results:
  # Currently, this is backed by an LZMA-compressed tar archive, where each
  # file within is either JSON-encoded or NumPy-encoded.
  # However, switching the implementation so it's backed by HDF5 instead should
  # be easy if I decide to do so -- the API as I've designed it here should
  # allow this.

  def __init__(self, fn):
    self._fn = fn
    self._to_add = {}

    if self._file_exists():
      with self._open() as F:
        self._names = set([self._resolve_name(fullname) for fullname in F.getnames()])
    else:
      self._names = set()

  def _resolve_name(self, fullname):
    return fullname.rsplit('.', 1)[0]

  def _file_exists(self):
    return os.path.exists(self._fn)

  def _open(self, mode='r'):
    assert mode in ('r', 'w')
    return tarfile.open(self._fn, '%s:xz' % mode, format=tarfile.PAX_FORMAT)

  def has(self, name):
    return name in self._names

  def save(self):
    if self._file_exists():
      with self._open() as F:
        for mem in F.getmembers():
          fullname = mem.name
          name = self._resolve_name(fullname)
          if name in self._to_add:
            continue
          self._to_add[name] = {
            'full_name': fullname,
            'mtime': mem.mtime,
            'bytes': F.extractfile(fullname).read(),
          }

    with self._open('w') as F:
      for name, data in self._to_add.items():
        ti = tarfile.TarInfo(data['full_name'])
        ti.mtime = data['mtime']
        ti.size = len(data['bytes'])
        F.addfile(ti, io.BytesIO(data['bytes']))

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
      'mtime': time.time(),
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
    reader = F.extractfile(full_name)
    if data_type == 'npy':
      # Hack. See https://github.com/numpy/numpy/issues/7989#issuecomment-340921579
      # This hack likely means that two copies of the array will exist in memory temporarily. Bleh.
      bio = io.BytesIO(reader.read())
      bio.seek(0)
      return np.load(bio)
    elif data_type == 'json':
      return json.load(reader)
    else:
      raise Exception('Unknown data type: %s' % data_type)

  def get(self, name):
    return self.get_many((name,))[name]

  def get_many(self, names):
    results = {}

    with self._open() as F:
      present = set(F.getnames())
      for name in names:
        for data_type in ('npy', 'json'):
          full_name = '%s.%s' % (name, data_type)
          if full_name in present:
            results[name] = self._load(full_name, data_type, F)
            break
        else:
          results[name] = None

    return results
