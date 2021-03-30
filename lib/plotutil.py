import os

def read_file(fn, basedir=None):
  if not basedir:
    if 'PLOTRESOURCES' in os.environ and os.path.exists(os.path.join(os.environ['PLOTRESOURCES'], fn)):
      basedir = os.environ['PLOTRESOURCES']
    else:
      basedir = os.path.join(os.path.dirname(__file__), '..', 'plot_resources')

  with open(os.path.join(basedir, fn)) as F:
    return F.read()

def js_on_load(js):
  return '<script type="text/javascript">document.addEventListener("DOMContentLoaded", () => { %s });</script>' % js

def hide_samples(sampnames, hidden_samples):
  if hidden_samples is None:
    return None
  hidden = set(hidden_samples)
  sampset = set(sampnames)
  assert hidden.issubset(sampset) and len(hidden) < len(sampset)
  visible_sampidxs = [idx for idx, samp in enumerate(sampnames) if samp not in hidden]
  return visible_sampidxs
