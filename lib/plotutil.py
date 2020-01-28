import os

def read_file(fn):
  basedir = os.path.join(os.path.dirname(__file__), '..', 'plot_resources')
  with open(os.path.join(basedir, fn)) as F:
    return F.read()

def js_on_load(js):
  return '<script type="text/javascript">document.addEventListener("DOMContentLoaded", () => { %s });</script>' % js
