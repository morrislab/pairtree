CORRECTIONS = {
  'SJBALL022612': {
    's26': 1/2.,
  },
}

def correct_vafs(sampid, variants):
  for varid in variants.keys():
    variants[varid]['vaf_correction'] = CORRECTIONS[sampid][varid] if (sampid in CORRECTIONS and varid in CORRECTIONS[sampid]) else 1.

def has_corrections(sampid):
  return sampid in CORRECTIONS

def corrected_vars(sampid):
  if sampid in CORRECTIONS:
    return list(CORRECTIONS[sampid].keys())
  else:
    return []
