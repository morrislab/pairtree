import argparse
import json
import common
import csv

def convert_clusters(scresultsfn, varid_map):
  clusters = []
  garbage = []
  with open(scresultsfn) as F:
    R = csv.DictReader(F, delimiter='\t')
    for row in R:
      chrom, pos, cluster = row['chr'], int(row['st']), row['cluster']
      varid = varid_map['%s_%s' % (chrom, pos)]
      if cluster == 'NA':
        garbage.append(varid)
      else:
        cluster = int(cluster)
        # This will leave first cluster empty, as the handbuilt JSON format requires.
        while len(clusters) < cluster + 1:
          clusters.append([])
        clusters[cluster].append(varid)
  return (clusters, garbage)


def build_variant_to_varid_map(variants):
  return {'%s_%s' % (V['chrom'], V['pos']): int(V['id'][1:]) for V in variants.values()}

def make_structure(clusters):
  N = len(clusters)
  # Make linear tree.
  return {idx: [idx + 1] for idx in range(len(clusters) - 1)}

def add_missing_sex_variants_to_garbage(variants, clusters, garbage):
  # I run SciClone without sex variants, since I don't know how to specify
  # total numbers of locus according to their inputs -- maybe I need to make a
  # quasi-CNA covering all of X and Y in males, but I looked into this and
  # couldn't figure it out. As such, just mark all sex variants as garbage.
  existing = set([V for C in clusters for V in C] + list(garbage))
  vids = sorted([int(V[1:]) for V in variants.keys()])
  for vid, var in variants.items():
    vid = int(vid[1:])
    if vid in existing:
      continue
    assert var['chrom'] in ('X', 'Y')
    garbage.append(vid)

def write_results(clusters, garbage, tree_type, handbuilt_outfn):
  J = {
    'clusters': clusters,
    'garbage': garbage,
    'structure': make_structure(clusters),
    'colourings': [{'left': 'D', 'right': 'R1'}],
    'samporders': [],
  }
  J = {'handbuilt.%s' % tree_type: J}
  with open(handbuilt_outfn, 'w') as F:
    json.dump(J, F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sampid')
  parser.add_argument('tree_type')
  parser.add_argument('ssm_fn')
  parser.add_argument('scresults_fn')
  parser.add_argument('handbuilt_out_fn')
  args = parser.parse_args()

  variants = common.parse_ssms(args.sampid, args.ssm_fn)
  varid_map = build_variant_to_varid_map(variants)
  clusters, garbage = convert_clusters(args.scresults_fn, varid_map)
  add_missing_sex_variants_to_garbage(variants, clusters, garbage)
  write_results(clusters, garbage, args.tree_type, args.handbuilt_out_fn)

main()
