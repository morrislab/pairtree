import argparse
import json
import common
import csv

def convert_clusters(scresultsfn, handbuiltfn, varid_map):
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

def write_results(clusters, garbage, handbuilt_outfn):
  J = {
    'clusters': clusters,
    'garbage': garbage,
    'structure': make_structure(clusters),
    'colourings': [{'left': 'D', 'right': 'R1'}],
    'samporders': [],
  }
  tree_types = ('xeno', 'patient')
  J = {'handbuilt.%s' % tt: J for tt in tree_types}
  with open(handbuiltfn, 'w') as F:
    json.dump(J, F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('sampid')
  parser.add_argument('ssm_fn')
  parser.add_argument('scresults_fn')
  parser.add_argument('handbuilt_existing_fn')
  parser.add_argument('handbuilt_out_fn')
  args = parser.parse_args()

  variants = common.parse_ssms(args.sampid, args.ssm_fn)
  varid_map = build_variant_to_varid_map(variants)
  clusters, garbage = convert_clusters(args.scresults_fn, args.handbuilt_out_fn, varid_map)
  write_results(clusters, garbage, args.handbuilt_fn)

main()
