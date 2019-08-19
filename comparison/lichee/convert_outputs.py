import re
import sys
import os
import argparse
import json

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import common
import inputparser
import evalutil

def write_mutrels(adjls, scores, clusters, lichee_garbage, true_garbage, weight_trees_by, mutrel_fn):
  clusters = [[]] + list(clusters)
  clusterings = [clusters for _ in range(len(adjls))]
  assert len(clusterings) == len(adjls) == len(scores)
  adjms = [common.convert_adjlist_to_adjmatrix(adjl) for adjl in adjls]

  mrel = evalutil.calc_mutrel_from_trees(adjms, scores, clusterings, weight_trees_by)
  mrel = evalutil.add_garbage(mrel, lichee_garbage)
  mrel = evalutil.add_garbage(mrel, true_garbage)
  evalutil.save_sorted_mutrel(mrel, mutrel_fn)

def _parse_node(line, snv_indices):
  tokens = line.split('\t')
  nidx, vaf_mask, vaf_present = tokens[:3]
  nidx = int(nidx)
  muts = tokens[3:]

  assert vaf_present.startswith('[') and vaf_present.endswith(']')
  vaf_present = [float(token) for token in vaf_present[1:-1].strip().split()]
  vaf = []

  for sidx in vaf_mask:
    if sidx == '0':
      vaf.append(0)
    elif sidx == '1':
      vaf.append(vaf_present.pop())
    else:
      raise Exception('Unknown sidx %s' % sidx)
  assert len(vaf_present) == 0
  assert len(vaf) == len(vaf_mask)

  cluster = []
  for mut in muts:
    assert mut.startswith('snv')
    snvidx = int(mut[3:])
    cluster.append(snv_indices[snvidx])

  return (nidx, vaf, cluster)

def _parse_tree(F):
  for line in F:
    line = line.strip()
    header = re.search(r'^\*\*\*\*Tree (\d+)\*\*\*\*$', line)
    if header is not None:
      break
    elif line == 'SNV info:':
      return None
  else:
    raise Exception('Premature EOF')

  tidx = int(header.group(1))

  adjl = {}
  for line in F:
    edge = re.search(r'^(\d+) -> (\d+)$', line.strip())
    if not edge:
      break
    A, B = [int(node) for node in edge.group(1, 2)]
    if A not in adjl:
      adjl[A] = []
    adjl[A].append(B)
  else:
    raise Exception('Premature EOF')

  score = line.strip()
  assert score.startswith('Error score: ')
  score = float(score.split(': ', 1)[1])

  return (tidx, adjl, score)

def _parse_nodes(F, snv_indices):
  assert next(F).strip() == 'Nodes:'
  nodes = []
  for line in F:
    line = line.strip()
    if line == '':
      break
    nodes.append(_parse_node(line, snv_indices))
  else:
    raise Exception('Premature EOF')
  return nodes

def _parse_trees(F):
  trees = []
  while True:
    tree = _parse_tree(F)
    if tree is None:
      return trees
    trees.append(tree)

def load_results(results_fn, snv_indices):
  with open(results_fn) as F:
    nodes = _parse_nodes(F, snv_indices)
    trees = _parse_trees(F)

  nidxs, vafs, clusters = zip(*nodes)
  assert list(nidxs) == list(range(1, len(nidxs) + 1))
  tidxs, adjls, scores = zip(*trees)
  assert list(tidxs) == list(range(len(tidxs)))

  return (clusters, adjls, scores)

def find_garbage(lichee_clusters, true_clusters):
  lichee_clustered = set([V for C in lichee_clusters for V in C])
  true_clustered = set([V for C in true_clusters for V in C])
  assert lichee_clustered.issubset(true_clustered)
  garbage = true_clustered - lichee_clustered
  return common.sort_vids(garbage)

def write_params(clusters, garbage, adjls, sampnames, scores, params_fn):
  params = {
    'clusters': clusters,
    'garbage': garbage,
    'structures': adjls,
    'scores': scores,
    'samples': sampnames,
  }
  with open(params_fn, 'w') as F:
    json.dump(params, F)

def parse_snv_indices(snvfn):
  snv_indices = {}
  with open(snvfn) as F:
    header = next(F)
    for idx, line in enumerate(F):
      tokens = line.split()
      # `desc` is the VID.
      chrom, pos, desc = tokens[:3]
      snv_indices[idx + 1] = desc
  return snv_indices

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--weight-trees-by', choices=('llh', 'uniform'), default='uniform')
  parser.add_argument('--mutrel', dest='mutrel_fn')
  parser.add_argument('--structures', dest='struct_fn')
  parser.add_argument('lichee_snv_fn')
  parser.add_argument('trees_fn')
  parser.add_argument('pairtree_params_fn')
  args = parser.parse_args()

  snv_indices = parse_snv_indices(args.lichee_snv_fn)
  clusters, adjls, scores = load_results(args.trees_fn, snv_indices)
  params = inputparser.load_params(args.pairtree_params_fn)
  true_clusters = params['clusters']
  true_garbage = params['garbage']
  lichee_garbage = find_garbage(clusters, true_clusters)

  if args.mutrel_fn is not None:
    write_mutrels(adjls, scores, clusters, lichee_garbage, true_garbage, args.weight_trees_by, args.mutrel_fn)
  if args.struct_fn is not None:
    write_params(
      clusters,
      true_garbage + lichee_garbage,
      adjls,
      params['samples'],
      scores,
      args.struct_fn,
    )

if __name__ == '__main__':
  main()
