from common import Models

def fit_phis():
  pass

def build_tree(relations, clusters, cidxs):
  N = len(relations)
  assert relations.shape == (N, N)
  assert len(clusters) == len(cidxs) == N
  assert relations[1,0] == Models.B_A

  adj = {0: [1], 1: []}
  for I in range(2, N):
    adj[I] = []
    I_placed = False

    for J in reversed(range(I)):
      if relations[I,J] == Models.B_A:
        if I_placed is False:
          adj[J].append(I)
          I_placed = True
      elif relations[I,J] == Models.diff_branches:
        pass
      else:
        raise Exception('Unexpected relation for (%s,%s): %s' % (I, J, Models._all[relations[I,J]]))
  return adj

