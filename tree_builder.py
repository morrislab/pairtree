from common import Models
import numpy as np
from phi_fitter import fit_all_phis

# Matrices: M mutations, K clusters
#   A: MxK, A[m,k]=1 iff mut m is in cluster k

def extract_mut_info(clusters, cidxs, variants):
  # Renumber variants so their indices are contiguous. They may not be
  # contiguous, e.g., when singleton clusters are removed.
  used_vars = set([V for C in clusters for V in C])
  var_map = {old: new for (new, old) in enumerate(sorted(used_vars))}

  S = len(list(variants.values())[0]['total_reads']) # Number of samples
  M = len(used_vars)
  K = len(cidxs)
  assert set(cidxs) == set(range(K))

  A = np.zeros((M, K))
  ref_reads = np.zeros((M, S))
  var_reads = np.zeros((M, S))

  for cidx, C in zip(cidxs, clusters):
    for V in C:
      # Copy variant so we don't modify original dict.
      variant = variants['s%s' % V]
      varidx = var_map[V]
      A[varidx,cidx] = 1
      ref_reads[varidx] = variant['ref_reads']
      var_reads[varidx] = variant['var_reads']

  return (A, ref_reads, var_reads)

def make_adj(relations):
  N = len(relations)
  assert relations.shape == (N, N)
  assert relations[1,0] == Models.B_A

  adj = np.eye(N)
  adj[0,1] = 1

  for I in range(2, N):
    I_placed = False

    for J in reversed(range(I)):
      if relations[I,J] == Models.B_A:
        if I_placed is False:
          adj[J,I] = 1
          I_placed = True
      elif relations[I,J] == Models.diff_branches:
        pass
      else:
        raise Exception('Unexpected relation for (%s,%s): %s' % (I, J, Models._all[relations[I,J]]))

  return adj

def build_tree(relations, clusters, cidxs, variants):
  adj = make_adj(relations)
  A, ref_reads, var_reads = extract_mut_info(clusters, cidxs, variants)
  phi = fit_all_phis(adj, A, ref_reads, var_reads)


  print(adj, phi)
