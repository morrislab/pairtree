import common
import numpy as np
np.set_printoptions(linewidth=200)

def make_mut_ancestry_from_cluster_adj(cluster_adj, clusters):
  cluster_anc = common.make_ancestral_from_adj(cluster_adj)
  M = sum([len(clus) for clus in clusters])
  K = len(clusters)
  mut_anc = np.eye(M)

  for k in range(K):
    self_muts = clusters[k]
    desc_clusters = np.flatnonzero(cluster_anc[k])
    desc_muts = [midx for cidx in desc_clusters for midx in clusters[cidx]]
    mut_anc[np.ix_(self_muts, desc_muts)] = 1

  return mut_anc

def make_ancestry_probs(model_probs):
  ancestry_probs = np.copy(model_probs[:,:,common.Models.A_B])
  np.fill_diagonal(ancestry_probs, 1)
  return ancestry_probs

def calc_llh(ancestry_probs, cluster_anc):
  probs = 1 - np.abs(ancestry_probs - cluster_anc)
  # Prevent log of zero.
  probs = np.maximum(1e-20, probs)
  llh = np.sum(np.log(probs))
  return llh

def init_cluster_adj(K):
  cluster_adj = np.eye(K)
  for k in range(1, K):
    cluster_adj[k-1,k] = 1
  return cluster_adj

def permute_adj(adj):
  adj = np.copy(adj)
  K = len(adj)

  assert np.array_equal(np.diag(adj), np.ones(K))
  # Diagonal should be 1, and every node except one of them should have a parent.
  assert np.sum(adj) == K + (K - 1)
  # Every column should have two 1s in it corresponding to self & parent,
  # except for column denoting root.
  assert np.array_equal(np.sort(np.sum(adj, axis=0)), np.array([1] + (K - 1)*[2]))

  anc = common.make_ancestral_from_adj(adj)
  A, B = np.random.choice(K, size=2, replace=False)
  np.fill_diagonal(adj, 0)

  if anc[B,A]:
    adj_BA = adj[B,A]
    assert anc[A,B] == adj[A,B] == 0
    if adj_BA:
      adj[B,A] = 0

    # Swap position in tree of A and B. I need to modify both the A and B
    # columns.
    acol, bcol = np.copy(adj[:,A]), np.copy(adj[:,B])
    arow, brow = np.copy(adj[A,:]), np.copy(adj[B,:])
    adj[A,:], adj[B,:] = brow, arow
    adj[:,A], adj[:,B] = bcol, acol

    if adj_BA:
      adj[A,B] = 1
  else:
    # Move B so it becomes child of A. I don't need to modify the A column.
    adj[:,B] = 0
    adj[:,B][A] = 1

  np.fill_diagonal(adj, 1)
  permute_adj.blah.add((A, B, np.array2string(adj)))
  #print(A, B, np.array2string(adj), sep='\n')
  return adj
permute_adj.blah = set()

def sample_trees(model_probs, clusters):
  ancestry_probs = make_ancestry_probs(model_probs)
  K = len(clusters)

  cluster_adj = [init_cluster_adj(K)]
  cluster_anc = make_mut_ancestry_from_cluster_adj(cluster_adj[0], clusters)
  llh = [calc_llh(ancestry_probs, cluster_anc)]

  for I in range(1000):
    old_llh, old_adj = llh[-1], cluster_adj[-1]
    new_adj = permute_adj(old_adj)
    cluster_anc = make_mut_ancestry_from_cluster_adj(new_adj, clusters)
    new_llh = calc_llh(ancestry_probs, cluster_anc)
    if False or new_llh - old_llh > np.log(np.random.uniform()):
      # Accept.
      cluster_adj.append(new_adj)
      llh.append(new_llh)
      print(I, llh[-1], 'accept', sep='\t')
    else:
      # Reject.
      cluster_adj.append(old_adj)
      llh.append(old_llh)
      print(I, llh[-1], 'reject', sep='\t')

  return (cluster_adj, llh)
