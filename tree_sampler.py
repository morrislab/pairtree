import common
import numpy as np
from tqdm import tqdm
import concurrent.futures
import multiprocessing
Models = common.Models

def make_mutrel_tensor_from_cluster_adj(cluster_adj, clusters):
  cluster_anc = common.make_ancestral_from_adj(cluster_adj)
  # In determining A_B relations, don't want to set mutaitons (i,j), where i
  # and j are in same cluster, to 1.
  np.fill_diagonal(cluster_anc, 0)

  M = sum([len(clus) for clus in clusters])
  K = len(clusters)
  mutrel = np.zeros((M, M, len(Models._all)))

  for k in range(K):
    self_muts = np.array(clusters[k])
    desc_clusters = np.flatnonzero(cluster_anc[k])
    desc_muts = np.array([midx for cidx in desc_clusters for midx in clusters[cidx]])

    if len(self_muts) > 0:
      mutrel[self_muts[:,None,None], self_muts[None,:,None], Models.cocluster] = 1
    if len(self_muts) > 0 and len(desc_muts) > 0:
      mutrel[self_muts[:,None,None], desc_muts[None,:,None], Models.A_B] = 1

  mutrel[:,:,Models.B_A] = mutrel[:,:,Models.A_B].T
  existing = (Models.cocluster, Models.A_B, Models.B_A)
  already_filled = np.sum(mutrel[:,:,existing], axis=2)
  mutrel[already_filled == 0,Models.diff_branches] = 1
  assert np.array_equal(np.ones((M,M)), np.sum(mutrel, axis=2))

  return mutrel

def calc_llh(data_mutrel, tree_mutrel):
  probs = 1 - np.abs(data_mutrel - tree_mutrel)
  # Prevent log of zero.
  probs = np.maximum(1e-20, probs)
  llh = np.sum(np.log(probs))
  return llh

def init_cluster_adj_linear(K):
  cluster_adj = np.eye(K)
  for k in range(1, K):
    cluster_adj[k-1,k] = 1
  return cluster_adj

def init_cluster_adj_branching(K):
  cluster_adj = np.eye(K)
  cluster_adj[0,1] = 1
  cluster_adj[1,range(2,K)] = 1
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
  if B == 0 and anc[B,A]:
    # Don't permit cluster 0 to become non-root node, since it corresponds to
    # normal cell population.
    return adj
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
  return adj
permute_adj.blah = set()

def run_chain(data_mutrel, clusters, nsamples, progress_queue):
  # Ensure each chain gets a new random state.
  np.random.seed()
  assert nsamples > 0
  K = len(clusters)

  init_choices = (init_cluster_adj_linear, init_cluster_adj_branching)
  init_cluster_adj = init_choices[np.random.choice(len(init_choices))]
  cluster_adj = [init_cluster_adj(K)]
  tree_mutrel = make_mutrel_tensor_from_cluster_adj(cluster_adj[0], clusters)
  llh = [calc_llh(data_mutrel, tree_mutrel)]
  progress_queue.put(0)

  for I in range(1, nsamples):
    progress_queue.put(I)
    old_llh, old_adj = llh[-1], cluster_adj[-1]
    new_adj = permute_adj(old_adj)
    tree_mutrel = make_mutrel_tensor_from_cluster_adj(new_adj, clusters)
    new_llh = calc_llh(data_mutrel, tree_mutrel)
    if new_llh - old_llh > np.log(np.random.uniform()):
      # Accept.
      cluster_adj.append(new_adj)
      llh.append(new_llh)
      #print(I, llh[-1], 'accept', sep='\t')
    else:
      # Reject.
      cluster_adj.append(old_adj)
      llh.append(old_llh)
      #print(I, llh[-1], 'reject', sep='\t')

  #print('Last LLH:', llh[-1])
  return (cluster_adj[-1], llh[-1])

def sample_trees(data_mutrel, clusters, nsamples, nchains=1):
  jobs = []
  total = nchains * nsamples

  manager = multiprocessing.Manager()
  progress_queue = manager.Queue()
  with tqdm(total=total, desc='Sampling trees', unit=' trees') as progress_bar:
    with concurrent.futures.ProcessPoolExecutor(max_workers=nchains) as ex:
      for C in range(nchains):
        jobs.append(ex.submit(run_chain, data_mutrel, clusters, nsamples, progress_queue))
      for _ in range(total):
        progress_queue.get()
        progress_bar.update()
  results = [J.result() for J in jobs]

  results.sort(key = lambda result: result[1]) # Sort by LLH
  adj, llh = zip(*results)
  return (adj, llh)
