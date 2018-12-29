import numpy as np
import scipy.stats
import common
import mutrel
from progressbar import progressbar
import phi_fitter
Models = common.Models
debug = common.debug

MIN_FLOAT = np.finfo(np.float).min

def calc_llh_phi(data_mutrel, supervars, superclusters, cluster_adj, fit_phis=True):
  phi, eta = phi_fitter.fit_phis(cluster_adj, superclusters, supervars, iterations=100, parallel=0)
  K, S = phi.shape
  alpha, beta = calc_beta_params(supervars)
  assert alpha.shape == beta.shape == (K-1, S)
  assert np.allclose(1, phi[0])
  phi_llh = scipy.stats.beta.logpdf(phi[1:,:], alpha, beta)
  phi_llh = np.sum(phi_llh)

  # I had NaNs creep into my LLH when my alpha and beta params were invalid
  # (i.e., when I had elements of beta that were <= 0).
  assert not np.isnan(phi_llh)
  # Prevent LLH of -inf.
  phi_llh = np.maximum(phi_llh, MIN_FLOAT)
  return phi_llh

def calc_llh_mutrel(data_mutrel, supervars, superclusters, cluster_adj):
  tree_mutrel = mutrel.make_mutrel_tensor_from_cluster_adj(cluster_adj, superclusters)
  mutrel_fit = 1 - np.abs(data_mutrel.rels - tree_mutrel.rels)
  # Prevent log of zero.
  mutrel_fit = np.maximum(common._EPSILON, mutrel_llh)
  mutrel_llh = np.sum(np.log(mutrel_fit))
  return mutrel_llh

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

def init_cluster_adj_random(K):
  # Parents for nodes [1, ..., K-1].
  parents = []
  # Note this isn't truly random, since node i can only choose a parent <i.
  # This prevents cycles.
  for idx in range(1, K):
    parents.append(np.random.randint(0, idx))
  cluster_adj = np.eye(K)
  cluster_adj[parents, range(1,K)] = 1
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
  #debug(adj)
  debug((A,B))
  if B == 0 and anc[B,A]:
    # Don't permit cluster 0 to become non-root node, since it corresponds to
    # normal cell population.
    debug('do nothing')
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
    debug('swapping', A, B)
  else:
    # Move B so it becomes child of A. I don't need to modify the A column.
    adj[:,B] = 0
    adj[A,B] = 1
    debug('moving', B, 'under', A)

  np.fill_diagonal(adj, 1)
  permute_adj.blah.add((A, B, np.array2string(adj)))
  return adj
permute_adj.blah = set()

def calc_beta_params(supervars):
  svids = common.extract_vids(supervars)
  V = np.array([supervars[svid]['var_reads'] for svid in svids])
  R = np.array([supervars[svid]['ref_reads'] for svid in svids])
  omega_v = np.array([supervars[svid]['omega_v'] for svid in svids])
  assert np.all(omega_v == 0.5)

  # Since these are supervars, we can just take 2*V and disregard omega_v, since
  # supervariants are always diploid (i.e., omega_v = 0.5).
  alpha = 2*V + 1
  # Must ensure beta is > 0.
  beta = np.maximum(1, R - V + 1)
  assert np.all(alpha > 0) and np.all(beta > 0)
  return (alpha, beta)

def run_chain(data_mutrel, supervars, superclusters, nsamples, progress_queue=None):
  # Ensure each chain gets a new random state.
  np.random.seed()
  assert nsamples > 0
  K = len(superclusters)

  if progress_queue is not None:
    progress_queue.put(0)
  init_choices = (init_cluster_adj_linear, init_cluster_adj_branching, init_cluster_adj_random)
  # Particularly since clusters may not be ordered by mean VAF, a branching
  # tree in which every node comes off the root is the least biased
  # initialization, as it doesn't require any steps that "undo" bad choices, as
  # in the linear or random (which is partly linear, given that later clusters
  # aren't allowed to be parents of earlier ones) cases.
  init_choices = (init_cluster_adj_branching,)
  init_cluster_adj = init_choices[np.random.choice(len(init_choices))]
  cluster_adj = [init_cluster_adj(K)]
  llh = [calc_llh_phi(data_mutrel, supervars, superclusters, cluster_adj[0])]

  for I in range(1, nsamples):
    if progress_queue is not None:
      progress_queue.put(I)
    old_llh, old_adj = llh[-1], cluster_adj[-1]
    new_adj = permute_adj(old_adj)
    new_llh = calc_llh_phi(data_mutrel, supervars, superclusters, new_adj)
    if new_llh - old_llh >= np.log(np.random.uniform()):
      # Accept.
      cluster_adj.append(new_adj)
      llh.append(new_llh)
      debug(I, llh[-1], 'accept', sep='\t')
    else:
      # Reject.
      cluster_adj.append(old_adj)
      llh.append(old_llh)
      debug(I, llh[-1], 'reject', sep='\t')

  return choose_best_tree(cluster_adj, llh)

def choose_best_tree(adj, llh):
  best_llh = -np.inf
  best_adj = None
  for A, L in zip(adj, llh):
    if L > best_llh:
      best_llh = L
      best_adj = A
  return (best_adj, best_llh)

def sample_trees(data_mutrel, supervars, superclusters, trees_per_chain, nchains, parallel):
  jobs = []
  total = nchains * trees_per_chain

  # Don't use (hard-to-debug) parallelism machinery unless necessary.
  if parallel > 0:
    import concurrent.futures
    import multiprocessing
    manager = multiprocessing.Manager()
    # What is stored in progress_queue doesn't matter. The queue is just used
    # so that child processes can signal when they've sampled a tree, allowing
    # the main process to update the progress bar.
    progress_queue = manager.Queue()
    with progressbar(total=total, desc='Sampling trees', unit='tree', dynamic_ncols=True) as pbar:
      with concurrent.futures.ProcessPoolExecutor(max_workers=parallel) as ex:
        for C in range(nchains):
          jobs.append(ex.submit(run_chain, data_mutrel, supervars, superclusters, trees_per_chain, progress_queue))

        # Exactly `total` items will be added to the queue. Once we've
        # retrieved that many items from the queue, we can assume that our
        # child processes are finished sampling trees.
        for _ in range(total):
          # Block until there's something in the queue for us to retrieve,
          # indicating a child process has sampled a tree.
          progress_queue.get()
          pbar.update()

    results = [J.result() for J in jobs]
  else:
    results = []
    for C in range(nchains):
      results.append(run_chain(data_mutrel, supervars, superclusters, trees_per_chain))

  adj, llh = choose_best_tree(*zip(*results))
  return ([adj], [llh])
