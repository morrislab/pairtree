import numpy as np
import scipy.stats
import common
import mutrel
from progressbar import progressbar
import phi_fitter
Models = common.Models
debug = common.debug

from collections import namedtuple
InnerSample = namedtuple('InnerSample', (
  'adj',
  'anc',
  'depth_frac',
  'llh',
  'logfit',
  'progress',
  'W_subtree',
))


def _calc_llh_phi(phi, V, N, omega_v, epsilon=1e-5):
  K, S = phi.shape
  for arr in V, N, omega_v:
    assert arr.shape == (K-1, S)

  assert np.allclose(1, phi[0])
  P = omega_v * phi[1:]
  P = np.maximum(P, epsilon)
  P = np.minimum(P, 1 - epsilon)

  phi_llh = scipy.stats.binom.logpmf(V, N, P)
  phi_llh = np.sum(phi_llh)
  assert not np.isnan(phi_llh)
  assert not np.isinf(phi_llh)
  return phi_llh

def _calc_llh_mutrel(cluster_adj, data_mutrel, superclusters):
  tree_mutrel = mutrel.make_mutrel_tensor_from_cluster_adj(cluster_adj, superclusters)
  # mutrel_error is symmetric, so every off-diagonal element is counted twice
  # in the sum that becomes the LLH. I could avoid this by making the matrix
  # (upper or lower) triangular, but for `logfit` to be useful, I don't want to
  # do this.
  mutrel_error = np.abs(data_mutrel.rels - tree_mutrel.rels)
  assert np.all(mutrel_error <= 1)

  mutrel_fit = 1 - mutrel_error
  # Prevent log of zero.
  mutrel_fit = np.maximum(common._EPSILON, mutrel_fit)
  logfit = np.log(mutrel_fit)

  logfit = np.sum(logfit, axis=(1,2))
  llh = np.sum(logfit)
  return (llh, logfit)

def _init_cluster_adj_linear(K):
  cluster_adj = np.eye(K, dtype=np.int)
  for k in range(1, K):
    cluster_adj[k-1,k] = 1
  return cluster_adj

def _init_cluster_adj_branching(K):
  cluster_adj = np.eye(K, dtype=np.int)
  # Every node comes off node 0, which will always be the tree root. Note that
  # we don't assume that the first cluster (node 1, cluster 0) is the clonal
  # cluster -- it's not treated differently from any other nodes/clusters.
  cluster_adj[0,:] = 1
  return cluster_adj

def _init_cluster_adj_random(K):
  # Parents for nodes [1, ..., K-1].
  parents = []
  # Note this isn't truly random, since node i can only choose a parent <i.
  # This prevents cycles.
  for idx in range(1, K):
    parents.append(np.random.randint(0, idx))
  cluster_adj = np.eye(K, dtype=np.int)
  cluster_adj[parents, range(1,K)] = 1
  return cluster_adj

def _modify_tree(adj, anc, A, B):
  '''If `B` is ancestral to `A`, swap nodes `A` and `B`. Otherwise, move
  subtree `B` under `A`.

  `B` can't be 0 (i.e., the root node), as we want always to preserve the
  property that node zero is root.'''
  K = len(adj)
  # Ensure `B` is not zero.
  assert 0 <= A < K and 0 < B < K

  if A == B:
    return adj

  adj = np.copy(adj)
  anc = np.copy(anc)

  assert np.array_equal(np.diag(adj), np.ones(K))
  # Diagonal should be 1, and every node except one of them should have a parent.
  assert np.sum(adj) == K + (K - 1)
  # Every column should have two 1s in it corresponding to self & parent,
  # except for column denoting root.
  assert np.array_equal(np.sort(np.sum(adj, axis=0)), np.array([1] + (K - 1)*[2]))

  np.fill_diagonal(adj, 0)
  np.fill_diagonal(anc, 0)

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
    #debug('tree_permute', (A,B), 'swapping', A, B)
  else:
    # Move B so it becomes child of A. I don't need to modify the A column.
    adj[:,B] = 0
    adj[A,B] = 1
    #debug('tree_permute', (A,B), 'moving', B, 'under', A)

  np.fill_diagonal(adj, 1)
  return adj

def _choose_nodes_uniformly(K, anc):
  while True:
    A, B = np.random.choice(K, size=2, replace=False)
    assert A != B
    # Prevent `B` from being the root node.
    # Also prevent `B` from being ancestral to `A`.
    if B != 0 and not anc[B,A]:
      break
  return (A, B)

def _calc_depth_frac(adj):
  K = len(adj)
  root = 0
  Z = np.copy(adj)
  np.fill_diagonal(Z, 0)

  depth = np.zeros(K)
  stack = [root]
  while len(stack) > 0:
    P = stack.pop()
    C = np.flatnonzero(Z[P])
    if len(C) == 0:
      continue
    depth[C] = depth[P] + 1
    stack += list(C)

  assert depth[root] == 0
  assert np.all(depth[:root] > 0) and np.all(depth[root+1:] > 0)

  depth_frac = depth / np.max(depth)
  return depth_frac

def calc_binom_params(supervars):
  svids = common.extract_vids(supervars)
  V = np.array([supervars[svid]['var_reads'] for svid in svids])
  R = np.array([supervars[svid]['ref_reads'] for svid in svids])
  omega_v = np.array([supervars[svid]['omega_v'] for svid in svids])
  assert np.all(omega_v == 0.5)

  N = V + R
  return (V, N, omega_v)

def _find_parents(adj):
  adj = np.copy(adj)
  np.fill_diagonal(adj, 0)
  return np.argmax(adj[:,1:], axis=0)

def _find_parent(node, adj):
  col = np.copy(adj[:,node])
  col[node] = 0
  parents = np.flatnonzero(col)
  assert len(parents) == 1
  return parents[0]

def _make_W_subtree(adj, depth_frac, logfit, progress):
  # Hyperparams:
  # * tau: weight of depth term
  # * rho: weight of mutrel fit term
  # * psi: how strongly peaked depth term is
  tau = 1
  rho = 5
  psi = 3

  K = len(adj)

  assert np.all(logfit <= 0)
  assert adj.shape == (K, K)

  assert 0 <= progress <= 1
  A = psi * progress + 1
  B = psi * (1 - progress) + 1
  depth_frac[1:] = np.minimum(0.99, depth_frac[1:])
  depth_frac[1:] = np.maximum(0.01, depth_frac[1:])
  weights_depth = depth_frac**(A-1) * (1 - depth_frac)**(B-1)
  weights_depth[0] = 0
  weights_depth /= np.sum(weights_depth)

  assert np.all(logfit <= 0)
  if np.all(logfit == 0):
    weights_fit = np.ones(K-1)
  else:
    weights_fit = np.maximum(1e-5, logfit)
  weights_fit /= np.sum(weights_fit)

  weights      = tau * weights_depth
  weights[1:] += rho * weights_fit
  assert weights[0] == 0 and np.all(weights[1:] >= 0) and np.any(weights > 0)

  norm_weights = weights / np.sum(weights)
  _make_W_subtree.debug = (progress, depth_frac, weights_depth, weights_fit, logfit, weights, norm_weights)
  return norm_weights

def _make_W_parents(anc, depth_frac, subtree_idx, mutrel):
  # theta: weight of `B_A` pairwise probabilities
  # kappa: weight of depth_frac
  theta = 8
  kappa = 1

  K = len(anc)
  cluster_idx = subtree_idx - 1
  assert mutrel.vids[cluster_idx] == 'S%s' % cluster_idx

  weights = np.zeros(K)
  weights[1:] += theta * mutrel.rels[cluster_idx,:,Models.B_A]
  weights     += kappa * depth_frac

  # We can use `anc` to put zero mass on choosing a parent that is already a
  # descendant of the node. But, for now, we do not.
  #subtree_desc = np.flatnonzero(anc[subtree_idx])
  #weights[subtree_desc] = 0

  weights /= np.sum(weights)
  # Use normalized weights to set root weight.
  # Intutiion: if none of the non-root nodes are high probability parents, then
  # the root should be high probability.
  weights[0] = max(0.001, 1 - np.max(weights[1:]))
  # Allow all nodes to be chosen with some small probability ...
  weights[1:] = np.maximum(1e-10, weights[1:])
  # ... but don't allow a node to be its own parent.
  weights[subtree_idx] = 0
  # Renormalize.
  weights /= np.sum(weights)

  return weights

def _sample_cat(W):
  return np.random.choice(len(W), p=W)

#def _load_truth(data_mutrel, superclusters):
#  truthfn = '/home/q/qmorris/jawinter/work/pairtree/scratch/inputs/sims.pairtree/sim_K3_S3_T50_M300_G300_run1.data.pickle'
#  import pickle
#  with open(truthfn, 'rb') as F:
#    truth = pickle.load(F)
#  true_adjm = truth['adjm']
#  true_phi = truth['phi']
#  true_parents = _find_parents(true_adjm)
#  true_llh, true_logfit = _calc_llh_mutrel(true_adjm, data_mutrel, superclusters)
#  return (true_parents, true_llh, true_logfit, true_adjm, true_phi)

def _run_inner_mh_chain(init_adj, nsamples, data_mutrel, superclusters):
  init_llh, init_logfit = _calc_llh_mutrel(init_adj, data_mutrel, superclusters)
  init_anc = common.make_ancestral_from_adj(init_adj)
  init_depth_frac = _calc_depth_frac(init_adj)
  init_progress = 0
  samps = [InnerSample(
    adj = init_adj,
    anc = init_anc,
    depth_frac = init_depth_frac,
    llh = init_llh,
    logfit = init_logfit,
    progress = init_progress,
    W_subtree = _make_W_subtree(init_adj, init_depth_frac, init_logfit, init_progress),
  )]

  assert nsamples > 0
  accepted = 0
  for I in range(1, nsamples):
    progress = I / nsamples
    old_samp = samps[-1]
    subtree_idx = _sample_cat(old_samp.W_subtree)
    assert subtree_idx != 0
    # Don't carry `old_W_parents` in `old_samp`, since it's dependent not just
    # on the current tree, but also on the subtree selected.
    old_W_parents = _make_W_parents(old_samp.anc, old_samp.depth_frac, subtree_idx, data_mutrel)
    dbg = tuple(_make_W_subtree.debug)
    new_parent = _sample_cat(old_W_parents)
    new_adj = _modify_tree(old_samp.adj, old_samp.anc, new_parent, subtree_idx)
    new_llh, new_logfit = _calc_llh_mutrel(new_adj, data_mutrel, superclusters)
    new_depth_frac = _calc_depth_frac(new_adj)
    new_samp = InnerSample(
      adj = new_adj,
      anc = common.make_ancestral_from_adj(new_adj),
      depth_frac = new_depth_frac,
      llh = new_llh,
      logfit = new_logfit,
      progress = progress,
      W_subtree = _make_W_subtree(new_adj, new_depth_frac, new_logfit, progress),
    )

    # This is p(subtree_idx, new_parent | old_state) = p(subtree_idx | old_state)p(new_parent | subtree_idx, old_state)
    old_parent = _find_parent(subtree_idx, old_samp.adj)
    log_p_new_given_old = np.log(old_samp.W_subtree[subtree_idx]) + np.log(old_W_parents[new_parent])
    # This is p(subtree_idx, old_parent | new_state) = p(subtree_idx | new_state)p(old_parent | subtree_idx, new_state)
    new_W_parents = _make_W_parents(new_samp.anc, new_samp.depth_frac, subtree_idx, data_mutrel)
    log_p_old_given_new = np.log(new_samp.W_subtree[subtree_idx]) + np.log(new_W_parents[old_parent])

    log_p_transition = (new_samp.llh - old_samp.llh) + (log_p_old_given_new - log_p_new_given_old)
    U = np.random.uniform()

    if log_p_transition >= np.log(U):
      action = 'accept'
      samps.append(new_samp)
      accepted += 1
    else:
      action = 'reject'
      samps.append(old_samp)

    #true_parents, true_llh, true_logfit, true_adjm, true_phi = _load_truth(data_mutrel, superclusters)
    #norm = -1 * (len(old_samp.adj) - 1)**2 * np.log(2)
    #old_mutrel_error = old_samp.llh / norm
    #new_mutrel_error = new_samp.llh / norm
    #true_mutrel_error = true_llh / norm
    #cols = ('progress', 'depth_frac', 'weights_depth', 'weights_fit', 'logfit', 'weights', 'norm_weights', 'nodes', 'old_W_subtree', 'new_W_subtree', 'old_W_parents', 'new_W_parents', 'old_parents', 'new_parents', 'true_parents', 'old_llh', 'new_llh', 'p_old_given_new', 'p_new_given_old', 'p_state_delta', 'old_mutrel_error', 'new_mutrel_error', 'true_mutrel_error', 'action')
    #new_dbg = dbg + ((old_parent, new_parent, subtree_idx), old_samp.W_subtree, new_samp.W_subtree, old_W_parents, new_W_parents, _find_parents(old_samp.adj), _find_parents(new_samp.adj), true_parents, '%.3f' % old_samp.llh, '%.3f' % new_samp.llh, '%.3f' % log_p_old_given_new, '%.3f' % log_p_new_given_old, '%.3f' % (log_p_old_given_new - log_p_new_given_old),  '%.3f' % old_mutrel_error, '%.3f' % new_mutrel_error, '%.3f' % true_mutrel_error, action)
    #debug(*['%s=%s' % (K, V) for K, V in zip(cols, new_dbg)], sep='\t')

  # Before, I called `choose_best_tree` to choose the tree to return from the
  # inner loop. Now, I just select the last tree sampled, on the assumption
  # this is a valid sample from the posterior (i.e., the chain has burned in).
  #
  # Before, when choosing the tree with the highest likelihood, I would often
  # always propose the same tree, missing reasonable structures altogether.
  accept_rate = accepted / (nsamples - 1)
  return (samps[-1].adj, accept_rate)

def _run_inner_mh_chain_simple(init_adj, nsamples, data_mutrel, superclusters):
  init_llh, init_logfit = _calc_llh_mutrel(init_adj, data_mutrel, superclusters)
  adj = [init_adj]
  llh = [init_llh]

  assert nsamples > 0
  accepted = 0
  for I in range(1, nsamples):
    old_adj = adj[-1]
    old_llh = llh[-1]
    old_anc = common.make_ancestral_from_adj(old_adj)
    A, B = _choose_nodes_uniformly(len(old_adj), old_anc)
    new_adj = _modify_tree(old_adj, old_anc, A, B)
    new_llh, new_logfit = _calc_llh_mutrel(new_adj, data_mutrel, superclusters)

    U = np.random.uniform()
    if new_llh - old_llh >= np.log(U):
      action = 'accept'
      accepted += 1
      adj.append(new_adj)
      llh.append(new_llh)
    else:
      action = 'reject'
      adj.append(old_adj)
      llh.append(old_llh)
  accept_rate = accepted / (nsamples - 1)
  return (adj[-1], accept_rate)

def _run_outer_mh_chain(data_mutrel, supervars, superclusters, nsamples, phi_method, phi_iterations, tree_perturbations, seed, progress_queue=None):
  # Ensure each chain gets a new random state. I add chain index to initial
  # random seed to seed a new chain, so I must ensure that the seed is still in
  # the valid range [0, 2**32).
  np.random.seed(seed % 2**32)

  if progress_queue is not None:
    progress_queue.put(0)
  K = len(superclusters)
  V, N, omega_v = calc_binom_params(supervars)

  def _calc_phi(adj):
    phi, eta = phi_fitter.fit_phis(adj, superclusters, supervars, method=phi_method, iterations=phi_iterations, parallel=0)
    return phi
  def _calc_llh(adj, phi):
    return _calc_llh_phi(phi, V, N, omega_v)

  # Particularly since clusters may not be ordered by mean VAF, a branching
  # tree in which every node comes off the root is the least biased
  # initialization, as it doesn't require any steps that "undo" bad choices, as
  # in the linear or random (which is partly linear, given that later clusters
  # aren't allowed to be parents of earlier ones) cases.
  cluster_adj = [_init_cluster_adj_branching(K)]
  phi = [_calc_phi(cluster_adj[0])]
  llh = [_calc_llh(cluster_adj[0], phi[0])]

  assert nsamples > 0
  accepted = 0
  for I in range(1, nsamples):
    if progress_queue is not None:
      progress_queue.put(I)
    old_llh, old_adj, old_phi = llh[-1], cluster_adj[-1], phi[-1]
    new_adj, inner_accept_rate = _run_inner_mh_chain(old_adj, tree_perturbations, data_mutrel, superclusters)
    new_phi = _calc_phi(new_adj)
    new_llh = _calc_llh(new_adj, new_phi)

    U = np.random.uniform()
    if new_llh - old_llh >= np.log(U):
      action = 'accept'
      accepted += 1
      cluster_adj.append(new_adj)
      phi.append(new_phi)
      llh.append(new_llh)
    else:
      action = 'reject'
      cluster_adj.append(old_adj)
      phi.append(old_phi)
      llh.append(old_llh)

    #true_parents, true_llh, true_logfit, true_adjm, true_phi = _load_truth(data_mutrel, superclusters)
    #_, _logfit = _calc_llh_mutrel(new_adj, data_mutrel, superclusters)
    #norm = -old_phi.size * np.log(2)
    #true_mutphi_llh = _calc_llh(true_adjm, true_phi) 
    #cols = ('loop', 'iter', 'action', 'old_llh', 'new_llh', 'llh_delta', 'old_mutphi', 'new_mutphi', 'true_mutphi', 'U', 'old_parents', 'new_parents', 'true_parents', 'inner_accept', 'logfit', 'true_logfit')
    #vals = (
    #  'mh_outer',
    #  I,
    #  action,
    #  '%.3f' % old_llh,
    #  '%.3f' % new_llh,
    #  '%.3f' % (new_llh - old_llh),
    #  '%.3f' % (old_llh / norm),
    #  '%.3f' % (new_llh / norm),
    #  '%.3f' % (true_mutphi_llh / norm),
    #  '%.3f' % U,
    #  _find_parents(old_adj),
    #  _find_parents(new_adj),
    #  true_parents,
    #  '%.3f' % inner_accept_rate,
    #  _logfit,
    #  true_logfit,
    #)
    #debug(*['%s=%s' % (K, V) for K, V in zip(cols, vals)], sep='\t')

  debug('outer_accept_rate=%s' % (accepted/(nsamples - 1)))
  return (cluster_adj, phi, llh)

def use_existing_structures(adjms, supervars, superclusters, phi_method, phi_iterations, parallel=0):
  V, N, omega_v = calc_binom_params(supervars)
  phis = []
  llhs = []

  for adjm in adjms:
    phi, eta = phi_fitter.fit_phis(adjm, superclusters, supervars, method=phi_method, iterations=phi_iterations, parallel=parallel)
    llh = _calc_llh_phi(phi, V, N, omega_v)
    phis.append(phi)
    llhs.append(llh)
  return (np.array(adjms), np.array(phis), np.array(llhs))

def choose_best_tree(adj, llh):
  best_llh = -np.inf
  best_idx = None
  for idx, (A, L) in enumerate(zip(adj, llh)):
    if L > best_llh:
      best_llh = L
      best_idx = idx
  return best_idx

def sample_trees(data_mutrel, supervars, superclusters, trees_per_chain, burnin_per_chain, nchains, phi_method, phi_iterations, tree_perturbations, seed, parallel):
  assert nchains > 0
  jobs = []
  total_per_chain = trees_per_chain + burnin_per_chain
  total = nchains * total_per_chain

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
          # Ensure each chain's random seed is different from the seed used to
          # seed the initial Pairtree invocation, yet nonetheless reproducible.
          jobs.append(ex.submit(_run_outer_mh_chain, data_mutrel, supervars, superclusters, total_per_chain, phi_method, phi_iterations, tree_perturbations, seed + C + 1, progress_queue))

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
      results.append(_run_outer_mh_chain(data_mutrel, supervars, superclusters, total_per_chain, phi_method, phi_iterations, tree_perturbations, seed + C + 1))

  merged_adj = []
  merged_phi = []
  merged_llh = []
  for A, P, L in results:
    merged_adj += A[burnin_per_chain:]
    merged_phi += P[burnin_per_chain:]
    merged_llh += L[burnin_per_chain:]
  return (merged_adj, merged_phi, merged_llh)
