import numpy as np
import scipy.stats
import common
import mutrel
import tqdm
import phi_fitter
import hyperparams as hparams
import math
import util
from numba import njit

Models = common.Models
debug = common.debug
from common import Models, debug, NUM_MODELS
Mutrel = mutrel.Mutrel

from collections import namedtuple
TreeSample = namedtuple('TreeSample', (
  'adj',
  'anc',
  'phi',
  'llh_phi',
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

def _init_cluster_adj_linear(K):
  cluster_adj = np.eye(K, dtype=np.int32)
  for k in range(1, K):
    cluster_adj[k-1,k] = 1
  return cluster_adj

def _init_cluster_adj_branching(K):
  cluster_adj = np.eye(K, dtype=np.int32)
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
  cluster_adj = np.eye(K, dtype=np.int32)
  cluster_adj[parents, range(1,K)] = 1
  return cluster_adj

def _init_cluster_adj_mutrels(data_logmutrel):
  K = len(data_logmutrel.rels) + 1
  adj = np.eye(K, dtype=np.int32)
  in_tree = set((0,))
  remaining = set(range(1, K))

  W_nodes = np.zeros(K)

  while len(remaining) > 0:
    nodeidxs = np.array(sorted(remaining))
    relidxs = nodeidxs - 1
    assert np.all(relidxs >= 0)
    anc_logprobs = data_logmutrel.rels[np.ix_(relidxs, relidxs)][:,:,Models.A_B]
    # These values are currently -inf, but we must replace these so that joint
    # probabilities of being ancestral to remaining don't all become 0.
    np.fill_diagonal(anc_logprobs, 0)

    assert anc_logprobs.shape == (len(remaining), len(remaining))
    log_W_nodes_remaining = np.sum(anc_logprobs, axis=1)
    # Use really "soft" softmax.
    W_nodes[nodeidxs] = _scaled_softmax(log_W_nodes_remaining)

    # Root should never be selected.
    assert W_nodes[0] == 0
    assert np.isclose(1, np.sum(W_nodes))
    assert np.all(W_nodes[list(in_tree)] == 0)
    nidx = _sample_cat(W_nodes)

    log_W_parents = np.full(K, -np.inf)
    others = np.array(sorted(remaining - set((nidx,))))
    truncated_logmutrel = mutrel.remove_variants_by_vidx(data_logmutrel, others - 1)
    for parent in in_tree:
      new_adj = np.copy(adj)
      new_adj[parent,nidx] = 1
      truncated_adj = util.remove_rowcol(new_adj, others)
      tree_logmutrel = _calc_tree_logmutrel(truncated_adj, truncated_logmutrel)
      log_W_parents[parent] = np.sum(np.triu(tree_logmutrel))
    W_parents = _scaled_softmax(log_W_parents)
    assert np.all(W_parents[nodeidxs] == 0)
    pidx = _sample_cat(W_parents)
    adj[pidx,nidx] = 1

    remaining.remove(nidx)
    in_tree.add(nidx)
    W_nodes[nidx] = 0

  assert np.all(W_nodes == 0)
  assert len(in_tree) == K
  assert len(remaining) == 0
  return adj

def _modify_tree(adj, anc, A, B):
  '''If `B` is ancestral to `A`, swap nodes `A` and `B`. Otherwise, move
  subtree `B` under `A`.

  `B` can't be 0 (i.e., the root node), as we want always to preserve the
  property that node zero is root.'''
  K = len(adj)
  # Ensure `B` is not zero.
  assert 0 <= A < K and 0 < B < K
  assert A != B

  adj = np.copy(adj)

  assert np.array_equal(np.diag(adj), np.ones(K))
  # Diagonal should be 1, and every node except one of them should have a parent.
  assert np.sum(adj) == K + (K - 1)
  # Every column should have two 1s in it corresponding to self & parent,
  # except for column denoting root.
  assert np.array_equal(np.sort(np.sum(adj, axis=0)), np.array([1] + (K - 1)*[2]))

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
    #debug('tree_permute', (A,B), 'swapping', A, B)
  else:
    # Move B so it becomes child of A. I don't need to modify the A column.
    adj[:,B] = 0
    adj[A,B] = 1
    #debug('tree_permute', (A,B), 'moving', B, 'under', A)

  np.fill_diagonal(adj, 1)
  return adj

def calc_binom_params(supervars):
  svids = common.extract_vids(supervars)
  V = np.array([supervars[svid]['var_reads'] for svid in svids])
  R = np.array([supervars[svid]['ref_reads'] for svid in svids])
  omega_v = np.array([supervars[svid]['omega_v'] for svid in svids])
  assert np.all(omega_v == 0.5)

  N = V + R
  return (V, N, omega_v)

def _find_parent(node, adj):
  col = np.copy(adj[:,node])
  col[node] = 0
  parents = np.flatnonzero(col)
  assert len(parents) == 1
  return parents[0]

def _scaled_softmax(A, R=100):
  # Ensures `max(softmax(A)) / min(softmax(A)) <= R`.
  #
  # Typically, I use this as a "softer softmax", ensuring that the largest
  # element in the softmax is at most 100x the magnitude of the smallest.
  # Otherwise, given large differences between the minimum and maximum values,
  # the softmax becomes even more sharply peaked, with one element absorbing
  # effectively all mass.
  noninf = np.logical_not(np.isinf(A))
  if np.sum(noninf) == 0:
    return util.softmax(A)
  delta = np.max(A[noninf]) - np.min(A[noninf])
  if np.isclose(0, delta):
    return util.softmax(A)
  B = min(1, np.log(R) / delta)
  return util.softmax(B*A)

def _make_W_nodes_mutrel(adj, anc, data_logmutrel):
  K = len(adj)
  assert adj.shape == (K, K)

  tree_logmutrel = _calc_tree_logmutrel(adj, data_logmutrel)
  pair_error = 1 - np.exp(tree_logmutrel)
  #pair_error *= 1 - anc

  assert np.allclose(0, np.diag(pair_error))
  assert np.allclose(0, pair_error[0])
  assert np.allclose(0, pair_error[:,0])
  pair_error = np.maximum(common._EPSILON, pair_error)
  node_error = np.sum(np.log(pair_error), axis=1)
  if common.debug.DEBUG:
    _make_W_nodes_mutrel.node_error = node_error

  weights = np.zeros(K)
  weights[1:] += _scaled_softmax(node_error[1:])
  weights[1:] += common._EPSILON
  weights /= np.sum(weights)
  assert weights[0] == 0

  return weights

def _make_W_nodes_uniform(adj, anc):
  K = len(adj)
  weights = np.ones(K)
  weights[0] = 0
  weights /= np.sum(weights)
  return weights

def _make_data_logmutrel(mutrel):
  K = len(mutrel.rels)
  valid_models = (Models.A_B, Models.B_A, Models.diff_branches)
  invalid_models = (Models.cocluster, Models.garbage)

  alpha = common._EPSILON
  logrels = np.full(mutrel.rels.shape, np.nan)
  logrels[:,:,invalid_models] = -np.inf
  logrels[:,:,valid_models] = np.log(mutrel.rels[:,:,valid_models] + alpha)

  logrels[range(K),range(K),:] = -np.inf
  logrels[range(K),range(K),Models.cocluster] = 0

  logrels[:,:,valid_models] -= np.log(1 + len(valid_models)*alpha)
  assert np.allclose(0, scipy.special.logsumexp(logrels, axis=2))
  assert not np.any(np.isnan(logrels))

  logmutrel = Mutrel(rels=logrels, vids=mutrel.vids)
  return logmutrel

@njit
def _calc_tree_logmutrel(adj, data_logmutrel):
  node_rels = util.compute_node_relations(adj)
  K = len(node_rels)
  assert node_rels.shape == (K, K)
  assert data_logmutrel.rels.shape == (K-1, K-1, NUM_MODELS)

  clust_rels = node_rels[1:,1:]
  # First row and column of `tree_logmutrel` will always be zero.
  tree_logmutrel = np.zeros((K,K))
  rng = range(K-1)
  for J in rng:
    for K in rng:
      JK_clustrel = node_rels[J+1,K+1]
      tree_logmutrel[J+1,K+1] = data_logmutrel.rels[J,K,JK_clustrel]

  assert np.array_equal(tree_logmutrel, tree_logmutrel.T)
  assert np.all(tree_logmutrel <= 0)
  return tree_logmutrel

def _make_W_dests_mutrel(subtree_head, curr_parent, adj, anc, data_logmutrel):
  assert subtree_head > 0
  assert adj[curr_parent,subtree_head] == 1
  cluster_idx = subtree_head - 1
  assert data_logmutrel.vids[cluster_idx] == 'S%s' % (cluster_idx + 1)
  K = len(adj)

  logweights = np.full(K, -np.inf)
  for dest in range(K):
    if dest == curr_parent:
      continue
    if dest == subtree_head:
      continue
    new_adj = _modify_tree(adj, anc, dest, subtree_head)
    tree_logmutrel = _calc_tree_logmutrel(new_adj, data_logmutrel)
    logweights[dest] = np.sum(np.triu(tree_logmutrel))

  assert not np.any(np.isnan(logweights))
  valid_logweights = np.delete(logweights, (curr_parent, subtree_head))
  assert not np.any(np.isinf(valid_logweights))

  weights = _scaled_softmax(logweights)
  # Since we end up taking logs, this can't be exactly zero. If the logweight
  # is extremely negative, then this would otherwise be exactly zero.
  weights += common._EPSILON
  weights[curr_parent] = 0
  weights[subtree_head] = 0
  weights /= np.sum(weights)
  return weights

def _make_W_dests_uniform(subtree_head, curr_parent, adj, anc):
  K = len(adj)
  weights = np.ones(K)
  weights[subtree_head] = 0
  weights[curr_parent] = 0
  weights /= np.sum(weights)
  return weights

def _sample_cat(W):
  assert np.all(W >= 0) and np.isclose(1, np.sum(W))
  choice = np.random.choice(len(W), p=W)
  assert W[choice] > 0
  return choice

def _load_truth(truthfn):
  if hasattr(common, '_true_adjm') and hasattr(common, '_true_phi'):
    return
  import pickle
  with open(truthfn, 'rb') as F:
    truth = pickle.load(F)
  common._true_adjm = truth['adjm']
  common._true_phi = truth['phi']

def _init_chain(seed, data_logmutrel, __calc_phi, __calc_llh_phi):
  # Ensure each chain gets a new random state. I add chain index to initial
  # random seed to seed a new chain, so I must ensure that the seed is still in
  # the valid range [0, 2**32).
  np.random.seed(seed % 2**32)

  if np.random.uniform() < hparams.iota:
    init_adj = _init_cluster_adj_mutrels(data_logmutrel)
  else:
    # Particularly since clusters may not be ordered by mean VAF, a branching
    # tree in which every node comes off the root is the least biased
    # initialization, as it doesn't require any steps that "undo" bad choices, as
    # in the linear or random (which is partly linear, given that later clusters
    # aren't allowed to be parents of earlier ones) cases.
    K = len(data_logmutrel.rels) + 1
    init_adj = _init_cluster_adj_branching(K)
  common.ensure_valid_tree(init_adj)

  init_anc = util.make_ancestral_from_adj(init_adj)
  init_phi = __calc_phi(init_adj)

  init_samp = TreeSample(
    adj = init_adj,
    anc = init_anc,
    phi = init_phi,
    llh_phi = __calc_llh_phi(init_adj, init_phi),
  )
  return init_samp

def _make_W_nodes_combined(adj, anc, data_logmutrel):
  W_nodes_uniform = _make_W_nodes_uniform(adj, anc)
  W_nodes_mutrel = _make_W_nodes_mutrel(adj, anc, data_logmutrel)
  return np.vstack((W_nodes_uniform, W_nodes_mutrel))

def _make_W_dests_combined(subtree_head, adj, anc, data_logmutrel):
  curr_parent = _find_parent(subtree_head, adj)
  W_dests_uniform = _make_W_dests_uniform(subtree_head, curr_parent, adj, anc)
  W_dests_mutrel = _make_W_dests_mutrel(subtree_head, curr_parent, adj, anc, data_logmutrel)
  return np.vstack((W_dests_uniform, W_dests_mutrel))

def _generate_new_sample(old_samp, data_logmutrel, __calc_phi, __calc_llh_phi):
  K = len(old_samp.adj)
  # When a tree consists of two nodes (i.e., one mutation cluster), proceeding with
  # the normal sample-generating process will produce an error (specifically,
  # when we try to divide by zero in _make_W_dests_uniform). Circumvent this by
  # returning the current (trivial) tree structure.
  if K == 2:
    return (old_samp, 0., 0.)

  # mode == 0: make uniform update
  # mode == 1: make mutrel-informed update
  mode_node_weights = np.array([1 - hparams.gamma, hparams.gamma])
  mode_dest_weights = np.array([1 - hparams.zeta,  hparams.zeta])
  mode_node = _sample_cat(mode_node_weights)
  mode_dest = _sample_cat(mode_dest_weights)

  W_nodes_old = _make_W_nodes_combined(old_samp.adj, old_samp.anc, data_logmutrel)
  B = _sample_cat(W_nodes_old[mode_node])
  W_dests_old = _make_W_dests_combined(
    B,
    old_samp.adj,
    old_samp.anc,
    data_logmutrel,
  )

  A = _sample_cat(W_dests_old[mode_dest])
  #A = _find_parent(B, common._true_adjm)
  new_adj = _modify_tree(old_samp.adj, old_samp.anc, A, B)
  new_phi = __calc_phi(new_adj)
  new_samp = TreeSample(
    adj = new_adj,
    anc = util.make_ancestral_from_adj(new_adj),
    phi = new_phi,
    llh_phi = __calc_llh_phi(new_adj, new_phi),
  )

  # `A_prime` and `B_prime` correspond to the node choices needed to reverse
  # the tree perturbation.
  if old_samp.anc[B,A]:
    # If `B` is ancestral to `A`, the tree perturbation swaps the nodes. Thus,
    # simply invert the swap to reverse the move.
    A_prime = B
    B_prime = A
  else:
    # If `B` isn't ancestral to `A`, the tree perturbation moves the subtree
    # headed by `B` so that `A` becomes its parent. To reverse the move, move
    # the `B` subtree back under its old parent.
    A_prime = _find_parent(B, old_samp.adj)
    B_prime = B

  W_nodes_new = _make_W_nodes_combined(new_samp.adj, new_samp.anc, data_logmutrel)
  W_dests_new = _make_W_dests_combined(
    B_prime,
    new_samp.adj,
    new_samp.anc,
    data_logmutrel,
  )

  if common.debug.DEBUG:
    true_parent = _find_parent(B, common._true_adjm)
    old_parent = _find_parent(B, old_samp.adj)
    _generate_new_sample.debug = (
      (
        B,
        (A, true_parent, old_parent),
      ),
      (mode_node, mode_dest),
      W_nodes_old[0],
      W_nodes_old[1],
      W_dests_old[0],
      W_dests_old[1],
      '%.3f' % np.max(W_nodes_old[1]),
      '%.3f' % np.max(W_dests_old[1]),
    )

  log_p_B_new_given_old = np.log(np.dot(mode_node_weights, W_nodes_old[:,B]))
  log_p_A_new_given_old = np.log(np.dot(mode_dest_weights, W_dests_old[:,A]))
  # The need to use `A_prime` and `B_prime` here rather than `A` and `B`
  # becomes apparent when you consider the case when `B` is ancestral to `A` in
  # the old tree.
  log_p_B_old_given_new = np.log(np.dot(mode_node_weights, W_nodes_new[:,B_prime]))
  log_p_A_old_given_new = np.log(np.dot(mode_dest_weights, W_dests_new[:,A_prime]))

  log_p_new_given_old = log_p_B_new_given_old + log_p_A_new_given_old
  log_p_old_given_new = log_p_B_old_given_new + log_p_A_old_given_new
  return (new_samp, log_p_new_given_old, log_p_old_given_new)

def _run_chain(data_logmutrel, supervars, superclusters, nsamples, thinned_frac, phi_method, phi_iterations, seed, progress_queue=None):
  assert nsamples > 0

  V, N, omega_v = calc_binom_params(supervars)
  def __calc_phi(adj):
    phi, eta = phi_fitter.fit_phis(adj, superclusters, supervars, method=phi_method, iterations=phi_iterations, parallel=0)
    return phi
  def __calc_llh_phi(adj, phi):
    return _calc_llh_phi(phi, V, N, omega_v)

  samps = [_init_chain(seed, data_logmutrel, __calc_phi, __calc_llh_phi)]
  accepted = 0
  if progress_queue is not None:
    progress_queue.put(0)

  assert 0 < thinned_frac <= 1
  record_every = round(1 / thinned_frac)
  # Why is `expected_total_trees` equal to this?
  #
  # We always taken the first tree, since `0%k = 0` for all `k`. There remain
  # `nsamples - 1` samples to take, of which we record every `record_every`
  # one.
  #
  # This can give somewhat weird results, since you intuitively expect
  # approximately `thinned_frac * nsamples` trees to be returned. E.g., if
  # `nsamples = 3000` and `thinned_frac = 0.3`, you expect `0.3 * 3000 = 900`
  # trees, but you actually get 1000. To not be surprised by this, try to
  # choose `thinned_frac` such that `1 / thinned_frac` is close to an integer.
  # (I.e., `thinned_frac = 0.5` or `thinned_frac = 0.3333333` generally give
  # results as you'd expect.
  expected_total_trees = 1 + math.floor((nsamples - 1) / record_every)

  old_samp = samps[0]
  for I in range(1, nsamples):
    if progress_queue is not None:
      progress_queue.put(I)

    new_samp, log_p_new_given_old, log_p_old_given_new = _generate_new_sample(
      old_samp,
      data_logmutrel,
      __calc_phi,
      __calc_llh_phi,
    )
    log_p_transition = (new_samp.llh_phi - old_samp.llh_phi) + (log_p_old_given_new - log_p_new_given_old)
    U = np.random.uniform()
    accept = log_p_transition >= np.log(U)
    if accept:
      samp = new_samp
    else:
      samp = old_samp

    def _print_debug():
      if not common.debug.DEBUG:
        return
      true_adj, true_phi = common._true_adjm, common._true_phi
      norm_phi_llh = -old_samp.phi.size * np.log(2)
      cols = (
        'iter',
        'action',
        'old_llh',
        'new_llh',
        'true_llh',
        'p_new_given_old',
        'p_old_given_new',
        'old_parents',
        'new_parents',
        'true_parents',
        'node_error',
        'nodes',
        'sample_mode',
        'W_nodes_old_0',
        'W_nodes_old_1',
        'W_dests_old_0',
        'W_dests_old_1',
        'max_W_nodes_old_1',
        'max_W_dests_old_1',
      )
      vals = (
        I,
        'accept' if accept else 'reject',
        '%.3f' % (old_samp.llh_phi / norm_phi_llh),
        '%.3f' % (new_samp.llh_phi / norm_phi_llh),
        '%.3f' % (__calc_llh_phi(true_adj, true_phi) / norm_phi_llh),
        '%.3f' % log_p_new_given_old,
        '%.3f' % log_p_old_given_new,
        util.find_parents(old_samp.adj),
        util.find_parents(new_samp.adj),
        util.find_parents(true_adj),
        _make_W_nodes_mutrel.node_error,
      )
      vals = vals + _generate_new_sample.debug
      print(*['%s=%s' % (K, V) for K, V in zip(cols, vals)], sep='\t')
    _print_debug()

    if I % record_every == 0:
      samps.append(samp)
    old_samp = samp
    if accept:
      accepted += 1

  if nsamples > 1:
    accept_rate = accepted / (nsamples - 1)
  else:
    accept_rate = 1.
  assert len(samps) == expected_total_trees
  return (
    [S.adj     for S in samps],
    [S.phi     for S in samps],
    [S.llh_phi for S in samps],
    accept_rate,
  )

def use_existing_structures(adjms, supervars, superclusters, phi_method, phi_iterations, parallel=0):
  V, N, omega_v = calc_binom_params(supervars)
  K = len(supervars)
  phis = []
  llhs = []

  for adjm in adjms:
    assert adjm.shape == (K+1, K+1)
    phi, eta = phi_fitter.fit_phis(adjm, superclusters, supervars, method=phi_method, iterations=phi_iterations, parallel=parallel)
    llh = _calc_llh_phi(phi, V, N, omega_v)
    phis.append(phi)
    llhs.append(llh)
  return (np.array(adjms), np.array(phis), np.array(llhs))

def sample_trees(data_mutrel, supervars, superclusters, trees_per_chain, burnin, nchains, thinned_frac, phi_method, phi_iterations, seed, parallel):
  assert nchains > 0
  assert trees_per_chain > 0
  assert 0 <= burnin <= 1
  assert 0 < thinned_frac <= 1

  if common.debug.DEBUG:
    _load_truth(common.debug._truthfn)

  jobs = []
  total = nchains * trees_per_chain
  data_logmutrel = _make_data_logmutrel(data_mutrel)

  # Don't use (hard-to-debug) parallelism machinery unless necessary.
  if parallel > 0:
    import concurrent.futures
    import multiprocessing
    import queue
    import time
    import sys

    manager = multiprocessing.Manager()
    # What is stored in progress_queue doesn't matter. The queue is just used
    # so that child processes can signal when they've sampled a tree, allowing
    # the main process to update the progress bar.
    progress_queue = manager.Queue()

    with tqdm.tqdm(total=total, desc='Sampling trees', unit='tree', dynamic_ncols=True) as pbar:
      with concurrent.futures.ProcessPoolExecutor(max_workers=parallel) as ex:
        for C in range(nchains):
          # Ensure each chain's random seed is different from the seed used to
          # seed the initial Pairtree invocation, yet nonetheless reproducible.
          jobs.append(ex.submit(_run_chain, data_logmutrel, supervars, superclusters, trees_per_chain, thinned_frac, phi_method, phi_iterations, seed + C + 1, progress_queue))

        while True:
          finished = 0
          last_check = time.perf_counter()

          for J in jobs:
            if J.done():
              exception = J.exception(timeout=0.001)
              if exception is not None:
                # Ideally, if an exception occurs in a child process, we'd like
                # to re-raise the exception from the parent process and cause
                # the application to crash immediately. However,
                # `concurrent.futures` will wait until all child processes
                # terminate (either because they finished successfully or
                # suffered an exception) to raise the exception from the parent
                # process. If an exception occurs in the child process, it
                # would normally be raised automatically when we call
                # `J.result()` below; however, since we're in a loop, we must
                # break out of the loop, which we can do so by explicitly
                # re-raising the exception.
                #
                # Note that the `print()` call will happen immediately, so the
                # user will be notified as soon as the error occurs. Actually
                # raising the exception and crashing the application will not
                # occur, however, until all child processes finish. When the
                # `raise` happens below, the parent process will effectively
                # freeze at that statement until all child processes finish,
                # such that the progress bar stops updating. Ieally, we could
                # terminate all child processes when we detect the exception,
                # but we would have to do this manually by sending SIGTERM to
                # the processes. It's evidently impossible to get the PIDs of
                # the child processes without installing `psutil`, so we will
                # just allow the application to wait until child processes
                # finish before crashing.
                print('Exception occurred in child process:', exception, file=sys.stderr)
                raise exception
              else:
                finished += 1

          if finished == nchains:
            break
          while time.perf_counter() - last_check < 5:
            try:
              # If there's something in the queue for us to retrieve, a child
              # process has sampled a tree.
              progress_queue.get(timeout=5)
              pbar.update()
            except queue.Empty:
              pass

    results = [J.result() for J in jobs]
  else:
    results = []
    for C in range(nchains):
      results.append(_run_chain(data_logmutrel, supervars, superclusters, trees_per_chain, thinned_frac, phi_method, phi_iterations, seed + C + 1))

  merged_adj = []
  merged_phi = []
  merged_llh = []
  accept_rates = []
  for A, P, L, accept_rate in results:
    assert len(A) == len(P) == len(L) == len(results[0][0])
    discard_first = round(burnin * len(A))
    merged_adj += A[discard_first:]
    merged_phi += P[discard_first:]
    merged_llh += L[discard_first:]
    accept_rates.append(accept_rate)
  assert len(merged_adj) == len(merged_phi) == len(merged_llh)
  return (merged_adj, merged_phi, merged_llh, accept_rates)

def compute_posterior(adjms, phis, llhs, sort_by_llh=True):
  unique = {}

  for A, P, L in zip(adjms, phis, llhs):
    parents = util.find_parents(A)
    H = hash(parents.tobytes())
    if H in unique:
      assert np.isclose(L, unique[H]['llh'])
      # Use relaxed `atol`, or sometimes the phis (at least when computed by
      # `projection`) won't be close for two tree samples with the same
      # adjacency matrix. Identical tree structures with (slightly) different
      # phis can arise despite the caching mechanism that stores phis for each
      # tree structure. This occurs because different chains running on
      # different cores might sample the same tree structure, but the caching
      # mechanism is chain-specific. `projection` is not entirely
      # deterministic, so it may compute slightly different phis for the same
      # tree structure.
      assert np.allclose(P, unique[H]['phi'], atol=1e-5)
      assert np.array_equal(parents, unique[H]['struct'])
      unique[H]['count'] += 1
    else:
      unique[H] = {
        'struct': parents,
        'phi': P,
        'llh': L,
        'count': 1,
      }

  if sort_by_llh:
    unique = sorted(unique.values(), key = lambda T: -(np.log(T['count']) + T['llh']))
  else:
    unique = list(unique.values())
  unzipped = {key: np.array([U[key] for U in unique]) for key in unique[0].keys()}
  unzipped['prob'] = util.softmax(np.log(unzipped['count']) + unzipped['llh'])

  return (
    unzipped['struct'],
    unzipped['count'],
    unzipped['phi'],
    unzipped['llh'],
    unzipped['prob'],
  )
