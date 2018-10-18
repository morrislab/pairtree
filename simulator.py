import numpy as np
import scipy.stats

def make_ancestral_from_adj(adj):
  K = len(adj)
  Z = np.zeros((K,K))
  adj = np.copy(adj)
  np.fill_diagonal(adj, 1)

  def _find_desc(I, vec):
    # Base case: if I have no children, my ancestor vec is just myself.
    if np.sum(vec) == 1:
      return vec
    else:
      children = np.array([_find_desc(idx, adj[idx]) for (idx, val) in enumerate(vec) if idx != I and val == 1])
      self_and_child = vec + np.sum(children, axis=0)
      self_and_child[self_and_child > 1] = 1
      return self_and_child

  for k in range(K):
    # If we know `adj` is topologically sorted, we can reduce the complexity of
    # this -- we would start at leaves and work our way upward, eliminating
    # need for recursive DFS. But since we don't expect `K` to be large, we can
    # write more general version that works for non-toposorted trees.
    Z[k] = _find_desc(k, adj[k])

  return Z

def make_parents(K):
  # Determine parents of nodes [1, 2, ..., K].
  parents = []
  # mu is the probability of extending the current branch.
  mu = 0.75
  for idx in range(K - 1):
    U = np.random.uniform()
    if U < mu:
      parents.append(len(parents))
    else:
      parents.append(np.random.randint(0, idx + 1))
  return parents

def make_adjm(K):
  parents = make_parents(K)
  # Make adjacency matrix from parents vector.
  adjm = np.zeros((K, K))
  adjm[parents, range(1, K)] = 1
  return adjm

def generate_tree(K, S):
  adjm = make_adjm(K) # KxK
  #leaves = np.flatnonzero(np.sum(adjm, axis=1) == 0)
  Z = make_ancestral_from_adj(adjm) # KXK

  eta = np.random.dirichlet(alpha = K*[1e0], size = S).T # KxS
  # In general, we want etas on leaves to be more "peaked" -- that is, only a
  # few subclones come to dominate, so they should have large etas relative to
  # internal nodes. We accomplish this by using a smaller alpha for these.
  #eta[leaves] += np.random.dirichlet(alpha = len(leaves)*[1e0], size = S).T
  eta /= np.sum(eta, axis=0)

  phi = np.dot(Z, eta) # KxS
  return (adjm, phi)

def generate_read_counts(phi, omega_v, T):
  M, S = phi.shape
  # T: total reads. Broadcast operation ensures V and T are same shape.
  T = np.broadcast_to(T, (M,S))
  # V: variant reads
  V = scipy.stats.binom.rvs(n=T, p=omega_v*phi)
  return (V, T)

def add_noise(mat, sigma=0.09):
  noisy = np.random.normal(loc=mat, scale=sigma)
  capped = np.maximum(0, np.minimum(1, noisy))
  return capped

def add_garbage(mat, N):
  K, S = mat.shape
  garbage = np.random.uniform(size=(N, S))
  return np.vstack((mat, garbage))

def choose_mutass(K, M):
  # Ensure every cluster has at least one mutation.
  assert M >= K
  first_mutass = np.arange(K)
  if M == K:
    return first_mutass

  # P: probability vector indicating what proportion of mutations should be
  # assigned to each cluster.
  P = np.random.dirichlet(alpha = K*[1e0])
  remaining_mutass = np.random.choice(K, p=P, size=(M - K))
  
  mutass = np.concatenate((first_mutass, remaining_mutass))
  np.random.shuffle(mutass)
  return mutass

def generate_data(K, S, T, M):
  # K: number of clusters
  # S: number of samples
  # T: reads per mutation
  # M: total number of mutations
  adjm, phi = generate_tree(K, S)
  mutass = choose_mutass(K, M)
  phi_mutations = np.array([phi[cluster] for cluster in mutass])
  omega_v = 0.5
  V, T = generate_read_counts(phi_mutations, omega_v, T)
  
  return {
    'adjm': adjm,
    'phi': phi,
    'mutass': mutass,
    'omega_v': np.broadcast_to(omega_v, M),
    'V': V,
    'T': T,
    'sampnames': ['Sample %s' % (sidx + 1) for sidx in range(S)],
  }
