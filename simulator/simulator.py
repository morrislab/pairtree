import numpy as np
import scipy.stats
from collections import defaultdict, OrderedDict
import common

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
  return np.array(parents)

def make_adjm(K, tree_type):
  assert tree_type in (None, 'monoprimary', 'polyprimary')
  while True:
    parents = make_parents(K)
    root_children = np.sum(parents == 0)
    if tree_type is None:
      break
    elif tree_type == 'monoprimary' and root_children == 1:
      assert parents[0] == 0 and not np.any(parents[1:] == 0)
      break
    elif tree_type == 'polyprimary' and root_children > 1:
      break

  # Make adjacency matrix from parents vector.
  adjm = np.eye(K)
  adjm[parents, range(1, K)] = 1
  return adjm

def generate_tree(K, S, tree_type):
  adjm = make_adjm(K, tree_type) # KxK
  #leaves = np.flatnonzero(np.sum(adjm, axis=1) == 0)
  Z = common.make_ancestral_from_adj(adjm) # KXK

  eta = np.random.dirichlet(alpha = K*[1e0], size = S).T # KxS
  # In general, we want etas on leaves to be more "peaked" -- that is, only a
  # few subclones come to dominate, so they should have large etas relative to
  # internal nodes. We accomplish this by using a smaller alpha for these.
  #eta[leaves] += np.random.dirichlet(alpha = len(leaves)*[1e0], size = S).T
  eta /= np.sum(eta, axis=0)

  phi = np.dot(Z, eta) # KxS
  assert np.allclose(1, phi[0])
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

def make_clusters(mutass):
  clusters = defaultdict(list)
  for midx, cidx in enumerate(mutass):
    clusters[cidx].append('s%s' % midx)
  assert set(clusters.keys()) == set(range(1, len(clusters) + 1))

  clusters = [clusters[cidx] for cidx in sorted(clusters.keys())]
  return clusters

def make_variants(V, T, omega_v):
  variants = OrderedDict()
  for midx in range(len(omega_v)):
    variant = {
      'id': 's%s' % midx,
      'name': 'S_%s' % midx,
      'var_reads': V[midx],
      'total_reads': T[midx],
      'omega_v': omega_v[midx],
    }
    variant['ref_reads'] = variant['total_reads'] - variant['var_reads']
    variant['vaf'] = variant['var_reads'] / variant['total_reads']
    variants[variant['id']] = variant
  return variants

def generate_data(K, S, T, M, G, tree_type):
  # K: number of clusters (excluding normal root)
  # S: number of samples
  # T: reads per mutation
  # M: total number of mutations
  # G: number of (additional) garbage mutations
  adjm, phi = generate_tree(K + 1, S, tree_type)
  # Add 1 to each mutation's assignment to account for normal root.
  mutass = choose_mutass(K, M) + 1 # Mx1
  clusters = make_clusters(mutass)

  phi_good_mutations = np.array([phi[cidx] for cidx in mutass]) # MxS
  phi_garbage = np.random.uniform(size=(G,S))
  phi_mutations = np.vstack((phi_good_mutations, phi_garbage))

  omega_v = np.broadcast_to(0.5, (M + G, S))
  V, T = generate_read_counts(phi_mutations, omega_v, T)

  variants = make_variants(V, T, omega_v)
  vids_good = ['s%s' % vidx for vidx in range(M)]
  vids_garbage = ['s%s' % vidx for vidx in range(M, M + G)]
  assert set(vids_good) == set([V for C in clusters for V in C])

  return {
    'adjm': adjm,
    'phi': phi,
    'clusters': clusters,
    'variants': variants,
    'vids_good': vids_good,
    'vids_garbage': vids_garbage,
    'sampnames': ['Sample %s' % (sidx + 1) for sidx in range(S)],
  }
