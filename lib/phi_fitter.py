import common
import numpy as np

def fit_phis(adj, superclusters, supervars, method, iterations, parallel):
  if method == 'debug':
    # Bypass cache when debugging.
    return _fit_phis(adj, superclusters, supervars, method, iterations, parallel)
  key = (hash(adj.tobytes()), iterations)
  if key not in fit_phis.cache:
    fit_phis.cache[key] = _fit_phis(adj, superclusters, supervars, method, iterations, parallel)
    fit_phis.cache_misses += 1
  else:
    fit_phis.cache_hits += 1
  return fit_phis.cache[key]

fit_phis.cache = {}
fit_phis.cache_hits = 0
fit_phis.cache_misses = 0

# Used only for `rprop_cached`.
last_eta = ['mle']

def _calc_llh(phi, V, N, omega_v, epsilon=1e-5):
  import scipy
  K, S = phi.shape
  for arr in V, N, omega_v:
    assert arr.shape == (K-1, S)

  assert np.allclose(1, phi[0])
  P = omega_v * phi[1:]
  P = np.maximum(P, epsilon)
  P = np.minimum(P, 1 - epsilon)

  phi_llh = scipy.stats.binom.logpmf(V, N, P) / np.log(2)
  assert not np.any(np.isnan(phi_llh))
  assert not np.any(np.isinf(phi_llh))

  llh_per_sample = -np.sum(phi_llh, axis=0) / K
  nlglh = np.sum(llh_per_sample) / S
  return (phi_llh, llh_per_sample, nlglh)

def lpdist(A, B, p=1):
  return np.sum(np.abs(A - B)**p)**(1/p)

def _fit_phis(adj, superclusters, supervars, method, iterations, parallel):
  # Calling `import` on each function call should be cheap, as Python caches a
  # reference to the module after the first load.
  if method in ('graddesc_old', 'rprop_old'):
    import phi_fitter_iterative
    eta = phi_fitter_iterative.fit_etas(adj, superclusters, supervars, method[:-4], iterations, parallel)

  elif method == 'rprop':
    import phi_fitter_lol
    eta = phi_fitter_lol.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel, eta_init='mle')

  elif method == 'projection':
    import phi_fitter_projection
    eta = phi_fitter_projection.fit_etas(adj, superclusters, supervars)

  elif method == 'proj_rprop':
    import phi_fitter_projection
    import phi_fitter_lol
    eta_proj = phi_fitter_projection.fit_etas(adj, superclusters, supervars)
    eta = phi_fitter_lol.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel, eta_init=eta_proj)

  elif method == 'debug':
    import phi_fitter_iterative
    import phi_fitter_projection
    import phi_fitter_lol
    import time
    fitters = {
      #'rprop_init_mle': lambda: phi_fitter_iterative.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel, eta_init=None),
      'lol_init_mle': lambda: phi_fitter_lol.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel, eta_init='mle'),
      #'lol_init_dirichlet': lambda: phi_fitter_lol.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel, eta_init='dirichlet'),
      'projection': lambda: phi_fitter_projection.fit_etas(adj, superclusters, supervars),
    }
    #fitters['lol_init_proj'] = lambda: phi_fitter_lol.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel, eta_init=fitters['projection']())
    #fitters['lol_init_prev'] = lambda: phi_fitter_lol.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel, eta_init=last_eta[0])

    Z = common.make_ancestral_from_adj(adj)
    svids = common.extract_vids(supervars)
    total_reads = np.array([supervars[svid]['total_reads'] for svid in svids])
    var_reads = np.array([supervars[svid]['var_reads'] for svid in svids])
    omega = np.array([supervars[svid]['omega_v'] for svid in svids])

    etas = {}
    scores = {}
    times = {}
    zeros = {}
    l1_dists = {}
    l2_dists = {}
    for name, F in fitters.items():
      time_start = time.perf_counter_ns()
      etas[name] = F()
      time_end = time.perf_counter_ns()
      phi = np.dot(Z, etas[name])
      scores[name] = _calc_llh(phi, var_reads, total_reads, omega)
      times[name] = (time_end - time_start)/1e6
      zeros[name] = np.sum(phi == 0)
      l1_dists[name] = lpdist(var_reads/(total_reads * omega), phi[1:], p=1)
      l2_dists[name] = lpdist(var_reads/(total_reads * omega), phi[1:], p=2)

    eta = etas['lol_init_mle']
    last_eta[0] = np.copy(eta)

    names = sorted(etas.keys())
    sep = '\t'
    if True and not hasattr(_fit_phis, 'printed_header'):
      print(*names, sep=sep)
      _fit_phis.printed_header = True
    print(
      *['%.3f' % scores[name][2] for name in names],
      np.nan,
      *['%.3f' % times[name] for name in names],
      np.nan,
      *[zeros[name] for name in names],
      np.nan,
      *['%.3f' % l1_dists[name] for name in names],
      np.nan,
      *['%.3f' % l2_dists[name] for name in names],
      sep=sep,
      flush=True
    )

  else:
    raise Exception('Unknown phi fitter %s' % method)

  assert np.allclose(1, np.sum(eta, axis=0))
  Z = common.make_ancestral_from_adj(adj)
  phi = np.dot(Z, eta)
  return (phi, eta)
