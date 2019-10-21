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
last_eta = [None]

def _fit_phis(adj, superclusters, supervars, method, iterations, parallel):
  # Calling `import` on each function call should be cheap, as Python caches a
  # reference to the module after the first load.
  if method in ('graddesc', 'rprop'):
    import phi_fitter_iterative
    eta = phi_fitter_iterative.fit_etas(adj, superclusters, supervars, method, iterations, parallel)

  elif method == 'projection':
    import phi_fitter_projection
    eta = phi_fitter_projection.fit_etas(adj, superclusters, supervars)

  elif method == 'proj_rprop':
    import phi_fitter_projection
    import phi_fitter_iterative
    eta_proj = phi_fitter_projection.fit_etas(adj, superclusters, supervars)
    eta = phi_fitter_iterative.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel, eta_init=eta_proj)

  elif method == 'debug':
    import phi_fitter_iterative
    import phi_fitter_projection
    import phi_fitter_lol
    import tree_sampler
    E = [
      ('rprop_init_mle', phi_fitter_iterative.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel, eta_init=None)),
      #('rprop_init_dirichlet', phi_fitter_iterative.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel, eta_init='dirichlet')),
      #('rprop_init_lasteta', phi_fitter_iterative.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel, eta_init=last_eta[0])),
      ('lol_init_mle', phi_fitter_lol.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel)),
      #('graddesc_init_mle', phi_fitter_iterative.fit_etas(adj, superclusters, supervars, 'graddesc', iterations, parallel, eta_init=None)),
      #('graddesc_init_dirichlet', phi_fitter_iterative.fit_etas(adj, superclusters, supervars, 'graddesc', iterations, parallel, eta_init='dirichlet')),
      #('graddesc_init_lasteta', phi_fitter_iterative.fit_etas(adj, superclusters, supervars, 'graddesc', iterations, parallel, eta_init=last_eta[0])),
      #('graddesc_numerical', phi_fitter_iterative.fit_etas(adj, superclusters, supervars, 'graddesc_numerical', iterations, parallel, eta_init=None)),
      ('projection', phi_fitter_projection.fit_etas(adj, superclusters, supervars)),
    ]
    #E.insert(3, ('rprop_init_proj', phi_fitter_iterative.fit_etas(adj, superclusters, supervars, 'rprop', iterations, parallel, eta_init=E[-1][1])))
    eta = E[0][1]
    last_eta[0] = np.copy(eta)

    Z = common.make_ancestral_from_adj(adj)
    svids = common.extract_vids(supervars)
    ref_reads = np.array([supervars[svid]['ref_reads'] for svid in svids])
    var_reads = np.array([supervars[svid]['var_reads'] for svid in svids])
    omega = np.array([supervars[svid]['omega_v'] for svid in svids])
    # `phi_fitter_iterative` assumes `omega` is 0.5, I think. I should fix this.
    assert np.all(omega == 0.5)

    K, S = eta.shape
    nlglh = [(name, -tree_sampler._calc_llh_phi(np.dot(Z, E_prime), var_reads, ref_reads + var_reads, omega) / (K*S*np.log(2))) for name, E_prime in E]
    names, scores = zip(*nlglh)
    if True and not hasattr(_fit_phis, 'printed_header'):
      print(*names, sep=',')
      _fit_phis.printed_header = True
    print(*scores, np.argsort(scores), sep=',', flush=True)

  else:
    raise Exception('Unknown phi fitter %s' % method)

  assert np.allclose(1, np.sum(eta, axis=0))
  Z = common.make_ancestral_from_adj(adj)
  phi = np.dot(Z, eta)
  return (phi, eta)
