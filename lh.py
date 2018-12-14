from common import Models, Variant, _EPSILON, _LOGEPSILON
import scipy.stats
import scipy.special
import scipy.integrate
import numpy as np
import util
import warnings

def generate_logprob_phi(N):
  prob = {}
  for modelidx in (Models.cocluster, Models.A_B, Models.B_A, Models.diff_branches):
    prob[modelidx] = np.zeros((N, N))

  for i in range(N):
    # cocluster
    prob[Models.cocluster][i,i] = 1./N
    for j in range(N):
      # B_A
      if i <= j:
        prob[Models.B_A][i,j] = 1. / (N*(N - 1)/2 + N)
      # A_B
      if i >= j:
        prob[Models.A_B][i,j] = 1. / (N*(N - 1)/2 + N)
      # diff_branches
      if i + j < N:
        prob[Models.diff_branches][i,j] = 1. / (N*(N - 1)/2 + N)

  logprob = {}
  for M in prob.keys():
    assert np.isclose(np.sum(prob[M]), 1)
    logprob[M] = np.zeros(prob[M].shape)
    logprob[M][prob[M] != 0] = np.log(prob[M][prob[M] != 0])
    logprob[M][prob[M] == 0] = _LOGEPSILON
    assert np.sum(logprob[M] == 0) == 0

  return logprob

def calc_lh_grid(var1, var2):
  grid_step = 0.001
  N = 3*int(1/grid_step + 1) # N: number of grid points
  G = np.linspace(start=0, stop=1, num=N)[:,np.newaxis] # Nx1

  S = len(var1.total_reads) # S
  logprob_phi = generate_logprob_phi(N) # MxNxN
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM

  for s in range(S):
    for modelidx in logprob_phi.keys():
      pv1, pv2 = [scipy.stats.binom.logpmf(V.var_reads[s], V.total_reads[s], V.omega_v[s]*G) for V in (var1, var2)] # Nx1
      p1 = np.tile(pv1, N)   # pv1 vector tiled as columns
      p2 = np.tile(pv2, N).T # pv2 vector tiled as rows
      P = p1 + p2 + logprob_phi[modelidx]
      logprob_models[s,modelidx] = scipy.special.logsumexp(P)

  return logprob_models

def _gen_samples(modelidx, S):
  U = scipy.stats.uniform.rvs(loc=0, scale=1, size=S)
  if modelidx == Models.A_B:
    phi1 = np.sqrt(U)
    phi2 = scipy.stats.uniform.rvs(loc=0, scale=phi1, size=S)
    area = 0.5
  elif modelidx == Models.B_A:
    phi2 = np.sqrt(U)
    phi1 = scipy.stats.uniform.rvs(loc=0, scale=phi2, size=S)
    area = 0.5
  elif modelidx == Models.diff_branches:
    phi1 = 1 - np.sqrt(1 - U)
    phi2 = scipy.stats.uniform.rvs(loc=0, scale=(1 - phi1), size=S)
    area = 0.5
  elif modelidx == Models.cocluster:
    phi1 = scipy.stats.uniform.rvs(loc=0, scale=1, size=S)
    phi2 = np.copy(phi1)
    area = 1
  else:
    raise Exception('Unknown model: %s' % model)
  return (phi1, phi2, area)

def calc_lh_mc_2D(var1, var2):
  S = len(var1.total_reads) # S
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM
  for s in range(S):
    for modelidx in (Models.cocluster, Models.A_B, Models.B_A, Models.diff_branches):
      mcsamps = int(1e6)
      phi1, phi2, area = _gen_samples(modelidx, S=mcsamps)

      logP = []
      for V, phi in ((var1, phi1), (var2, phi2)):
        logP.append(scipy.stats.binom.logpmf(V.var_reads[s], V.total_reads[s], V.omega_v[s]*phi))
      logprob_models[s,modelidx] = scipy.special.logsumexp(logP[0] + logP[1]) - np.log(mcsamps)

  return logprob_models

def calc_lh_mc_2D_dumb(var1, var2):
  S = len(var1.total_reads) # S
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM
  for s in range(S):
    for modelidx in (Models.cocluster, Models.A_B, Models.B_A, Models.diff_branches):
      mcsamps = int(1e6)
      phi1 = scipy.stats.uniform.rvs(loc=0, scale=1, size=mcsamps)
      phi2 = scipy.stats.uniform.rvs(loc=0, scale=1, size=mcsamps)

      if modelidx == Models.A_B:
        prior = phi1 >= phi2
        area = 0.5
      elif modelidx == Models.B_A:
        prior = phi2 >= phi1
        area = 0.5
      elif modelidx == Models.diff_branches:
        prior = phi1 + phi2 <= 1
        area = 0.5
      elif modelidx == Models.cocluster:
        phi2 = phi1
        prior = np.array(mcsamps * [True])
        area = 1
      else:
        raise Exception('Unknown model: %s' % model)

      logP = []
      for V, phi in ((var1, phi1), (var2, phi2)):
        logP.append(scipy.stats.binom.logpmf(V.var_reads[s], V.total_reads[s], V.omega_v[s]*phi))
      logprob_models[s,modelidx] = scipy.special.logsumexp((logP[0] + logP[1] - np.log(area))[prior]) - np.log(mcsamps)

  return logprob_models

def _calc_garbage_smart(V1, V2):
  S = len(V1.total_reads) # S
  evidence = np.nan * np.ones(S)

  for sidx in range(S):
    logP = [-np.log(V1.omega_v[sidx]), -np.log(V2.omega_v[sidx])]
    for V in (V1, V2):
      A, B = V.var_reads[sidx] + 1, V.ref_reads[sidx] + 1
      logP.append(util.log_N_choose_K(V.total_reads[sidx], V.var_reads[sidx]))
      logP.append(np.log(scipy.special.betainc(A, B, V.omega_v[sidx])))
      logP.append(scipy.special.betaln(A, B)) # Denormalization factor for beta
    evidence[sidx] = np.sum(logP)

  return evidence

def _calc_garbage_dumb(V1, V2):
  S = len(V1.total_reads) # S
  evidence = np.nan * np.ones(S)
  mcsamps = int(1e5)

  for sidx in range(S):
    logP = []
    for V in (V1, V2):
      phi = scipy.stats.uniform.rvs(loc=0, scale=1, size=mcsamps)
      binom = scipy.stats.binom.logpmf(V.var_reads[sidx], V.total_reads[sidx], V.omega_v[sidx]*phi)
      logP.append(scipy.special.logsumexp(binom) - np.log(mcsamps))
    evidence[sidx] = np.sum(logP)

  return evidence

def calc_lh_mc_1D(V1, V2):
  S = len(V1.total_reads) # S
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM
  for sidx in range(S):
    for modelidx in (Models.cocluster, Models.A_B, Models.B_A, Models.diff_branches):
      mcsamps = int(1e5)
      phi1 = scipy.stats.uniform.rvs(loc=0, scale=1, size=mcsamps)

      if modelidx == Models.cocluster:
        logP = []
        for V in (V1, V2):
          logP.append(scipy.stats.binom.logpmf(V.var_reads[sidx], V.total_reads[sidx], V.omega_v[sidx]*phi1))
        logprob_models[sidx,modelidx] = scipy.special.logsumexp(logP[0] + logP[1]) - np.log(mcsamps)
      else:
        logP = scipy.stats.binom.logpmf(V1.var_reads[sidx], V1.total_reads[sidx], V1.omega_v[sidx]*phi1)

        lower = _make_lower(phi1, modelidx)
        upper = _make_upper(phi1, modelidx)
        A = V2.var_reads[sidx] + 1
        B = V2.ref_reads[sidx] + 1

        betainc = [scipy.special.betainc(A, B, V2.omega_v[sidx]*limit) for limit in (upper, lower)]
        logdenorm = scipy.special.betaln(A, B)
        # Add epsilon to ensure we don't take log of zero.
        # TODO: check if betainc[0] is close to betainc[1]. If so, the log will
        # likely produce -inf, so we can just set the LLH to -inf.
        logP += np.log(betainc[0] - betainc[1] + _EPSILON) + logdenorm

        logP = scipy.special.logsumexp(logP) - np.log(mcsamps)
        logP += np.log(2)
        logP += util.log_N_choose_K(V2.total_reads[sidx], V2.var_reads[sidx])
        logP -= np.log(V2.omega_v[sidx])
        logprob_models[sidx,modelidx] = logP

  return logprob_models

def _make_lower(phi1, midx):
  ret = {
    Models.A_B: 0,
    Models.B_A: phi1,
    Models.diff_branches: 0,
  }[midx]
  return np.broadcast_to(ret, np.array(phi1).shape)

def _make_upper(phi1, midx):
  ret = {
    Models.A_B: phi1,
    Models.B_A: 1,
    Models.diff_branches: 1 - phi1,
  }[midx]
  return np.broadcast_to(ret, np.array(phi1).shape)

# SciPy's implementation is really slow. This is called hundreds of thousands
# of times, so I need something fast.
def binom_logpmf(X, N, P):
  assert 0 <= X <= N and 0 <= P <= 1, ('bad_binom X=%s N=%s P=%s' % (X, N, P))
  if np.isclose(0, P):
    if X == 0:
      return 0
    else:
      return -np.inf
  if np.isclose(1, P):
    if X == N:
      return 0
    else:
      return -np.inf
  return util.log_N_choose_K(N, X) + X*np.log(P) + (N - X)*np.log(1 - P)

def _integral_separate_clusters(phi1, V1, V2, sidx, midx, logsub=0):
  logP = binom_logpmf(
    V1.var_reads[sidx],
    V1.total_reads[sidx],
    V1.omega_v[sidx]*phi1
  )
  lower = _make_lower(phi1, midx)
  upper = _make_upper(phi1, midx)

  A = V2.var_reads[sidx] + 1
  B = V2.ref_reads[sidx] + 1
  betainc_upper = scipy.special.betainc(A, B, V2.omega_v[sidx] * upper)
  betainc_lower = scipy.special.betainc(A, B, V2.omega_v[sidx] * lower)
  if np.isclose(betainc_upper, betainc_lower):
    return 0
  logP += np.log(betainc_upper - betainc_lower)
  logP -= logsub

  return np.exp(logP)

def _integral_same_cluster(phi1, V1, V2, sidx, logsub=0):
  binom = [binom_logpmf(V.var_reads[sidx], V.total_reads[sidx], V.omega_v[sidx]*phi1) for V in (V1, V2)]
  logP = binom[0] + binom[1]
  logP -= logsub
  return np.exp(logP)

def quad(*args, **kwargs):
  with warnings.catch_warnings():
    warnings.simplefilter('ignore', category=scipy.integrate.IntegrationWarning)
    return scipy.integrate.quad(*args, **kwargs)

def calc_lh_quad(V1, V2):
  max_splits = 50
  S = len(V1.total_reads) # S
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM

  for sidx in range(S):
    V1_phi_mle = V1.vaf[sidx] / V1.omega_v[sidx]
    V1_phi_mle = np.minimum(1, V1_phi_mle)
    V1_phi_mle = np.maximum(0, V1_phi_mle)

    for modelidx in (Models.cocluster, Models.A_B, Models.B_A, Models.diff_branches):
      if modelidx == Models.cocluster:
        logmaxP = np.log(_integral_same_cluster(V1_phi_mle, V1, V2, sidx) + _EPSILON)
        P, P_error = quad(_integral_same_cluster, 0, 1, args=(V1, V2, sidx, logmaxP), limit=max_splits)
        P = np.maximum(_EPSILON, P)
        logP = np.log(P) + logmaxP
      else:
        logmaxP = np.log(_integral_separate_clusters(V1_phi_mle, V1, V2, sidx, modelidx) + _EPSILON)
        P, P_error = quad(_integral_separate_clusters, 0, 1, args=(V1, V2, sidx, modelidx, logmaxP), limit=max_splits)
        logdenorm = scipy.special.betaln(V2.var_reads[sidx] + 1, V2.ref_reads[sidx] + 1)
        P = np.maximum(_EPSILON, P)
        logP = np.log(P) + logmaxP + logdenorm + np.log(2) + util.log_N_choose_K(V2.total_reads[sidx], V2.var_reads[sidx]) - np.log(V2.omega_v[sidx])
      logprob_models[sidx,modelidx] = logP
  return logprob_models

def _find_bad_samples(V1, V2):
  read_threshold = 3
  omega_threshold = 1e-3

  # In cases where we have zero variant reads for both variants, the garbage
  # models ends up taking considerable posterior mass, but we have no information
  # to judge whether the pair should be garbage. The contribution from lots of
  # samples with zero variant reads on both variants can end up overwhelming
  # informative samples with non-zero variant reads, such that the pair across
  # samples jointly is deemed garbage. This is undesireable behaviour, clearly,
  # and so ignore such (0, 0) cases by zeroing out their evidence. This
  # effectively states that we have a uniform posterior over the relations of
  # those mutations in that sample, which is the sensible thing to do.
  #
  # See this discussion on Slack for more details:
  # https://morrislab.slack.com/archives/CCE5HNVSP/p1540392641000100
  #
  # This also seems to apply to cases where both variants have only a couple
  # reads, meaning most of their 2D binomial mass is close to the origin. To
  # reudce the garbage model's tendency to win in these cases, don't allow them
  # to contribute to the posterior, either.
  too_few_var_reads = np.logical_and(V1.var_reads < read_threshold, V2.var_reads < read_threshold)

  # If omega_v is below threhsold, we conclude that all variant reads are
  # noise, and that this variant thus confers no evidence concerning its
  # relationship with other variants within the sample in question. In such
  # instances, set uniform evidence for this relationship.
  uninformative_omega = np.logical_or(V1.omega_v < omega_threshold, V2.omega_v < omega_threshold)

  bad_samples = np.logical_and(too_few_var_reads, uninformative_omega)
  return bad_samples

def _filter_samples(V, samp_filter):
  return Variant(
    id = V.id,
    ref_reads = V.ref_reads[samp_filter],
    var_reads = V.var_reads[samp_filter],
    total_reads = V.total_reads[samp_filter],
    vaf = V.vaf[samp_filter],
    omega_v = V.omega_v[samp_filter],
  )

calc_garbage = _calc_garbage_smart
#calc_garbage = _calc_garbage_dumb

def calc_lh(V1, V2, _calc_lh=None):
  if _calc_lh is None:
    _calc_lh = calc_lh_quad

  S = len(V1.omega_v)
  evidence_per_sample = np.zeros((S, len(Models._all)))

  if V1.id == V2.id:
    # If they're the same variant, they should cocluster with certainty.
    evidence_per_sample += -np.inf
    evidence_per_sample[:,Models.cocluster] = 0
    evidence = evidence_per_sample[0]
    return (evidence, evidence_per_sample)

  bad_samples = _find_bad_samples(V1, V2)
  good_samples = np.logical_not(bad_samples)
  V1 = _filter_samples(V1, good_samples)
  V2 = _filter_samples(V2, good_samples)

  evidence_per_sample[good_samples] = _calc_lh(V1, V2)
  evidence_per_sample[good_samples,Models.garbage] = calc_garbage(V1, V2)

  evidence = np.sum(evidence_per_sample, axis=0)
  return (evidence, evidence_per_sample)
