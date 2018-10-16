from common import Models
import scipy.stats
import scipy.special
import scipy.integrate
import numpy as np
import util

_LOGEPSILON = -400
_EPSILON    = np.exp(_LOGEPSILON)

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

  S = len(var1['total_reads']) # S
  logprob_phi = generate_logprob_phi(N) # MxNxN
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM

  for s in range(S):
    for modelidx in logprob_phi.keys():
      pv1, pv2 = [scipy.stats.binom.logpmf(V['var_reads'][s], V['total_reads'][s], (1 - V['mu_v'])*G) for V in (var1, var2)] # Nx1
      p1 = np.tile(pv1, N)   # pv1 vector tiled as columns
      p2 = np.tile(pv2, N).T # pv2 vector tiled as rows
      P = p1 + p2 + logprob_phi[modelidx]
      logprob_models[s,modelidx] = scipy.special.logsumexp(P)

  return logprob_models

def _gen_samples(modelidx, S):
  U = scipy.stats.uniform.rvs(loc=0, scale=1, size=S)
  if modelidx == Models.A_B:
    phi2 = np.sqrt(U)
    phi1 = scipy.stats.uniform.rvs(loc=0, scale=phi2, size=S)
    area = 0.5
  elif modelidx == Models.B_A:
    phi1 = np.sqrt(U)
    phi2 = scipy.stats.uniform.rvs(loc=0, scale=phi1, size=S)
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
  S = len(var1['total_reads']) # S
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM
  for s in range(S):
    for modelidx in (Models.cocluster, Models.A_B, Models.B_A, Models.diff_branches):
      mcsamps = int(1e6)
      phi1, phi2, area = _gen_samples(modelidx, S=mcsamps)

      logP = []
      for V, phi in ((var1, phi1), (var2, phi2)):
        logP.append(scipy.stats.binom.logpmf(V['var_reads'][s], V['total_reads'][s], (1 - V['mu_v'])*phi))
      logprob_models[s,modelidx] = scipy.special.logsumexp(logP[0] + logP[1]) - np.log(mcsamps)

  return logprob_models

def calc_lh_mc_2D_dumb(var1, var2):
  S = len(var1['total_reads']) # S
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
        logP.append(scipy.stats.binom.logpmf(V['var_reads'][s], V['total_reads'][s], (1 - V['mu_v'])*phi))
      logprob_models[s,modelidx] = scipy.special.logsumexp((logP[0] + logP[1] - np.log(area))[prior]) - np.log(mcsamps)

  return logprob_models

def calc_lh_mc_1D(V1, V2):
  S = len(V1['total_reads']) # S
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM
  for sidx in range(S):
    for modelidx in (Models.cocluster, Models.A_B, Models.B_A, Models.diff_branches):
      mcsamps = int(1e5)
      phi1 = scipy.stats.uniform.rvs(loc=0, scale=1, size=mcsamps)

      if modelidx == Models.cocluster:
        logP = []
        for V in (V1, V2):
          logP.append(scipy.stats.binom.logpmf(V['var_reads'][sidx], V['total_reads'][sidx], (1 - V['mu_v'])*phi1))
        logprob_models[sidx,modelidx] = scipy.special.logsumexp(logP[0] + logP[1]) - np.log(mcsamps)
      else:
        logP = scipy.stats.binom.logpmf(V1['var_reads'][sidx], V1['total_reads'][sidx], (1 - V1['mu_v'])*phi1)
        logP += np.log(4)
        logP += util.log_N_choose_K(V2['total_reads'][sidx], V2['var_reads'][sidx])

        lower = _make_lower(phi1, modelidx)
        upper = _make_upper(phi1, modelidx)
        A = V2['var_reads'][sidx] + 1
        B = V2['ref_reads'][sidx] + 1

        betainc = [scipy.special.betainc(A, B, 0.5*limit) for limit in (upper, lower)]
        logdenorm = scipy.special.betaln(V2['var_reads'][sidx] + 1, V2['ref_reads'][sidx] + 1)
        # Add delta to ensure we don't take log of zero.
        logP += np.log(betainc[0] - betainc[1] + _EPSILON) + logdenorm

        logprob_models[sidx,modelidx] = scipy.special.logsumexp(logP) - np.log(mcsamps)

  return logprob_models

def R(X, A, B, Z):
  K = int(0.5 * X)
  if X % 2 == 0:
    return np.log(K) + np.log(B - K) + np.log(Z) - np.log(A + 2*K - 1) - np.log(A + 2*K)
  else:
    return np.log(A + K) + np.log(A + B + K) + np.log(Z) - np.log(A + 2*K) - np.log(A + 2*K + 1)

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

def _integral_separate_clusters(phi1, V1, V2, sidx, midx, logsub=None):
  logP = scipy.stats.binom.logpmf(
    V1['var_reads'][sidx],
    V1['total_reads'][sidx],
    (1 - V1['mu_v'])*phi1
  )
  lower = _make_lower(phi1, midx)
  upper = _make_upper(phi1, midx)

  A = V2['var_reads'][sidx] + 1
  B = V2['ref_reads'][sidx] + 1
  betainc = [scipy.special.betainc(A, B, 0.5*limit) for limit in (upper, lower)]
  # Add delta to ensure we don't take log of zero.
  logP += np.log(betainc[0] - betainc[1] + _EPSILON)
  if logsub is not None:
    logP -= logsub

  return np.exp(logP)

def _integral_same_cluster(phi1, V1, V2, sidx, midx, logsub=None):
  binom = [scipy.stats.binom.logpmf( V['var_reads'][sidx], V['total_reads'][sidx], (1 - V['mu_v'])*phi1) for V in (V1, V2)]
  logP = binom[0] + binom[1]
  if logsub is not None:
    logP -= logsub
  return np.exp(logP)

def calc_lh_quad(V1, V2):
  max_splits = 50
  S = len(V1['total_reads']) # S
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM

  for sidx in range(S):
    for modelidx in (Models.cocluster, Models.A_B, Models.B_A, Models.diff_branches):
      if modelidx == Models.cocluster:
        logmaxP = np.log(_integral_same_cluster(V1['vaf'][sidx], V1, V2, sidx, modelidx) + _EPSILON)
        P, P_error = scipy.integrate.quad(_integral_same_cluster, 0, 1, args=(V1, V2, sidx, modelidx, logmaxP), limit=max_splits)
        P = np.maximum(_EPSILON, P)
        logP = np.log(P) + logmaxP
      else:
        logmaxP = np.log(_integral_separate_clusters(V1['vaf'][sidx], V1, V2, sidx, modelidx) + _EPSILON)
        P, P_error = scipy.integrate.quad(_integral_separate_clusters, 0, 1, args=(V1, V2, sidx, modelidx, logmaxP), limit=max_splits)
        logdenorm = scipy.special.betaln(V2['var_reads'][sidx] + 1, V2['ref_reads'][sidx] + 1)
        P = np.maximum(_EPSILON, P)
        # TODO: this factor of 2 is presumably wrong when mu_v != 1/2
        logP = np.log(P) + logmaxP + logdenorm + np.log(4) + util.log_N_choose_K(V2['total_reads'][sidx], V2['var_reads'][sidx])
      logprob_models[sidx,modelidx] = logP
  return logprob_models
