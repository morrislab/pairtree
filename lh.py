from common import Models
import scipy.stats
import scipy.special
import scipy.integrate
import numpy as np
import util

def generate_logprob_phi(N):
  prob = {}
  for modelidx, model in enumerate(Models._all):
    if modelidx == Models.garbage:
      continue
    prob[model] = np.zeros((N, N))

  for i in range(N):
    # cocluster
    prob['cocluster'][i,i] = 1./N
    for j in range(N):
      # B_A
      if i <= j:
        prob['B_A'][i,j] = 1. / (N*(N - 1)/2 + N)
      # A_B
      if i >= j:
        prob['A_B'][i,j] = 1. / (N*(N - 1)/2 + N)
      # diff_branches
      if i + j < N:
        prob['diff_branches'][i,j] = 1. / (N*(N - 1)/2 + N)

  logprob = {}
  for M in prob.keys():
    assert np.isclose(np.sum(prob[M]), 1)
    logprob[M] = np.zeros(prob[M].shape)
    logprob[M][prob[M] != 0] = np.log(prob[M][prob[M] != 0])
    logprob[M][prob[M] == 0] = -300
    assert np.sum(logprob[M] == 0) == 0

  return logprob

def calc_lh_binom_grid(var1, var2):
  grid_step = 0.001
  N = 3*int(1/grid_step + 1) # N: number of grid points
  G = np.linspace(start=0, stop=1, num=N)[:,np.newaxis] # Nx1

  S = len(var1['total_reads']) # S
  logprob_phi = generate_logprob_phi(N) # MxNxN
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM

  for s in range(S):
    for modelidx, model in enumerate(Models._all):
      if modelidx == Models.garbage:
        continue
      pv1, pv2 = [scipy.stats.binom.logpmf(V['var_reads'][s], V['total_reads'][s], np.minimum(1., (1 - V['mu_v'])*G)) for V in (var1, var2)] # Nx1
      p1 = np.tile(pv1, N)   # pv1 vector tiled as columns
      p2 = np.tile(pv2, N).T # pv1 vector tiled as rows
      P = p1 + p2 + logprob_phi[model]
      logprob_models[s,modelidx] = scipy.misc.logsumexp(P)

  return logprob_models

def _gen_samples(modelidx, S):
  if modelidx == Models.A_B:
    phi1 = scipy.stats.uniform.rvs(loc=0, scale=1, size=S)
    phi2 = scipy.stats.uniform.rvs(loc=0, scale=phi1, size=S)
    area = 0.5
  elif modelidx == Models.B_A:
    phi2 = scipy.stats.uniform.rvs(loc=0, scale=1, size=S)
    phi1 = scipy.stats.uniform.rvs(loc=0, scale=phi2, size=S)
    area = 0.5
  elif modelidx == Models.diff_branches:
    phi1 = scipy.stats.uniform.rvs(loc=0, scale=1, size=S)
    phi2 = scipy.stats.uniform.rvs(loc=0, scale=(1 - phi1), size=S)
    area = 0.5
  elif modelidx == Models.cocluster:
    phi1 = scipy.stats.uniform.rvs(loc=0, scale=1, size=S)
    phi2 = np.copy(phi1)
    area = 1
  else:
    raise Exception('Unknown model: %s' % model)
  return (phi1, phi2, area)

def calc_lh_binom_mc_2D(var1, var2):
  S = len(var1['total_reads']) # S
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM
  for s in range(S):
    for modelidx, model in enumerate(Models._all):
      if modelidx == Models.garbage:
        continue
      mcsamps = 1000000
      phi1, phi2, area = _gen_samples(modelidx, S=mcsamps)

      logP = []
      for V, phi in ((var1, phi1), (var2, phi2)):
        logP.append(scipy.stats.binom.logpmf(V['var_reads'][s], V['total_reads'][s], (1 - V['mu_v'])*phi))
      logprob_models[s,modelidx] = scipy.misc.logsumexp(logP[0] + logP[1]) - np.log(mcsamps)

  return logprob_models

def calc_lh_binom_mc_1D(V1, V2):
  S = len(V1['total_reads']) # S
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM
  for sidx in range(S):
    for modelidx, model in enumerate(Models._all):
      if modelidx == Models.garbage or modelidx == Models.cocluster:
        continue
      mcsamps = 1000000
      phi1 = scipy.stats.uniform.rvs(loc=0, scale=1, size=mcsamps)

      logP = scipy.stats.binom.logpmf(V1['var_reads'][sidx], V1['total_reads'][sidx], (1 - V1['mu_v'])*phi1)
      logP += np.log(2)
      logP += util.log_N_choose_K(V2['total_reads'][sidx], V2['var_reads'][sidx])

      lower = _make_lower(phi1, modelidx)
      upper = _make_upper(phi1, modelidx)
      A = V2['var_reads'][sidx] + 1
      B = V2['ref_reads'][sidx] + 1
      betainc = [scipy.special.betainc(A, B, 0.5*limit) for limit in (upper, lower)]
      logdenorm = scipy.special.betaln(V2['var_reads'][sidx] + 1, V2['ref_reads'][sidx] + 1)
      # Add delta to ensure we don't take log of zero.
      logP += np.log(betainc[0] - betainc[1] + 1e-300) + logdenorm

      logprob_models[sidx,modelidx] = scipy.misc.logsumexp(logP) - np.log(mcsamps)

  return logprob_models

def _make_lower(phi1, midx):
  return {
    Models.A_B: 0,
    Models.B_A: phi1,
    Models.diff_branches: 0,
  }[midx]

def _make_upper(phi1, midx):
  return {
    Models.A_B: phi1,
    Models.B_A: 1,
    Models.diff_branches: 1 - phi1,
  }[midx]

def _integral(phi1, V1, V2, sidx, midx, logsub=None):
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
  logP += np.log(betainc[0] - betainc[1] + 1e-300)
  if logsub is not None:
    logP -= logsub

  return np.exp(logP)

def calc_lh_binom_quad(V1, V2):
  S = len(V1['total_reads']) # S
  logprob_models = np.nan * np.ones((S, len(Models._all))) # SxM
  for sidx in range(S):
    for modelidx, model in enumerate(Models._all):
      if modelidx == Models.garbage:
        continue
      elif modelidx == Models.cocluster:
        P, P_error = scipy.integrate.quad(lambda phi1: np.prod([scipy.stats.binom.pmf(
          V['var_reads'][sidx],
          V['total_reads'][sidx],
          (1 - V['mu_v'])*phi1) for V in (V1, V2)]), 0, 1)
        P = np.maximum(1e-100, P)
        logprob_models[sidx,modelidx] = np.log(P)
      else:
        logmaxP = np.log(_integral(V1['vaf'][sidx], V1, V2, sidx, modelidx) + 1e-300)
        P, P_error = scipy.integrate.quad(_integral, 0, 1, args=(V1, V2, sidx, modelidx, logmaxP))
        logdenorm = scipy.special.betaln(V2['var_reads'][sidx] + 1, V2['ref_reads'][sidx] + 1)
        P = np.maximum(np.exp(-300), P)
        # TODO: this factor of 2 is presumably wrong when mu_v != 1/2
        logP = np.log(P) + logmaxP + logdenorm + np.log(2) + util.log_N_choose_K(V2['total_reads'][sidx], V2['var_reads'][sidx])
        logprob_models[sidx,modelidx] = logP
  return logprob_models
