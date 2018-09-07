import argparse
import numpy as np
import scipy.stats, scipy.misc
import json
from common import parse_ssms, Models
import concurrent.futures
from tqdm import tqdm
import itertools

#from numba import jit

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
    logprob[M][prob[M] == 0] = -100
    assert np.sum(logprob[M] == 0) == 0

  return logprob

def _calc_posterior(evidence, include_garbage, include_cocluster):
  evidence = np.copy(evidence)
  if not include_garbage:
    evidence[Models.garbage] = -np.inf
  if not include_cocluster:
    # We still want the posterior for coclustering of node `i` with `i` to be
    # certain. To signal that we're computing the (i, i) pairwise posterior,
    # the evidence for coclustering will be set to 0, and all other entries
    # will be -inf. In this instance, don't set the coclustering evidence to
    # -inf -- instead, leave it at 0 so that the posterior will indicate
    # certainty about coclustering.
    if evidence[Models.cocluster] < 0:
      evidence[Models.cocluster] = -np.inf
  B = np.max(evidence)
  evidence -= B
  posterior = np.exp(evidence) / np.sum(np.exp(evidence))
  return posterior

def _calc_model_prob(var1, var2):
  if var1['id'] == var2['id']:
    # If they're the same variant, they should cocluster with certainty.
    evidence = -np.inf * np.ones(len(Models._all))
    evidence[Models.cocluster] = 0
    return evidence

  grid_step = 0.01
  N = int(1/grid_step + 1) # N: number of grid points
  G = np.linspace(start=0, stop=1, num=N)[:,np.newaxis] # Nx1

  S = len(var1['total_reads']) # S
  logprob_phi = generate_logprob_phi(N) # MxNxN
  logprob_models = np.zeros((S, len(Models._all))) # SxM

  for s in range(S):
    for modelidx, model in enumerate(Models._all):
      if modelidx == Models.garbage:
        continue
      pv1, pv2 = [scipy.stats.binom.logpmf(V['var_reads'][s], V['total_reads'][s], np.minimum(1., (1./V['vaf_correction'])*(1 - V['mu_v'])*G)) for V in (var1, var2)] # Nx1
      p1 = np.tile(pv1, N)   # pv1 vector tiled as columns
      p2 = np.tile(pv2, N).T # pv1 vector tiled as rows
      P = p1 + p2 + logprob_phi[model]
      logprob_models[s,modelidx] = scipy.misc.logsumexp(P)

  # Garbage model
  pv1, pv2 = [scipy.stats.randint.logpmf(V['var_reads'], 0, V['var_reads'] + V['total_reads'] + 1) for V in (var1, var2)]
  logprob_models[:,Models.garbage] = pv1 + pv2

  evidence = np.sum(logprob_models, axis=0)
  return evidence

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

def _calc_model_prob(var1, var2):
  if var1['id'] == var2['id']:
    # If they're the same variant, they should cocluster with certainty.
    evidence = -np.inf * np.ones(len(Models._all))
    evidence[Models.cocluster] = 0
    return evidence

  S = len(var1['total_reads']) # S
  logprob_models = np.zeros((S, len(Models._all))) # SxM
  for s in range(S):
    for modelidx, model in enumerate(Models._all):
      if modelidx == Models.garbage:
        continue
      mcsamps = 1000000
      phi1, phi2, area = _gen_samples(modelidx, S=mcsamps)

      logP = []
      for V, phi in ((var1, phi1), (var2, phi2)):
        alpha = 2*V['var_reads'][s] + 1
        beta = np.maximum(1, V['ref_reads'][s] - V['var_reads'][s] + 1)
        logP.append(scipy.stats.beta.logpdf(phi, alpha, beta))
      logprob_models[s,modelidx] = scipy.misc.logsumexp(logP[0] + logP[1] - np.log(area)) - np.log(mcsamps)

  garb1, garb2 = [scipy.stats.randint.logpmf(V['var_reads'], 0, V['var_reads'] + V['total_reads'] + 1) for V in (var1, var2)]
  logprob_models[:,Models.garbage] = garb1 + garb2
  evidence = np.sum(logprob_models, axis=0)
  return evidence

def _calc_dummy_prob(var1, var2):
  # Just specify uniform posterior over possible relations between `var1` and
  # `var2`. This makes model computation tractable in scenarios for which we
  # have too many variants to compute this efficiently (e.g., SJETV010).
  M = len(Models._all)
  posterior = np.ones(M) / M
  evidence = np.log(posterior)
  return (posterior, evidence)

def swap_evidence(evidence):
  swapped = np.zeros(len(evidence)) + np.nan
  for midx, M in enumerate(Models._all):
    if midx in (Models.garbage, Models.cocluster, Models.diff_branches):
      swapped[midx] = evidence[midx]
    elif midx == Models.A_B:
      swapped[midx] = evidence[Models.B_A]
    elif midx == Models.B_A:
      swapped[midx] = evidence[Models.A_B]
    else:
      raise Exception('Unknown model')
  assert not np.any(np.isnan(swapped))
  return swapped

def calc_posterior(variants, use_dummy=False, parallel=1, include_garbage_in_posterior=False, include_cocluster_in_posterior=False):
  N = len(variants)
  posterior = {}
  evidence = {}

  _calc_prob = _calc_dummy_prob if use_dummy else _calc_model_prob
  #_calc_prob(variants['C24'], variants['C26'])

  combos = list(itertools.combinations(variants.keys(), 2)) + [(K,K) for K in variants.keys()]
  # Don't bother invoking all the parallelism machinery if we're only fitting a
  # single sample at a time.
  if parallel > 1:
    pairs = []
    futures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=parallel) as ex:
      for vidx1, vidx2 in combos:
        pairs.append((vidx1, vidx2))
        futures.append(ex.submit(_calc_prob, variants[vidx1], variants[vidx2]))
      for F in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc='Computing relations', unit=' pairs', dynamic_ncols=True):
        pass
    for pair, F in zip(pairs, futures):
      evidence[pair] = F.result()
  else:
    for vidx1, vidx2 in combos:
      pair = (vidx1, vidx2)
      evidence[pair] = _calc_prob(variants[vidx1], variants[vidx2])

  # Duplicate evidence's keys, since we'll be modifying dictionary.
  for A, B in list(evidence.keys()):
    if A == B:
      continue
    evidence[(B,A)] = swap_evidence(evidence[(A,B)])
  for pair in evidence.keys():
    posterior[pair] = _calc_posterior(evidence[pair], include_garbage_in_posterior, include_cocluster_in_posterior)
  return (posterior, evidence)

def generate_results(posterior, evidence, variants):
  model_probs = {}
  model_evidence = {}
  for midx, model in enumerate(Models._all):
    model_probs[model]    = {'%s,%s' % (vid1, vid2): P[midx] for (vid1, vid2), P in posterior.items() }
    model_evidence[model] = {'%s,%s' % (vid1, vid2): P[midx] for (vid1, vid2), P in evidence.items()  }

  results = {
    'models': Models._all,
    'model_probs': model_probs,
    'model_evidence': model_evidence,
    'variants': { V: {'name': variants[V]['name']} for V in variants.keys() },
  }
  return results

def write_results(results, outfn):
  with open(outfn, 'w') as F:
    json.dump(results, F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--dummy', dest='use_dummy', action='store_true',
    help='Rather than calculating posterior and evidence, use uniform to save computation time')
  parser.add_argument('sampid')
  parser.add_argument('ssm_fn')
  parser.add_argument('out_fn')
  args = parser.parse_args()

  variants = parse_ssms(args.sampid, args.ssm_fn)
  posterior, evidence = calc_posterior(variants, args.use_dummy, include_garbage_in_posterior=False, include_cocluster_in_posterior=False)
  results = generate_results (posterior, evidence, variants)
  write_results(results, args.out_fn)

if __name__ == '__main__':
  main()
