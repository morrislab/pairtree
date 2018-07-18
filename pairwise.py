import argparse
import numpy as np
import scipy.stats, scipy.misc
import json
from common import parse_ssms, Models
import concurrent.futures
from tqdm import tqdm

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

#@jit
def _calc_model_prob(var1, var2):
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

  logpm = np.sum(logprob_models, axis=0)
  B = np.max(logpm)
  logpm -= B
  posterior = np.exp(logpm) / np.sum(np.exp(logpm))
  evidence = logpm + B
  return (posterior, evidence)

def _calc_dummy_prob(var1, var2):
  # Just specify uniform posterior over possible relations between `var1` and
  # `var2`. This makes model computation tractable in scenarios for which we
  # have too many variants to compute this efficiently (e.g., SJETV010).
  M = len(Models._all)
  posterior = np.ones(M) / M
  evidence = np.log(posterior)
  return (posterior, evidence)

#@jit
def calc_posterior(variants, use_dummy=False, parallel=1):
  N = len(variants)
  posterior = {}
  evidence = {}

  _calc_prob = _calc_dummy_prob if use_dummy else _calc_model_prob

  pairs = []
  futures = []
  with concurrent.futures.ProcessPoolExecutor(max_workers=parallel) as ex:
    for vidx1 in sorted(variants.keys()):
      for vidx2 in sorted(variants.keys()):
        pairs.append((vidx1, vidx2))
        futures.append(ex.submit(_calc_prob, variants[vidx1], variants[vidx2]))
    for F in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc='Computing relations', unit=' pairs'):
      pass
  for pair, F in zip(pairs, futures):
    posterior[pair], evidence[pair] = F.result()

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
  posterior, evidence = calc_posterior(variants, args.use_dummy)
  results = generate_results (posterior, evidence, variants)
  write_results(results, args.out_fn)

if __name__ == '__main__':
  main()
