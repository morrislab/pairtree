import argparse
import numpy as np
import scipy.stats, scipy.misc
import json
from common import parse_ssms, Models

np.seterr(divide='raise')
np.set_printoptions(threshold=np.nan, linewidth=120)

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
    logprob[M][prob[M] == 0] = -100
    logprob[M][prob[M] != 0] = np.log(prob[M][prob[M] != 0])
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
      pv1, pv2 = [scipy.stats.binom.logpmf(V['var_reads'][s], V['total_reads'][s], (1 - V['mu_v'])*G) for V in (var1, var2)] # Nx1
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
  normpm = np.exp(logpm) / np.sum(np.exp(logpm))
  return (normpm, logpm + B)

#@jit
def calc_posterior(variants):
  N = len(variants)
  done = 0
  posterior = {}
  evidence = {}

  for vidx1 in variants.keys():
    for vidx2 in variants.keys():
      posterior[(vidx1,vidx2)], evidence[(vidx1,vidx2)] = _calc_model_prob(variants[vidx1], variants[vidx2])
      done += 1
      print(vidx1, vidx2, '%.1f%%' % (100*(done / N**2)), posterior[(vidx1,vidx2)], sep='\t')

  return (posterior, evidence)

def write_posterior(posterior, evidence, variants, outfn):
  model_probs = {}
  model_evidence = {}
  for midx, model in enumerate(Models._all):
    model_probs[model]    = {'%s,%s' % (vid1, vid2): P[midx] for (vid1, vid2), P in posterior.items() }
    model_evidence[model] = {'%s,%s' % (vid1, vid2): P[midx] for (vid1, vid2), P in evidence.items()  }

  var_names = { V: variants[V]['name'] for V in variants.keys() }
  out = {
    'models': Models._all,
    'model_probs': model_probs,
    'model_evidence': model_evidence,
    'var_names': var_names,
  }
  with open(outfn, 'w') as F:
    json.dump(out, F)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('ssm_fn')
  parser.add_argument('out_fn')
  args = parser.parse_args()

  variants = parse_ssms(args.ssm_fn)
  posterior, evidence = calc_posterior(variants)
  write_posterior(posterior, evidence, variants, args.out_fn)

main()
