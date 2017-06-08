import argparse
import csv
import numpy as np
import scipy.stats, scipy.misc
import json
np.seterr(divide='raise')
#from numba import jit

def parse(ssmfn):
  vaf = []
  ssm_ids = []
  var_names = []
  variants = {}

  with open(ssmfn) as F:
    reader = csv.DictReader(F, delimiter='\t')
    for row in reader:
      if False and len(variants) >= 3:
        break
      variant = {
        'id': row['id'],
        'name': row['gene'],
        'ref_reads': np.array([float(V) for V in row['a'].split(',')]),
        'total_reads': np.array([float(V) for V in row['d'].split(',')]),
      }
      variant['var_reads'] = variant['total_reads'] - variant['ref_reads']
      variants[row['id']] = variant

  return variants

def generate_logprob_phi(N, models):
  prob = {}
  for M in models:
    prob[M] = np.zeros((N, N))

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
def _calc_model_prob(var1, var2, models):
  grid_step = 0.01
  N = int(1/grid_step + 1) # N: number of grid points
  G = 0.5 * np.linspace(start=0, stop=1, num=N)[:,np.newaxis] # Nx1

  S = len(var1['total_reads']) # S
  logprob_phi = generate_logprob_phi(N, models) # MxNxN
  logprob_models = np.zeros((S, len(models))) # SxM

  for s in range(S):
    for modelidx, model in enumerate(models):
      pv1, pv2 = [scipy.stats.binom.logpmf(V['var_reads'][s], V['total_reads'][s], G) for V in (var1, var2)] # Nx1
      p1 = np.tile(pv1, N)   # pv1 vector tiled as columns
      p2 = np.tile(pv2, N).T # pv1 vector tiled as rows
      P = p1 + p2 + logprob_phi[model]
      logprob_models[s,modelidx] = scipy.misc.logsumexp(P)

  logpm = np.sum(logprob_models, axis=0)
  logpm -= np.max(logpm)
  normpm = np.exp(logpm) / np.sum(np.exp(logpm))
  return normpm

#@jit
def calc_posterior(variants, models):
  N = len(variants)
  done = 0
  model_probs = {}

  for vidx1 in variants.keys():
    for vidx2 in variants.keys():
      model_probs[(vidx1,vidx2)] = _calc_model_prob(variants[vidx1], variants[vidx2], models)
      done += 1
      print(vidx1, vidx2, '%.1f%%' % (100*(done / N**2)), model_probs[(vidx1,vidx2)], sep='\t')

  return model_probs

def write_posterior(posterior, N, models, outfn):
  out = np.zeros((N, N, len(models)))
  used_pairs = set()

  for vid1, vid2 in posterior.keys():
    pair = (vid1, vid2)
    vidx1 = int(vid1[1:])
    vidx2 = int(vid2[1:])
    post = { models[idx]: posterior[(vid1,vid2)][idx] for idx in range(len(models)) }
    assert pair not in used_pairs
    used_pairs.add(pair)
    out[vidx1,vidx2] = post

  assert len(used_pairs) == N**2
  np.save(outfn, out)

def write_posterior(posterior, variants, models, outfn):
  model_probs = {}
  for midx, model in enumerate(models):
    model_probs[model] = {'%s,%s' % (vid1, vid2): P[midx] for (vid1, vid2), P in posterior.items() }

  var_names = { V: variants[V]['name'] for V in variants.keys() }
  out = {
    'models': models,
    'model_probs': model_probs,
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

  variants = parse(args.ssm_fn)
  models = ('cocluster', 'A_B', 'B_A', 'diff_branches')
  posterior = calc_posterior(variants, models)
  write_posterior(posterior, variants, models, args.out_fn)

main()
