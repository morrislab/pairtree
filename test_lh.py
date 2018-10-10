import numpy as np
import lh
import util
import common

def softmax(V):
  return np.exp(V) / np.sum(np.exp(V))

def main():
  np.set_printoptions(linewidth=400, precision=3, threshold=np.nan, suppress=True)
  np.seterr(divide='raise', invalid='raise')

  V1 = {'var_reads': np.array([1702]), 'total_reads': np.array([4069])}
  V2 = {'var_reads': np.array([2500]), 'total_reads': np.array([19100])}
  for V in (V1, V2):
    V['var_reads'] = (V['var_reads'] / 10).astype(np.int)
    V['total_reads'] = (V['total_reads'] / 10).astype(np.int)
    V['mu_v'] = 0.5
    V['ref_reads'] = V['total_reads'] - V['var_reads']
    V['vaf'] = V['var_reads'].astype(np.float) / V['total_reads']

  S = 1
  for V, var, total in ((V1, [500], [1000]), (V2, [100], [1000])):
    V['var_reads'] = np.array(S * var)
    V['total_reads'] = np.array(S * total)
    V['ref_reads'] = V['total_reads'] - V['var_reads']
    V['vaf'] = V['var_reads'] / V['total_reads']

  for M in (
    lh.calc_lh_binom_quad,
    lh.calc_lh_binom_mc_1D,
    lh.calc_lh_binom_mc_2D,
    lh.calc_lh_binom_grid,
  ):
    M_name = M.__name__
    M = util.time_exec(M)
    evidence_per_sample = M(V1, V2)
    evidence_per_sample[:,common.Models.garbage] = -np.inf
    evidence = np.sum(evidence_per_sample, axis=0)
    print(
      M_name,
      '%.3f ms' % util.time_exec._ms,
      evidence,
      softmax(evidence),
      sep='\t',
    )

main()
