import time
import numpy as np
import lh

def time_exec(f):
  def wrap(*args):
    time1 = time.time()
    ret = f(*args)
    time2 = time.time()
    ms = (time2-time1)*1000.0
    time_exec._ms = ms
    return ret
  return wrap
time_exec._ms = None

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
  for V, var, total in ((V1, [40], [100]), (V2, [10], [100])):
    V['var_reads'] = np.array(S * var) / 10
    V['total_reads'] = np.array(S * total) / 10
    V['ref_reads'] = V['total_reads'] - V['var_reads']
    V['vaf'] = V['var_reads'] / V['total_reads']

  for M in (
    lh.calc_lh_binom_quad,
    lh.calc_lh_binom_mc_1D,
    lh.calc_lh_binom_mc_2D,
    lh.calc_lh_binom_grid,
  ):
    M_name = M.__name__
    M = time_exec(M)
    model_prob = M(V1, V2)
    print(M_name, '%.3f ms' % time_exec._ms, np.sum(model_prob, axis=0), sep='\t')

main()
