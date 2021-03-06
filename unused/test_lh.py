import numpy as np
import lh
import util
import common
import sys

def create_vars():
  variants = {
    #'V1': {'var_reads': [18], 'total_reads': [100]},
    #'V2': {'var_reads': [10], 'total_reads': [100]},

    'V1': {'var_reads': [500], 'total_reads': [1000]},
    'V2': {'var_reads': [100], 'total_reads': [1000]},

    #'V1': {'var_reads': [1702], 'total_reads': [4069]},
    #'V2': {'var_reads': [2500], 'total_reads': [19100]},

    #'V1': {'var_reads': [0], 'total_reads': [200]},
    #'V2': {'var_reads': [179], 'total_reads': [356]},
  }
  print(sorted(variants.items()))

  S = 1
  for vid, V in variants.items():
    for K in ('var_reads', 'total_reads'):
      V[K] = np.array(S*V[K]).astype(np.int)
    V['id'] = vid
    V['ref_reads'] = V['total_reads'] - V['var_reads']
    V['vaf'] = V['var_reads'].astype(np.float) / V['total_reads']
    V['omega_v'] = np.array(S*[0.5])

  variants = {vid: common.convert_variant_dict_to_tuple(V) for vid, V in variants.items()}
  return (variants['V1'], variants['V2'])

def main():
  np.set_printoptions(linewidth=400, precision=3, threshold=sys.maxsize, suppress=True)
  np.seterr(divide='raise', invalid='raise')
  V1, V2 = create_vars()

  estimators = (
    # lh.calc_lh_quad will be slower than usual on its first invocation due to
    # Numba JIT compilation. Don't be alarmed by seemingly poor runtime from it
    # as a result.
    lh.calc_lh_quad,
    lh.calc_lh_mc_1D,
    lh.calc_lh_mc_2D,
    lh.calc_lh_mc_2D_dumb,
    lh.calc_lh_grid,
  )
  max_estimator_len = max([len(M.__name__) for M in estimators])

  for M in estimators:
    M_name = M.__name__
    M = util.time_exec(M)
    evidence_per_sample = M(V1, V2)
    evidence_per_sample[:,common.Models.garbage] = lh.calc_garbage(V1, V2)
    evidence = np.sum(evidence_per_sample, axis=0)
    print(
      M_name.ljust(max_estimator_len),
      '%.3f ms' % util.time_exec._ms,
      evidence,
      util.softmax(evidence),
      sep='\t',
    )

main()
