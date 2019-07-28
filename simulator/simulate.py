import pairwise
import numpy as np
import scipy.stats

def create_var(phihat, S, T):
  V = {
    'total_reads': np.array(S*[T]),
    'vaf_correction': 1,
    'mu_v': 0.499
  }
  V['var_reads'] = scipy.stats.binom.rvs(n=V['total_reads'], p=(V['mu_v'] * phihat))
  # Stupid SciPy makes this a scalar when S = 1.
  if S == 1:
    V['var_reads'] = np.array([V['var_reads']])
  return V

def summarize(vec):
  vec = np.array(vec)
  print('mean=%s\tstdev=%s' % (np.mean(vec, axis=0), np.std(vec, axis=0)))

def main():
  np.set_printoptions(suppress=True, threshold=np.nan, linewidth=300, precision=3)
  params = {
    'trials': 100,
    'phihat1': 0.9,
    'delta': 0.1,
    'T': 200,
    'S_lower': 1,
    'S_upper': 10,
  }
  params['phihat2'] = params['phihat1'] - params['delta']
  print(params)
  print()

  for S in range(params['S_lower'], params['S_upper'] + 1):
    posterior = []
    evidence = []
    for _ in range(params['trials']):
      V1 = create_var(params['phihat1'], S=S, T=params['T'])
      V2 = create_var(params['phihat2'], S=S, T=params['T'])
      P, E = pairwise._calc_model_prob(V1, V2)
      posterior.append(P)
      evidence.append(E)
    print('S=%s' % S)
    summarize(posterior)
    summarize(evidence)

main()
