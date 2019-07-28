explanations = {
  'rho': '''
    Weight of mutrel fit term when selecting node to move within tree, such
    that we prefer nodes with high mutrel error
  ''',

  'tau': '''
    Weight of depth term when selecting node to move within tree, such that we  prefer nodes deeper in tree
  ''',

  'theta': '''
    Weight of ancestral pairwise probabilities when determining potential
    parent probability distribution for selected node while initializing tree,
    such that nodes with high ancestral probability are preferred as parents
  ''',

  'kappa': '''
    Weight of tree depth when determining potential parent probability
    distribution for selected node while initializing tree, such that nodes
    deeper in the existing tree are preferred as parents
  ''',

  'gamma': '''
    Proportion of tree modifications that should use uniform rather than
    mutrle-informed perturbation
  ''',
}

defaults = {
  'rho': 4.,
  'tau': 1.,
  'theta': 4.,
  'kappa': 1.,
  'gamma': 0.2,
}

rho = None
tau = None
theta = None
kappa = None
gamma = None
