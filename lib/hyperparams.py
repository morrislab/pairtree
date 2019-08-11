explanations = {
  'gamma': '''
    Proportion of tree modifications that should use mutrel-informed choice for
    node to move, rather than uniform choice
  ''',

  'zeta': '''
    Proportion of tree modifications that should use mutrel-informed choice for
    destination to move node to, rather than uniform choice
  ''',

  'iota': '''
    Probability of initializing with mutrel-informed tree rather than fully
    branching tree when beginning chain
  '''
}

defaults = {
  'gamma': 0.7,
  'zeta': 0.7,
  'iota': 0.7,
}

assert set(explanations.keys()) == set(defaults.keys())
