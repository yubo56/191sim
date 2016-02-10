#!/usr/bin/env python
import lib.simulation as sim

s = sim.ContinuousSim(30)
# initial state
init = [0.6, 0.4, 0]
# reagants into init vector
reags = {'X':0, 'Y': 1, 'Z' : 2}
# reactions
rxns = [(['X', 'Y'], ['Z', 'Z'], 1),
        (['X', 'Z'], ['X', 'X'], 1),
        (['Z', 'Y'], ['Y', 'Y'], 1)]

if s.setState(init) and s.setReagants(reags) and s.setReactions(rxns):
    s.run()
    s.plot()
else:
    print('Failed to init!')
