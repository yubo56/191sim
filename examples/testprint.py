#!/usr/bin/env python
# 02/11/16 -- example initial continuous simulation
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

if s.setAll(init, reags, rxns):
    print(s.toString())
else:
    print('Failed to init!')
