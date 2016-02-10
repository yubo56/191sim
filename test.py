#!/usr/bin/env python
import lib.simulation as sim

s = sim.ContinuousSim(2)
init = [8, 0]                       # initial state
reags = {'X':0, 'A': 1}             # reagants/indicies into state vector
rxns = [(['X'], ['X', 'A'], 0.5),     # ([in], [out], rate) tuples
        (['A', 'A'], ['A'], 0.5)]

if s.setState(init) and s.setReagants(reags) and s.setReactions(rxns):
    s.run()
    s.plot()
else:
    print('Failed to init!')
