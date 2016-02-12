#!/usr/bin/env python
# 02/12/16 -- example initial stochastic simulation
import lib.simulation as sim

s = sim.StochasticSim(30)
# initial state
init = [1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0]
# reagants into init vector
reags = {
        'C1':0, 
        'C2':1, 
        'C3':2, 
        'C4':3, 
        'X1-0': 4,
        'X1-1': 5,
        'X2-0': 6,
        'X2-1': 7,
        'X3-0': 8,
        'X3-1': 9,
        'X4-0': 10,
        'X4-1': 11
        }
# reactions
rxns = [
        (['C1', 'X1-0'], ['C1', 'X1-1'], 0.3),
        (['C1', 'X1-1'], ['C2', 'X1-0'], 0.3),
        (['C2', 'X2-0'], ['C1', 'X2-1'], 0.3),
        (['C2', 'X2-1'], ['C3', 'X2-0'], 0.3),
        (['C3', 'X3-0'], ['C1', 'X3-1'], 0.3),
        (['C3', 'X3-1'], ['C4', 'X3-0'], 0.3),
        (['C4', 'X4-0'], ['C1', 'X4-1'], 0.3),
        (['C4', 'X4-1'], ['C1', 'X4-0'], 0.3)
        ]

if s.setState(init) and s.setReagants(reags) and s.setReactions(rxns):
    s.run(printTraj=True)
    s.plot()
else:
    print('Failed to init!')
