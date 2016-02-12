#!/usr/bin/env python
# 02/12/16 --   tests a system with multiple allowed branching reactions,
#               two-decay possibilities.
#               Also shows off new  syntax parsing for reagants, reactions
from lib.simulation import *
if __name__ == '__main__':
    stoch = StochasticSim(2)
    # same init, reagants, reactions for both!
    init = [49, 0, 0]
    reagants = ['X', 'R', 'C']
    reactions = [
            'X -0.1-> R',
            'X -0.01-> C'
            ]
    if stoch.setState(init) and stoch.setReagants(reagants) and\
            stoch.setReactions(reactions):
        stoch.run(printTraj=True)
        stoch.plot()
