#! /usr/bin/env python
# class defining a simulation 
# Changelog:
#   02/07/16 -- Set up base Simulation abstract class

from abc import ABCMeta, abstractmethod                 # for abstract classes

class Simulation(object):
    """
    Abstract class, represents a general simulaiton

    Idiomatically, should first setState, then setReagants, then setReactions
    """
    __metaclass__ = ABCMeta

    def __init__(self, tf,
            state=list(), 
            reagants=dict(), 
            reactions=list(),
            ignoreCase=True):
        """
        :ivar float tf: final time, length of simulation
        :ivar list state: list of quantity of each reagant. indicies are stored
            in dict reagants. Default: empty list()
        :ivar dict reagants: dict of reagants for the simulation, mapping
            strr reagant name to int index in state. Default: empty dict()
        :ivar list reactions: list of reactions for the simulation, implemented
            as list of 3-tuples, each tuple a pair of lists of in/out reagants
            and a reaction rate
            Default: empty list()
        :ivar bool ignoreCase: whether to ignore case in reagants/reactants.
            Default: True
        :ivar list traj: list of states, simulation trajectory
        :ivar list times: list of times for traj
        """
        self.tf = tf
        self.setState(state)
        self.setReagants(reagants)
        self.setReactions(reactions)
        self.ignoreCase = ignoreCase
        self.traj = list()
        self.times = list()

    def setVal(self, reagant, value):
        """
        Sets state of reagant to value. Should never need to be used
        :param str reagant: name of reagant
        :param float value: value to set reagant to
        :return: True if set, False if KeyError on reagants dict (not a valid
            reagant name)
        """
        try:
            self.state[self.reagants[reagant]] = value
            return True
        except KeyError:
            return False

    def setReagants(self, reagants):
        """
        sets the reagants dict
        type-checks for unique list of keys, values and checks range of values
        is equal to length of state vector
        :param dict reagants: dict of reagants
        :return: True if set, False if error
        """
        assert(isinstance(reagants, dict))
        if len(reagants.values()) == len(state)\
                and list(set(reagants.values())) == list(range(len(state)))\
                and list(set(keys)) == list(keys):
            # right length if len(values) = len(state), unique/right values if
            # set = range(len(state)). Then check uniqueness for keys
            self.reagants = reagants
            return True
        else:
            return False

    def setReactions(self, reactions):
        """
        set reactions list.
        Checks formatting of list whether all reagants are in reagants dict and
        whether reaction rates are >= 0
        :param list reactions: list of 3-tuples, each list() of reagants, list()
        of reagants and float reaction rate >= 0
        :return: True if set, False if error
        """
        assert(isinstance(reactions, list))
        for rxn in reactions:
            if len(rxn) != 3 or not instanceof(rxn[2], float) or rxn[2] < 0\
                    or not all([i in self.reagants for i in [rxn[0], rxn[1]]]):
                # checking length of tuple and reaction rate and whether
                # reagants in self.reagants
                return False
        self.reactions=reactions
        return True

    def setLen(self, tf):
        """
        Sets the termination time for the simulation
        Checks if nonnegative float
        :param float tf: final time
        :return: True if set, False if not
        """
        assert(isinstance(tf, float))
        if tf < 0:
            return False
        else:
            self.tf = tf
            return True

    def getResults(self):
        """
        Fetches simulation trajectory
        """
        return self.traj, self.times

    @abstractmethod
    def setState(self, state):
        """
        sets the state vector to input list
        asserts state is list
        implementations should type-check for list values (e.g. int for
        discrete, nonnegative float for continuous)
        :param list state: input state vector
        :return: True if passed type checks and set, False otherwise
        """
        assert(isinstance(state, list))

    @abstractmethod
    def run(self):
        """
        Should just run simulation and store results into self.traj, self.times
        Should be designed such that can modify self.tf and rerun
        """
        pass

class ContinuousSim(Simulation):
    """
    Represents a continuous-time simulation of a CRN, uses
    scipy.integrate.numint to perform numerical integration of mass action
    equations
    """
    def __init__(self, tf,
            state=list(), 
            reagants=dict(), 
            reactions=list(),
            ignoreCase=True):
        """
        Everything handled in super constructor Simulation()
        """
        super(ContinuousSim, self).__init__(tf, state, reagants, reactions,
                ignoreCase)

    def setState(self, state):
        """
        sets the state vector to input list
        asserts state is list, checks nonnegative floats and sets
        :param list state: input state vector
        :return: True if passed type checks and set, False otherwise
        """
        if super(ContinuousSim, self).setSate(state) and\
                all([isinstance(i, float) and i >= 0 for i in state]):
            # checks super (assert is list) and nonnegative float
            return True
        return False

    def run(self):
        """
        Runs the continuous-time simulation using numint
        """
        pass
