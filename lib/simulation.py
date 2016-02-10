#! /usr/bin/env python
# class defining a simulation 
# Changelog:
#   02/07/16 -- Set up base Simulation abstract class

from abc import ABCMeta, abstractmethod                 # for abstract classes
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

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
        self.state = state
        self.reagants = reagants 
        self.revreag = reagants 
        self.reactions = reactions
        self.ignoreCase = ignoreCase

        self.traj = None
        self.times = None
        self.setState(state)
        self.setReagants(reagants)
        self.setReactions(reactions)

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
        if len(reagants.values()) == len(self.state) and\
            list(set(reagants.values())) == list(range(len(self.state))) and\
            len(set(reagants.keys())) == len(reagants.keys()):
            # right length if len(values) = len(self.state), unique/right values if
            # set = range(len(self.state)). Then check uniqueness for keys
            self.reagants = reagants
            for i in reagants.keys():
                self.revreag[reagants[i]] = i
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
            if len(rxn) != 3 or not (isinstance(rxn[2], float) or\
                    isinstance(rxn[2], int)) or rxn[2] < 0 or\
                    not all([i in self.reagants for i in rxn[0]]) or\
                    not all([i in self.reagants for i in rxn[1]]):
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

    def plot(self, FN=None, title='Simulation Trajectory', xlabel='Time (s)',
            ylabel='Concentration (M)', loc='best'):
        """
        plot results if have
        :param str FN: to save file
        :param str title: title of graph
        :param str xlabel: xlabel of graph
        :param str ylabel: ylabel of graph
        :param str loc: legend location
        :return: handle to plot, Nonetype if no trajectory to plot
        """
        plt.style.use('ggplot')
        if self.traj is None:
            return None
        else: # plot traj
            for i, traj in enumerate(self.traj):
                plt.plot(self.times, traj, label=self.revreag[i])
            plt.axis([self.times.min(), self.times.max(), 0.95 *
                self.traj.min(), 1.05 * self.traj.max()])
            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.legend(loc=loc)
            if FN is None:
                plt.show()
            else:
                plt.savefig(FN)
            plt.clf()


    @abstractmethod
    def setState(self, state):
        """
        asserts state is list, implementations should set state
        implementations should type-check for list values (e.g. int for
        discrete, nonnegative float for continuous)
        :param list state: input state vector
        :return: True if passed type checks and set, else fails assert
        """
        assert(isinstance(state, list))
        return True

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
    scipy.integrate.odeint to perform numerical integration of mass action
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
        try:
            if super(ContinuousSim, self).setState(state) and\
                    all([i >= 0 for i in state]):
                # checks super (assert is list) and nonnegative float
                self.state = [float(i) for i in state]
                return True
            else:
                print('hi')
                return False
        except ValueError:
            return False

    def run(self):
        """
        Runs the continuous-time simulation using odeint
        """
        maxk = max([rxn[2] * len(rxn[0]) * np.sqrt(sum(self.state))
            for rxn in self.reactions]) # heuristic for 1/max timescale
        self.times = np.arange(0, self.tf, 1/(10 * maxk))
        self.traj = np.transpose(odeint(self._ds, self.state, self.times,
            args=(self.reagants, self.reactions)))

    @staticmethod
    def _ds(s, t, reagants, reactions):
        """
        helper function to compute derivatives
        """
        ds = [0] * len(reagants)   # differential state
        for rxn in reactions:
            rate = rxn[2]
            for x in rxn[0]:            # inputs, calc rate
                rate *= s[reagants[x]]
            for x in rxn[0]:            # subtracts from ds (in's)
                ds[reagants[x]] -= rate
            for x in rxn[1]:            # adds to ds (out's)
                ds[reagants[x]] += rate
        return ds
