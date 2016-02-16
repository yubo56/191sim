#! /usr/bin/env python
# class defining a simulation 
# Changelog:
#   02/11/16 -- Set up stochastic Simulation child class
#   02/09/16 -- Set up continuous Simulation child class
#   02/07/16 -- Set up base Simulation abstract class

from abc import ABCMeta, abstractmethod                 # for abstract classes
from scipy.integrate import odeint                      # ODE integrator
import random                                           # random.random()
import numpy as np
import matplotlib.pyplot as plt
import bisect                                           # binary search

def toRxn(rxn):
    """
    Turns a string reaction into a ([], [], k) form
    String format is "A + B + C -k-> D + E + F$
    """
    try:
        rxn = rxn.strip().split('-')
        k = float(rxn[1])
        ins = [i.strip() for i in rxn[0].split('+')]
        outs = [i.strip() for i in rxn[2][1: ].split('+')]
        return (ins, outs, k)
    except ValueError:
        raise ValueError

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
            02/11/16 -- added support for passing in list of reagants
        :ivar list reactions: list of reactions for the simulation, implemented
            as list of 3-tuples, each tuple a pair of lists of in/out reagants
            and a reaction rate
            Default: empty list()
        :ivar bool ignoreCase: whether to ignore case in reagants/reactants.
            Default: True
        :ivar list traj: list of states, simulation trajectory
        :ivar list times: list of times for traj
        """
        self._tf = tf
        self._state = list()
        self._reagants = dict()
        self._revreag = dict()
        self._reactions = list()
        self._ignoreCase = ignoreCase

        self._traj = None
        self._times = None
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
            self._state[self._reagants[reagant]] = value
            return True
        except KeyError:
            return False

    def setReagants(self, reagants):
        """
        sets the reagants dict from list or dict
        type-checks for unique list of keys, values and checks range of values
            is equal to length of state vector, or checks length
        :param dict reagants: dict of reagants
        :return: True if set, False if error
        """
        if isinstance(reagants, dict):
            if len(reagants.values()) == len(self._state) and\
                list(set(reagants.values())) == list(range(len(self._state))) and\
                len(set(reagants.keys())) == len(reagants.keys()):
                # right length if len(values) = len(self._state), unique/right values if
                # set = range(len(self._state)). Then check uniqueness for keys
                self._reagants = reagants
                for i in reagants.keys():
                    self._revreag[reagants[i]] = i
                return True
            else:
                return False
        elif isinstance(reagants, list):
            if len(reagants) == len(self._state):
                for i, j in enumerate(reagants):
                    self._reagants[j] = i
                    self._revreag[i] = j
                return True
            else:
                return False
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
        if isinstance(reactions, list):
            for index, rxn in enumerate(reactions):
                if isinstance(rxn, str):
                    try:
                        rxn = toRxn(rxn) # just convert it and typecheck below
                        reactions[index] = rxn  # update it for final step
                    except ValueError:
                        return False
                if len(rxn) != 3 or not (isinstance(rxn[2], float) or\
                        isinstance(rxn[2], int)) or rxn[2] < 0 or\
                        not all(map(lambda i : i in self._reagants, rxn[0])) or\
                        not all(map(lambda i : i in self._reagants, rxn[1])):
                    # checking length of tuple and reaction rate and whether
                    # reagants in self._reagants
                    return False
            self._reactions=reactions
            return True
        else:
            return False

    def setAll(self, states, reagants, reactions):
        return all([self.setState(states), self.setReagants(reagants),
                self.setReactions(reactions)])

    def setLen(self, tf):
        """
        Sets the termination time for the simulation
        Checks if nonnegative float
        :param float tf: final time
        :return: True if set, False if not
        """
        try:
            if tf < 0:
                return False
            else:
                self._tf = tf
                return True
        except TypeError: # usually when tf is not number
            return False

    def getTraj(self):
        """
        Fetches simulation trajectory
        """
        return self._times, self._traj

    def printTraj(self):
        """
        Prints Trajectory and times in pretty print"""
        for i, j in zip(self._times, np.transpose(self._traj)):
            print("%s\t%s" % (str(round(i,3)), str(j)))

    def plot(self, FN=None, title='Simulation Trajectory', xlabel='Time (s)',
            ylabel='Concentration (M)', loc='best', select=None):
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
        if self._traj is None:
            return None
        else: # plot traj
            maxrange = self._traj.min()
            minrange = self._traj.max()
            if select is not None:
                try:
                    for i in [self._reagants[sel] for sel in select]:
                        plt.plot(np.append(self._times, self._tf), 
                                np.append(self._traj[i], self._traj[i][-1]),
                                label=self._revreag[i])
                        maxrange = max(self._traj[i].tolist() + [maxrange])
                        minrange = min(self._traj[i].tolist() + [minrange])
                        print(sel + ': ' + str(self._traj[i][-1]))
                except KeyError:
                    print("Invalid Select! Returning")
                    return None
            else:
                maxrange = self._traj.max()
                minrange = self._traj.min()
                for i, traj in enumerate(self._traj):
                    plt.plot(np.append(self._times, self._tf), 
                            np.append(traj, traj[-1]), label=self._revreag[i])
            plt.axis([0, self._tf, 0.95 *
                minrange, 1.05 * maxrange])
            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.legend(loc=loc)
            print("Plotted " + str(len(self._times)) + " points!")
            if FN is None:
                plt.show()
            else:
                plt.savefig(FN)
            plt.clf()


    @abstractmethod
    def setState(self, state):
        """
        True if  state is list, implementations should set state
        implementations should type-check for list values (e.g. int for
        discrete, nonnegative float for continuous)
        :param list state: input state vector
        :return: True if passed type checks and set, else false
        """
        return isinstance(state, list)

    @abstractmethod
    def run(self):
        """
        Should just run simulation and store results into self._traj, self._times
        Should be designed such that can modify self._tf and rerun
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
        state is list, checks nonnegative number (castable to float) and
        sets state vector
        :param list state: input state vector
        :return: True if passed type checks and set, False otherwise
        """
        try:
            if super(ContinuousSim, self).setState(state) and\
                    all(map(lambda i:i >= 0, state)):
                # checks super (True if list) and nonnegative float
                self._state = list(map(lambda i: float(i), state))
                return True
            else:
                return False
        except ValueError:
            return False

    def run(self, printTraj=False):
        """
        Runs the continuous-time simulation using odeint
        """
        maxk = max(map(lambda rxn: rxn[2] * len(rxn[0]) * 
            np.sqrt(sum(self._state)), self._reactions)) # heuristic timescale
        self._times = np.arange(0, self._tf, 1/(10 * maxk))
        self._traj = np.transpose(odeint(self._ds, self._state, self._times,
            args=(self._reagants, self._reactions)))
        if printTraj == True:
            self.printTraj()

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

class StochasticSim(Simulation):
    """
    Represents a stochastic CRN, uses first principles to calculate propensities
    and choose among allowed reactions. No external packages are used to evolve
    forward in time
    """
    def __init__(self, tf,
            state=list(), 
            reagants=dict(), 
            reactions=list(),
            ignoreCase=True):
        """
        Everything handled in super constructor Simulation()
        """
        super(StochasticSim, self).__init__(tf, state, reagants, reactions,
                ignoreCase)

    def setState(self, state):
        """
        sets the state vector to input list
        checks state is list, checks nonnegative *ints* and sets
        :param list state: input state vector
        :return: True if passed type checks and set, False otherwise
        """
        try:
            if super(StochasticSim, self).setState(state) and\
                    all(map(lambda i:i >= 0, state)) and\
                    all(map(lambda i:isinstance(i, int), state)):
                # checks super (True if list) and nonnegative float
                self._state = list(state)
                return True
            else:
                return False
        except ValueError:
            return False

    def run(self, printTraj=False):
        """
        Runs stochastic CRN until tf
        """
        t = 0
        state = list(self._state)
        traj = list()
        traj.append(list(state))
        times = [0]
        while t < self._tf:
            dt, rxn = self._nextReaction(state)
            t += dt
            if t < self._tf:
                state = self._doReaction(state, rxn)
                traj.append(list(state))
                times.append(t)
        # no off by one since last reaction doesn't happen
        self._traj = np.transpose(traj)
        self._times = np.array(times)
        if printTraj == True:
            self.printTraj()

    def plot(self, FN=None, title='Simulation Trajectory', xlabel='Time (s)',
            ylabel='Number of Molecules', loc='best', select=None):
        super(StochasticSim, self).plot(FN=FN, title=title, xlabel=xlabel,
                ylabel=ylabel, loc=loc, select=select)

    def _nextReaction(self, state):
        """
        Perform a reaction and return time elapsed
        calculates using propensities
        :return float, reaction: next time, reaction computed
            returns float('inf'), None if no reaction possible
        """
        props = [(0,None)]              # partial sums of propensities
        for rxn in self._reactions:
            prop = rxn[2]
            for i in rxn[0]:            # multiply by counts of inputs
                prop *= state[self._reagants[i]]
            if prop != 0:               # append if nonzero
                props.append((props[-1][0] + prop, rxn))
        r = random.random() * props[-1][0]
        index = bisect.bisect(props, (r, None)) # Need tuple to work
                                                # in range (0, len - 1)
        if props[-1][0] == 0:
            return float('inf'), None
        else:
            return random.expovariate(props[-1][0]),\
                    props[index][1] # time taken

    def _doReaction(self, state, rxn):
        """
        helper function to change a state by a reaction
        """
        for i in rxn[0]:
            state[self._reagants[i]] -= 1
        for i in rxn[1]:
            state[self._reagants[i]] += 1
        return state
