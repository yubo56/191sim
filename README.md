# 191sim
Python implementation of simulations for Caltech class BE/CS/CNS/Bi 191a. 

Note: To run examples in `examples/`, run `python -m examples/<file>` or `cp` to root directory of projcet.

Dependencies:
* Python 3.5.1
* Python-Scipy 0.16.1-2
* Python-Numpy 1.10.4-1

List of current features:
* Continuous CRN
* Stochastic CRN
* matplotlib.pyplot Plotting to both file and screen
    * Print Trajectory
* Constructor takes either list/list(string) for reagants/reactions instead of dict/list(tuple).
    * Syntax parser for reactions
* Simulation setups can be toString() and fromFile()
    * Simulation results are not written to file to save space.
* Reversible Reactions
* Subsample for plotting
    * Python plotting very slow for large lists.

List of planned features (in rough order of priority, not necessarily difficulty):
* Enumerator for DSDs
* Better error messages when failure, more edge cases
* Interactive Front End
    * Automatically create reagants and set to 0 in add\_reactions?
    * Need incremental mutators rather than just setBLAH
* Loading/Saving simulation object to file
* Adaptive Stepsize

Changelog:
* 02/29/16 -- Reversibile reactions, subsampled plotting
* 02/27/16 -- Removed some unnecessary dict lookups
    * Load from file, print string (can save to file)
* 02/11/16 -- Set up syntax parser for both reactions and reagants.
    * Set up stochastic Simulation child class
* 02/09/16 -- Set up continuous Simulation child class
* 02/07/16 -- Set up base Simulation abstract class
