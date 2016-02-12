# 191sim
Python implementation of simulations for Caltech class BE/CS191a. 

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

List of planned features (in rough order of priority, not necessarily difficulty):
* Implement ignoreCase feature (already have flag)
* Interactive Front End
    * Automatically create reagants and set to 0 in add\_reactions?
    * Need incremental mutators rather than just setBLAH
* Loading/Saving simulation object to file
