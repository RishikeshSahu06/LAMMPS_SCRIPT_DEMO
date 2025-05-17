## How to run 

First generate the required polymer using $polymer.py$

First use run_simulations.py to run various simulations parallely for different seeds

Dumps are stored in $/dumps$ folder 

Using that we use $production.lammpstrj$ to compute required parameters using $analyze.cpp$ and $analyze.sh$ and plot the results using $analyzePlot.py$