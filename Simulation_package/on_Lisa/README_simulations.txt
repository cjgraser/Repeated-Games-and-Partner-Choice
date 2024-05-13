This directory contains all code that was used to run the simulations on the LISA computing cluster (https://www.surf.nl/en/lisa-computing-cluster-extra-computing-power-for-research).

the "basic.sh" script submits multiple simulations as a slurm job. 
the "analyze.sh" script performs some data cleaning steps and provides some readouts that require less storage than the complete simulation results

all other scripts are sourced in from these two. 
The main functionality of the simulations is in the "all_functions_5d.jl" file

