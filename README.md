# Repeated-Games-and-Partner-Choice

This repository contains Julia scripts we used for the simulations in (*insert line with link to paper*) in the "simulation_package" directory. Moreover, it contains some additional functionality to numerically compute stationary distributions of the process (1) for arbitrary finite strategy sets in the limit of rare mutations ("Leaving/Explicit_Markov_Chain"), and (2) for the four-strategy model discussed in the paper for arbitrary mutation rates ("Leaving/Four_strategy_model_Simplices"). 
FSA are represented as matrices. 
The first row represents the output (0,1,2) --> (Leave, Cooperate, Defect)
With errors, rows 2 to 4 represent the out-arrows for (CC CD DC DD). Without errors only two additional rows are needed, so rows 2 to 3 correspond to the out-arrows for (,C and ,D). 

To pre-allocate storage space, one can choose a maximal number of states for the automata (max_n). (We choose max_n so that for our choices of probabilities of adding 
vs deleting states, and an initial automaton with only a single state, max_n is never reached.)
To reduce required memory usage, we store strategies as quintuplets (or triplets) of a 3-adic and 2 (4) max_n-adic numbers.

The "data_type_back_and_forths.jl" script contains all conversions between these strategy representations.

The "Eq_Check_Module" contains the code for the best responder algorithm, and wrappers that use this algorithm to assess whether a population state is in equilibrium
