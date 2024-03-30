These are the Matlab codes to replicate the main results and tables in "Assignment Markets: Theory and Experiments". 
The data is in the online supplementary material to the paper.

-Code-

There are 3 separate stages to produce the results in the paper. Running them in this order will replicate the main tables in LaTeX format (printed to Matlab output), both for experimental subjects and simulated data points.

All files assume that YALMIP* is installed and added to Matlab path (see https://yalmip.github.io/) and a mixed-integer solver is setup. Currently all files assume gurobi is installed. For CPLEX or other solvers, please edit the sdpsettings lines in the files.

* Lofberg, J. (2004, September). YALMIP: A toolbox for modeling and optimization in MATLAB. In 2004 IEEE international conference on robotics and automation (IEEE Cat. No. 04CH37508) (pp. 284-289). IEEE.

1. main.m
This file processes the experimental data, producing the main tables (Table 2 and left sides of H.10-12) printing the results and saving them to experimentaldata.mat

2. simulations.m
This file runs the main simulations (Table F.7 and right sides of H.10-12), printing the results and saving them to powersimulations.mat

3. find_psi.m
Produces the Predictive Success Indices (Table F.8), assuming that powersimulations.mat and experimentaldata.mat were already produced by stages 1 and 2.

Remaining matlab files are function definitions used in the process.

Additionally, the bargaining solutions can also be simulated:

a. simulate_nbs.m
Additional simulations for the Nash bargaining solution.

b. simulate_pncore.m
Additional simulations for the Pairwise Nash core.

-Data-

The same experimental data is in matlab format (raw_data.mat) and in csv format (raw_data.csv). 
Every matched group of six participants is uniquely identified by a code in the total_group column, which consists of session id and group number within session.
Rounds from 1 to 15 are identified by subsession__round_number variable.
For bargaining treatments, the offers are recorded in offers_to_player and player_offer variables.
For auction treatments, the standing bids and asks (minimum prices for the sellers) at the end of the round are stored in group_asks and group_bids variables.