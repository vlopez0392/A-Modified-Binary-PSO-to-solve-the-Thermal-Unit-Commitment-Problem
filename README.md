# A-Modified-Binary-PSO-to-solve-the-Thermal-Unit-Commitment-Problem

This repository contains the results of my Master thesis as published in here:
https://etd.lis.nsysu.edu.tw/ETD-db/ETD-search/view_etd?URN=etd-0713118-191039

This project was written in the following language: MATLAB 

This project contains the folllowing files:

__## Scripts__ : 
- __Swarm_Generator.m__: Generates a swarm of binary particlesof dimensions N x T. This script is useful for any swarm-based metaheuristic that aims to tackle the classic Thermal Unit Scheduling problem.

- __Swarm Optimizer.m__: Implements a Binary Particle Swarm Optimization (BPSO) algorithm to solve the Thermal Unit Commitment Problem. Eight modified activations functions based on the popular sigmoid and tanh functions are considered. 


__## Helper Functions__:
- __AFLC.m__: Returns a sorted array in a descending fashion according to AFLC criteria.

- __ED_fmincon.m__: Finds the minimum of the constrained non-linear optimization objective (fmincon)

- __F_LIM_ED.m__: Enhanced Lambda Iteration (ELI) Method to solve the Economic Dispatch problem (EDP)

- __F_LIM_ED_RR.m__ : Same as F_LIM_ED.m but considers the ramping rate constraints in the TUCP problem. 

- __LIM_ED.m__ : Lambda Iteration Method to solve the Economic Dispatch problem (EDP)

- __LIM_ED_RR.m__ Same as F_LIM_ED.m but considers the ramping rate constraints in the TUCP problem. 

- __Recomm_swp.m__: Recommit units to satisfy System Power Demand at time t

- __SU_COST.m__: Computes the start up and shut down costs (to expand) of the UC schedule 

- __check_MUT_MDT.m__: Checks if the current schedule satisfies the MUT/MDT constraints are satisfied in the UC problem 

- __check_SR_PD.m__: Checks if the current schedule satisfies the SR constraint is satisfied in the UC problem 

- __constraint_repair.m__ : Schedule repair function based on the proposed pivot heuristic algorithm 

- __count_intervals.m__ :Returns a vector of ON/OFF intervals & their respective duration. 

- __mod_repair_MDT_MUT.m__: Modified straightforward repair for pivot heuristic

- __sf_repair_MDT_MUT.m__: Straightforward repair of MDT/MUT constraint in the Thermal UC problem 

- __until_zero__ : Returns the first or last occurring zero closest to the pivot.This function is part of the modified repair strategy by pivot heuristic



