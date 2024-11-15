# HMM-MCMC

The posterior distributions of the coupling matrix $A_{ij}$ and the effective population sizes $N_{{\rm eff},i}$ at each deme are computed by using MCMC (assuming 
flat priors).


# Test
In terminal, go to the folder HMM_MCMC, which is the same directory as main.cpp.
- Compile main.cc for example, by entering *g++  main.cc -std=c++11 -O3*
- Then, enter *./a.out -f demo -g res_demo -m 100 -d example -b 0 -n 0 -D DB*
(Minimally, the code can be run by *./a.out -f demo  -d example*).

If you can compile & run the code, an example data (HMM_MCMC/input/example/) is used, and the output will be created in output/example.

# Input files

- *counts_demo.csv*:
    A file of time-series count data for Nlins lineages, Ndeme demes, T timepoints. The input csv file is a T*Nlins-by-Ndemes matrix, which is created by arranging Nlins T-by-Ndemes matrices vertically.

- *shape_demo.csv*:
    A file of 3 lines, where T, Nlins, Ndeme are written.

# Output files
- *A_res_demo.csv* shows the ensembles of the coupling matrix A. For example, for Ndeme=3, each row of this file shows A00,A01,A02,A10,A11,A12,A20,A21,A22 at a ce*rtain Monte Carlo time. 
- *Ne_res_demo.csv* shows the ensembles of the effective population sizes Ne. For example, for Ndeme=3, each row of this file shows Ne_1, Ne_2, Ne_3 at a certain Monte Carlo time.
The number of recorded A and Ne in these files is set to numprint=1000 in main.cpp. 

- *logLH_res_demo.csv* records the value of log likelihood as a function of Monte Carlo time.
- *logfile_res_demo.csv* records the options specified by -D and -n and also records the accepted percentage as a function of Monte Carlo time.


# Options 
- *-m 10000*: Total Monte Carlo steps. The default value 2000 is small. Using at least 10000 is recommended. 100000 is sufficient in most cases.  

- *-f demo*: filename of the input. -f demo is used if the input files are "counts_demo.csv" and "shape_demo.csv".

- *-g res_demo*: filename of the output. If -g demo is used, A_demo.csv and Ne_demo.csv are created. The default value is 'inferred'.

- *-d dir*: The name of the folder that contains the input files. The output files are crated in output/dir/.

- *-n 0* (default), *-n 1* or *-n 2* : Noise assumption, specified by 0 or 1 or 2. The proportionality constants of the diagonal components of the covariance matrices are;
    - -n 0: $f_{t,i,l}(1-f_{t,i,l})$
    - -n 1: $f_l(1-f_l)$
    -  -n 2: $f_{i,l}(1-f_{i,l})$

, where $f_{t,i,l}$ is the observed frequency of lineage l in deme i  at time t. $f_l$ and $f_{i,l}$ are frequencies averaged over t&i and over t, respectively. 


- *-D nonDB* or *-D DB*: nonDB is used for  a general row-stochastic matrix ($\sum_j A_{ij} =1$). DB is used if the detailed balance is imposed. The default value is nonDB.

- *-b 0*: Fraction of burn-in (between 0 and 1). By default, b=0. You should throw away the initial part. *logLH_res_demo.csv* 


-  *-N 1000*: [Minor option] The value of the effective population size $N_{\rm eff}$ at the initial Monte Carlo step. The default value is N=1000 (for all demes). Usually you do not need to modify this value. 


# Example for simulated data
See usage_HMM_WF.ipynb.


# References
- Bishop, Christopher M., and Nasser M. Nasrabadi. Pattern recognition and machine learning. Vol. 4. No. 4. New York: springer, 2006.

- Noé, Frank. "Probability distributions of molecular observables computed from Markov models." The Journal of chemical physics 128.24 (2008): 244103.