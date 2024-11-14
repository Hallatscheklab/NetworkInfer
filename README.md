# Inferring inter-community disease transmission networks
Code for reproducing analyses in the manuscript:

Uncovering heterogenous inter-community disease transmission from neutral allele frequency time series

Takashi Okada, Giulio Isacchini, QinQin Yu, and Oskar Hallatschek

Description of Repository
---
This repository contains the computational methods we used in the above study to infer disease importation rates from allele frequency time series. Specifically, our approach uses a hidden Markov model to infer the fraction of infections a community imports from other communities based on how rapidly the allele frequencies in the focal community converge to those in the donor communities. 


HMM-EM
---

The HMM-EM method is validated using simulated data in usage_HMM_WF.ipynb.

* Dependencies
    * Numpy, Scipy, CVXOPT
    
* INPUT
    * Spatio-temporal data of allele (or lineage) counts.
    * Spatio-temporal data of the total number of sampled sequences.

* OUTPUT
    * Record of log likelihood across EM cycles.
    * Inferred importation-rate matrix ${\mathbf A}_{ij}$.
    * Inferred effective population size.
    * Least squares estimation of ${\mathbf A}_{ij}$.(noise ignored).
    * Inferred measurement noise overdispersion.

### modules/
Miscellaneous tools used for analysis.

### data/
Data from England and the USA used in the analysis.

 
