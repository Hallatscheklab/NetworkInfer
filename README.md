# Computational Tools for Inferring Inter-Community Disease Transmission

This repository contains the computational tools and scripts used in our study: **“Uncovering Heterogeneous Inter-Community Disease Transmission from Neutral Allele Frequency Time Series”** by Takashi Okada, Giulio Isacchini, QinQin Yu, and Oskar Hallatschek. The study uses genomic surveillance data and mathematical models to infer transmission patterns between communities, with applications to SARS-CoV-2 data.

Introduction
---

The COVID-19 pandemic emphasized the need for precise models of disease transmission to forecast pathogen spread and inform public health interventions. Our study introduces a novel, data-driven approach to infer inter-community disease transmission rates using allele frequency convergence. The approach leverages Hidden Markov Models (HMMs) to model time series data of allele frequencies, bypassing the need for traditional, tree-based phylogenetic methods.



Methods Overview
---

The computational method relies on analyzing the convergence of allele frequencies between regions to infer cross-community transmission rates. The main computational tool is a Hidden Markov Model (HMM) integrated with an Expectation-Maximization (EM) algorithm and Markov Chain Monte Carlo (MCMC) methods.

Key Algorithms and Techniques

    •   Hidden Markov Model (HMM): Used to infer allele frequency dynamics, treating true frequencies as hidden states.
    •   Kalman Filtering: Efficiently handles noise in frequency measurements by modeling genetic drift and sampling error.
    •   Expectation-Maximization (EM) Algorithm: Accelerates the inference of model parameters.
    •   Markov Chain Monte Carlo (MCMC): Provides posterior distributions for transmission rates and genetic drift parameters.
    •   Multidimensional Scaling (MDS): Used to visualize inferred transmission networks and compare them to geographical distances.


Preparation
---
Clone this repository and prepare a conda environment 

    •   git clone https://github.com/Hallatscheklab/NetworkInfer.git 
    •   cd NetworkInfer
    •   conda env create -f environment.yml

The HMM-MCMC method is implemented in C++. A C++ compiler that supports at least C++11 is required. In the directory `HMM_MCMC/` of the local repository, compile the main.cc file

    •   g++  main.cc -std=c++1 -o NI_MCMC

You can run the program by executing

    •   ./NI_MCMC -f WFsim


Note that the HMM-MCMC method utilizes the [Eigen](http://eigen.tuxfamily.org) library, which is a header-only library located in HMM_MCMC/src/ within this repository.

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


HMM-MCMC
---
................

License
---

This project is licensed under the MIT License - see the LICENSE file for details.

Acknowledgments
---

We would like to thank the COVID-19 Genomics UK Consortium (COG-UK) and GISAID for providing access to the SARS-CoV-2 genomic data used in this study.
 
