# Computational Tools for Inferring Inter-Community Disease Transmission

This repository contains the computational tools and scripts used in our study: **“Uncovering Heterogeneous Inter-Community Disease Transmission from Neutral Allele Frequency Time Series”** by Takashi Okada, Giulio Isacchini, QinQin Yu, and Oskar Hallatschek. The study uses genomic surveillance data and mathematical models to infer transmission patterns between communities, with applications to SARS-CoV-2 data.

## Introduction


The COVID-19 pandemic emphasized the need for precise models of disease transmission to forecast pathogen spread and inform public health interventions. Our study introduces a novel, data-driven approach to infer inter-community disease transmission rates using allele frequency convergence. The approach leverages Hidden Markov Models (HMMs) to model time series data of allele frequencies, bypassing the need for traditional, tree-based phylogenetic methods.



## Methods Overview


The computational method relies on analyzing the convergence of allele frequencies between regions to infer cross-community transmission rates. The main computational tool is a Hidden Markov Model (HMM) integrated with an Expectation-Maximization (EM) algorithm and Markov Chain Monte Carlo (MCMC) methods.

### Key Algorithms and Techniques

- **Hidden Markov Model (HMM)**: Used to infer allele frequency dynamics, treating true frequencies as hidden states.
- **Kalman Filtering**: Efficiently handles noise in frequency measurements by modeling genetic drift and sampling error.
- **Expectation-Maximization (EM) Algorithm**: Accelerates the inference of model parameters.
- **Markov Chain Monte Carlo (MCMC)**: Provides posterior distributions for transmission rates and genetic drift parameters.
- **Multidimensional Scaling (MDS)**: Used to visualize inferred transmission networks and compare them to geographical distances.


## Preparation

The **HMM-EM method** is implemented in Python (dependencies: Numpy, Scipy, CVXOPT). To prepare a Conda environment, clone this repository and execute the following commands:
```
git clone https://github.com/Hallatscheklab/NetworkInfer.git 
cd NetworkInfer
conda env create -f environment.yml
```
You can test the program by executing the following command in the root directory:
```
python HMMEM.py WFsim
```

The **HMM-MCMC method** is implemented in C++. A C++ compiler that supports at least C++11 is required. In the directory `HMM_MCMC/` of the local repository, compile the main.cc file
```
g++  main.cc -std=c++1 -o NI_MCMC
```
You can test the program by executing the following command in the directory `HMM_MCMC/`:
```
 ./NI_MCMC -f WFsim
```
<!-- The flag `-f` specifies the input filename (here, WFsim) and loads the input data located in `HMM_MCMC/input`. For more information on the list of flags, see below or refer to `HMM_MCMC/README_HMMMCMC.md`. -->

**Note**: The HMM-MCMC method utilizes the [Eigen](http://eigen.tuxfamily.org) library, which is a header-only library located in `HMM_MCMC/src/` within this repository. 

## HMM-EM
 
### INPUT
Located in `HMM_EM/input/`:
- `counts_filename.npy`: Spatio-temporal data of allele (or lineage) counts. A numpy array with shape ($T$, $N_{\rm allele}$, $N_{\rm deme}$), where $T$ denotes the number of timepoints. 
- `totcounts_filename.npy`: Spatio-temporal data of the total number of sampled sequences.  A numpy array with shape ($T$, $N_{\rm deme}$). 

### OUTPUT
Located in `HMM_EM/output/`:
- `LH_ridge0.0_filename.npy`: Record of log likelihood across EM cycles.
- `A_ridge0.0_filename.npy`: Inferred importation-rate matrix ${\mathbf{A}_{ij}}$.
- `Ne_ridge0.0_filename.npy`: Inferred effective population size.
- `Csn_ridge0.0_filename.npy`: Inferred measurement noise overdispersion.
- `A_LS_filename.npy`: Least squares estimation of ${\mathbf{A}_{ij}}$ (noise ignored).

Here `ridge0.0` in the output files denotes the strength of regularization, represented by $c_{\rm ridge}$ (see usage below).

### Usage

In the root directory, execute the following command:
```
python HMMEM.py filename 
```
For L2 regularizaiton, add the flag `--ridge 0.01` (here, $c_{\rm ridge}=0.01$; the default value is $c_{\rm ridge}=0.0$). When a nonzero $c_{\rm ridge}$ is specified, a numpy array $\Lambda_{ij}$ is loaded from `HMM_EM/input/ridgemat_filename.npy` and the term $c_{\rm ridge} \sum_{i,j}\Lambda_{ij} A_{ij}^2$ is introduce to the cost function. 


The HMM-EM method is validated using simulated data in `usage_HMM_WF.ipynb`.
   

## HMM-MCMC

    
### INPUT
Located in `HMM_MCMC/input/`:
- `counts_filename.csv`: Spatio-temporal data of allele (or lineage) counts. This CSV file has dimension ($T \times N_{\rm allele}$,  $N_{\rm deme}$). The count data needs to be concatenated along alleles, starting with the counts for the 1st allele, followed by the 2nd allele, and so on, such that the $T\times (i-1)+t$-th row of the csv file represents the allele frequencies of the $i$-th allele at the $t$-th timepoint across $N_{\rm deme}$ populations. 

- `totcounts_filename.csv`: Spatio-temporal data of the total number of sampled sequences.  This csv file has dimensions ($T $,  $N_{\rm deme}$)
- `shape_filename.csv` contains the dimensions $(T,  N_{\rm allele}$,  $N_{\rm deme})$, which are used to convert the 2-dim data of allele counts into 3-dimensional format in the program. 

### OUTPUT
Located in `HMM_MCMC/output/`:
- `A_filename.csv`: Inferred importation-rate matrix ${\mathbf A}_{ij}$.
- `Ne_filename.csv`: Inferred effective population size.
- `C_filename.csv`: Inferred measurement noise overdispersion.
- `logLH_filename.csv`: Recoded values of log likelihood along a sequence of MCMC 
- `logfile_filename.csv`: Log file

### Usage

In the directory `HMM_MCMC/`, execute the following command:
```
 ./NI_MCMC -f filename -m 100000 -D nonDB -N 1000 
```
The flags, other than `-f filename`, are optional. The `-m` flag specifies the total Monte Carlo steps. The `-D` flag accepts either nonDB or DB, which indicates whether detailed balance is imposed in the inference of Aij. The `-N` flag specifies the initial value of the effective population size. For other flags, refer to `HMM_MCMC/README_HMMMCMC.md`.

The HMM-MCMC method is validated using simulated data in `usage_HMM_WF.ipynb`.

---

### `modules/`
Miscellaneous tools used for analysis.

### `data/`
Data from England and the USA used in the analysis.

License
---

This project is licensed under the MIT License - see the LICENSE file for details.

Acknowledgments
---

We would like to thank the COVID-19 Genomics UK Consortium (COG-UK) and GISAID for providing access to the SARS-CoV-2 genomic data used in this study.
 
