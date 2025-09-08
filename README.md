# EMbru-code
Code for the implementation of the EMbru method, presented in the article “EMbru: A quick and accurate Bayesian inference method for Hawkes point process modelling”. This repository includes scripts to reproduce the results of the simulation study, as well as code to simulate data from a Hawkes process using the acceptance–rejection method. The scripts are written in R and organised as follows:

- "simulate_data_thinning.R" : This script simulates data from a Hawkes process using the acceptance–rejection method. It generates spatio-temporal events based on specified triggering functions.

- "synthetic_data.RData" : This file contains a simulated dataset generated using a Gaussian triggering function for the spatial component, and an exponential triggering function for the temporal component.

- "run_EMbru.R" : This is the main script for estimating the parameters of the Hawkes process using the Bayesian EMbru approach. It involves two stages:
  
  ⁃ Stage 1 (EM method): implemented in "code_em_stage1.R", which includes the necessary functions for the Expectation-Maximisation procedure. 
  
  ⁃ Stage 2 (inlabru): implemented in "code_inla_stage2.R", which provides the functions required to perform inference with inlabru, based on the output from Stage 1. 
  
- "compute_hessian_standard_errors.R" contains the code to compute the Hessian and standard errors from both stages of EMbru. To run it, the objects resultado_em from stage 1 of EMbru and fit_inlabru from stage 2 of EMbru are required (They can be requested from the author at bernabeu@uji.es).
