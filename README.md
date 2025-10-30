# IRG-gKSS - A Pure Hypothesis Test for Inhomogeneous Random Graph Models  

A kernelised Stein discrepancy (KSD) type test to evaluate the fit of inhomogeneous random graph (IRG) models. Implemented in R.  

## Overview  

This repository provides:  

- **`Source.R`**: Core functions to:  
  - Simulate networks from IRG models.  
  - Simulate networks with planted clique and planted hubs.  
  - Perform the IRG-gKSS test.  
- **`Example_Florentine_marriage_network.R`**: A worked example applying the test to the Florentine marriage network.  

## Prerequisites  

Ensure the following R libraries are installed:  
```r
install.packages(c("igraph", "HelpersMG", "graphkernels", "Matrix"))

## Version Control

### Software Versions
This project is tested with the following versions:

- **R**: `4.5.0`  
- **R Libraries**:
  ```r
  igraph      # v2.1.4
  HelpersMG   # v6.5
  graphkernels# v1.6.1

  Matrix      # v1.7.3
