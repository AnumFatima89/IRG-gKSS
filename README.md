# IRG-gKSS  
### A Pure Hypothesis Test for Inhomogeneous Random Graph Models  

This repository contains the implementation of the method proposed in the paper:  

**“A Pure Hypothesis Test for Inhomogeneous Random Graph Models Based on a Kernelised Stein Discrepancy.”**  

The code is implemented primarily in **R**, with an optional Python-based kernel implementation.

---

## Overview  

This repository provides:

- **`source-code.R`**  
  Core implementation including:
  - Simulation of inhomogeneous random graphs (IRGs)  
  - Simulation of networks with planted cliques and planted hubs  
  - Implementation of the IRG-gKSS test  

- **`source-code_graphlet_py.R`**  
  Extension of the method using **graph kernels from the Python package GraKeL** (via `reticulate`).

- **`Experiments/`**  
  A directory containing all R scripts required to reproduce the experimental results presented in the paper.

---

## Prerequisites  

### R packages  

Install the required R packages:

```r
install.packages(c(
  "igraph", "netrankr", "network", "graphkernels", "randnet",
  "arrangements", "igraphdata", "dplyr", "ggplot2", "ggpubr",
  "doParallel", "foreach", "qgraph", "sjmisc", "boot",
  "devtools", "remotes", "ggrepel", "digest", "Matrix", "tictoc"
))
```

## Python (optional, for GraKeL version)

If using `source-code_graphlet_py.R`, install:

- Python ≥ 3.x  
- `grakel` package  

and configure `reticulate` to point to your Python environment:

```r
Sys.setenv(RETICULATE_PYTHON = "path/to/your/python")
```

## Basic example

```r
source("source-code.R")

# Example group labels
C <- c(1, 1, 1, 2, 2, 2)

# Example probability matrix
p <- matrix(c(
  0.2, 0.05,
  0.05, 0.15
), nrow = 2, byrow = TRUE)

# Simulate a graph
G <- sample_IRG(p)

# Run IRG-gKSS test
result <- GOF_IRG(test_G = G, C = C, p_0 = p)

print(result$P_value)
```

## Version Information

This project has been tested with:

- **R**: 4.5.0  

### Key package versions

```r
igraph        # 2.1.4
HelpersMG     # 6.5
graphkernels  # 1.6.1
Matrix        # 1.7.3
reticulate    # 1.42.0
grakel        # 0.1.8
```

## Spectral Test Extensions

The spectral test extensions (utils.R) included in this repository are adapted from publicly available code:

https://www.dropbox.com/scl/fi/evx0mvcgmgfeixm9ubhz2/gof_test_sbm.zip

Originally developed for:

> Jing Lei (2016). *A goodness-of-fit test for stochastic block models.*  
> The Annals of Statistics, 44(1): 401–424.

The implementation here includes minor modifications to facilitate comparison with the IRG-gKSS framework.

## Data

This repository includes external datasets used for reproducibility of the experiments.

- **Lazega Lawyers dataset**  
  Sourced from: https://www.stats.ox.ac.uk/~snijders/siena/Lazega_lawyers_data.htm  

These datasets are **not owned** by the authors of this repository.  
Please refer to the corresponding data directory for full attribution and usage information.
