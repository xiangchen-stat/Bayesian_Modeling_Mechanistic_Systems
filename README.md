# Bayesian Modeling of Mechanistic Systems

This repository contains the code to reproduce the results presented in the paper *Dynamic Bayesian Learning for Spatiotemporal Mechanistic Models*.

### Repository Structure

The repository is organized as follows:

- `R/`: Contains R code implementing the methods.
- `src/`: Contains C/C++ code, optimized using **Rcpp** and **RcppEigen** for performance improvements.
- `data/`: Contains example data files used in paper.
- `figures/`: Contains figures generated by vignettes.

### Dependencies

Make sure that all the required libraries used in the vignettes are installed by running `Rscript install_dependencies.R`. You will be able to run the vignettes without modifying it.

### Run vignettes

To reproduce the results from the application section of the paper, go to the subdirectories within the `vignettes` folder. Each subdirectory contains a script that generates all the results for the corresponding section of the paper, whether via `Rscript [vignette_file.R]` or opening the R shell and running `source("[vignette_file.R]")`, where [vignette_file.R] is substituted for the appropriate vignette file. (e.g. vignettes/emulation/emulate_EP_MNIW.R for emulation plots.) The output will be saved in the `figures` directory.
