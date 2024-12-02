source("../R/gen_data.R")
# library(matrixsampling)
# library(matrixNormal)
# library(MBSP)
library(tidyverse)
library(microbenchmark)
# library(Rcpp)
# sourceCpp("../src/rmn.cpp")

seed <- 1234
set.seed(seed)
N <- 100 / 20
S <- 100 / 50
p <- 2
niter <- 10
m <- gen_ran_matrix(N, S)

# single sample
microbenchmark(cbind(m, m),
               bind_cols(m, m))
