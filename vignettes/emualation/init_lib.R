# library(MNIW)
# library(fst)
library(mniw)
library(LaplacesDemon) # dmatrixnorm()
library(matrixsampling) # rinvwishart()
library(matrixcalc) # is.positive.definite()
library(invgamma) # rinvgamma
library(matrixStats) # logSumExp()
library(loo)
# library(microbenchmark)
library(Rcpp)
library(here)
library(tidyverse)
options(tidyverse.quiet = TRUE)
# library(matrixNormal) # dmatnorm() may have underflow
library(here)
source(here("R/func_gen_data.R"))
source(here("R/func_MNIW_sampler.R"))
source(here("R/func_FF.R"))
source(here("R/func_BS.R"))
source(here("R/func_FFBS.R"))
source(here("R/func_generate_grid_Dan.R"))
sourceCpp(here("src/FF.cpp"))
sourceCpp(here("src/BS.cpp"))
sourceCpp(here("src/rmn.cpp"))