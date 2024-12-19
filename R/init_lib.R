# List of required packages
required_packages <- c(
  "mniw", "LaplacesDemon", "matrixsampling", "matrixcalc", "invgamma", 
  "matrixStats", "loo", "Rcpp", "here", "tidyverse"
)

# Function to install and load packages if not already installed
install_and_load <- function(pkg) {
  if (!library(pkg, logical.return = T)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Install and load the required packages
lapply(required_packages, install_and_load)

# library(mniw)
# library(LaplacesDemon) # dmatrixnorm()
# library(matrixsampling) # rinvwishart()
# library(matrixcalc) # is.positive.definite()
# library(invgamma) # rinvgamma
# library(matrixStats) # logSumExp()
# library(loo)
# library(Rcpp)
# library(here)
# library(tidyverse)
# library(here)
# options(tidyverse.quiet = TRUE)
# library(matrixNormal) # dmatnorm() may have underflow
# library(microbenchmark)
# library(MNIW)
# library(fst)

source(here("R/func_gen_data.R"))
source(here("R/func_MNIW_sampler.R"))
source(here("R/func_FF.R"))
source(here("R/func_BS.R"))
source(here("R/func_FFBS.R"))
source(here("R/func_generate_grid_Dan.R"))
sourceCpp(here("src/FF.cpp"))
sourceCpp(here("src/BS.cpp"))
sourceCpp(here("src/rmn.cpp"))