source("../R/gen_data.R")
source("../R/BS.R")
source("../R/pds.R")
source("../R/MNIW_sampler.R")
source("naiive_BS.R")
source("test_utils.R")
source("naiive_FF.R")
source("naiive_FFBS.R")
library(MNIW)
library(matrixsampling)
library(microbenchmark)
library(Rcpp)
sourceCpp("../src/FF.cpp")
sourceCpp("../src/BS.cpp")


seed <- 1234
set.seed(seed)
N <- 100
S <- 100
p <- 2
nT <- 10
nsam <- 10
Sys.time()
dt <- gen_ffbs_data(N = N, S = S, p = p, nT = nT)

## FFBS----
Sys.time()
res_ffbs <- FFBS_R_naiive(nsam = nsam, Y = dt$Y, F_ls = dt$para$F0, G_ls = dt$para$G0, 
                      W_ls = dt$para$W0, V_ls = dt$para$V0, 
                      m0 = dt$para$m0, M0 = dt$para$M0, 
                      n0 = dt$para$n0, D0 = dt$para$D0, 
                      nT = nT, delta = 1.0)
Sys.time()

# 1k by 1k, T=10, 20s to generate
# nsam = 10, 1min40s to FFBS
# temp <- res_ffbs$T0[,,1]
# for (i in 2:nsam) {
#   temp <- temp + res_ffbs$T0[,,i]
# }
# temp <- temp/nsam
# norm(temp - dt$para$m0, "F")
