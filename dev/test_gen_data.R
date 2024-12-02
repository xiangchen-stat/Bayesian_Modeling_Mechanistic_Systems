setwd("D:/Documents/UCLA/0-Administrative/GSR/MNIW/dev")
source("../R/gen_data.R")

seed <- 1234
set.seed(seed)
N <- 100
S <- 100
p <- 4
nT <- 20

dt <- gen_ffbs_data(N = N, S = S, p = p, nT = nT)
