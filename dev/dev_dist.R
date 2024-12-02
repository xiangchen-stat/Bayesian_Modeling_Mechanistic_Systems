source("../R/gen_data.R")
library(matrixsampling)
library(magrittr)
seed <- 1234
set.seed(seed)
N <- 100
S <- 100
p <- 2
nT <- 20
dt <- gen_ffbs_data(N = N, S = S, p = p, nT = nT)

x <- runif(S)
y <- runif(S)
loc <- cbind(x, y) %>% as.matrix()
dist <- stats::dist(loc, method = "euclidean", diag = T, upper = T) %>% as.matrix()

