source("../R/gen_data.R")
library(mniw)
library(LaplacesDemon)
library(MixMatrix)
library(matrixNormal)
library(magrittr)
library(microbenchmark)

seed <- 1234
set.seed(seed)
N <- 1
S <- 1
p <- 2
niter <- 10
Yt <- gen_ran_matrix(N, S)
M <- gen_pd_matrix(N)
Sigma <- gen_pd_matrix(S)
RM = chol(M)
RSigma = chol(Sigma)


mniw::dMNorm(X = Yt, Lambda = Yt,
             SigmaR = M, SigmaC = Sigma, log = TRUE)
LaplacesDemon::dmatrixnorm(X = Yt, M = Yt,
                           U = M, V = Sigma, log = TRUE)
matrixNormal::dmatnorm(A = Yt, M = Yt,
                       U = M, V = Sigma, log = TRUE)
MixMatrix::dmatrixnorm(x = Yt, mean = Yt,
                       L = M, R = Sigma, log = TRUE)

log(dnorm(Yt, Yt, sqrt(M * Sigma)))



# single sample
microbenchmark(mniw::dMNorm(X = Yt, Lambda = Yt,
                            SigmaR = M, SigmaC = Sigma, log=TRUE),
               LaplacesDemon::dmatrixnorm(X = Yt, M = Yt,
                                          U = M, V = Sigma, log = TRUE),
               matrixNormal::dmatnorm(A = Yt, M = Yt,
                                      U = M, V = Sigma, log = TRUE),
               MixMatrix::dmatrixnorm(x = Yt, mean = Yt,
                                      L = M, R = Sigma, log = TRUE)
)



