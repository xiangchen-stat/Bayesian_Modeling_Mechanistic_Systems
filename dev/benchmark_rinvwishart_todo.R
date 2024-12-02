source("../R/gen_data.R")
library(matrixsampling)
library(matrixNormal)
library(MBSP)
library(magrittr)
library(microbenchmark)
library(Rcpp)
sourceCpp("../src/rmn.cpp")

seed <- 1234
set.seed(seed)
N <- 100 / 4
S <- 100 / 2
p <- 2
niter <- 10
m <- gen_ran_matrix(N, S)
M <- gen_pd_matrix(N)
Sigma <- gen_pd_matrix(S)
RM = chol(M)
RSigma = chol(Sigma)

matrixsampling::rinvwishart()
LaplacesDemon::rinvwishart()
mniw::rwishart(inverse = T)

a1 <- rmatrixnormal(1, M = m, U = M, V = Sigma, checkSymmetry = F, keep = F)
a2 <- rmatnorm(1, M = m, U = M, V = Sigma)
a3 <- matrix_normal(M = m, U = M, V = Sigma)
a4 <- rmn_chol(m = m, RM = chol(M), RSigma = chol(Sigma))

# single sample
microbenchmark(rmatrixnormal(1, M = m, U = M, V = Sigma, checkSymmetry = F, keep = F),
               rmatnorm(1, M = m, U = M, V = Sigma),
               matrix_normal(M = m, U = M, V = Sigma),
               rmn_chol(m = m, RM = chol(M), RSigma = chol(Sigma)))

# single sample cpp
microbenchmark(rmn_chol(m = m, RM = chol(M), RSigma = chol(Sigma)),
               rmn_chol_cpp(m = m, RM = chol(M), RSigma = chol(Sigma)),
               rmn_chol_cpp(m = m, RM = RM, RSigma = RSigma),
               rmn_cpp(m = m, U = M, V = Sigma))

nsam <- 100
res1 <- array(dim = c(N, S, nsam))
res2 <- array(dim = c(N, S, nsam))
res3 <- array(dim = c(N, S, nsam))
for (i in 1:nsam) {
  res1[,,i] <- rmn_chol(m = m, RM = chol(M), RSigma = chol(Sigma))
  res2[,,i] <- rmn_chol_cpp(m = m, RM = chol(M), RSigma = chol(Sigma))
  res3[,,i] <- rmn_cpp(m = m, U = M, V = Sigma)
}
res4 <- rmatrixnormal(n = nsam, M = m, U = M, V = Sigma)

qqplot(x = as.vector(res1[,,20]), y = as.vector(res2[,,20]))
abline(0, 1, col = "red")

qqplot(x = as.vector(res2[,,20]), y = as.vector(res3[,,20]))
abline(0, 1, col = "red")

qqplot(x = as.vector(res4[,,20]), y = as.vector(res1[,,20]))
abline(0, 1, col = "red")

qqplot(x = as.vector(res4[,,20]), y = as.vector(res2[,,20]))
abline(0, 1, col = "red")

vec1 <- as.vector(res1)
vec2 <- as.vector(res2)
vec3 <- as.vector(res3)
vec4 <- as.vector(res4)
qqplot(x = vec1, y = vec2)
abline(0, 1, col = "red")

qqplot(x = vec2, y = vec3)
abline(0, 1, col = "red")

qqplot(x = vec4, y = vec1)
abline(0, 1, col = "red")

qqplot(x = vec4, y = vec2)
abline(0, 1, col = "red")

# multiple sample
microbenchmark(rmatrixnormal(50, M = m, U = M, V = Sigma, checkSymmetry = F, keep = F),
               rmn_chol_more(nsam = 50, m = m, RM = chol(M), RSigma = chol(Sigma)),
               {for (i in 1:50) {
                 rmn_cpp(m = m, U = M, V = Sigma)
               }},
               {for (i in 1:50) {
                 rmn_chol_cpp(m = m, RM = RM, RSigma = RSigma)
               }})
