# Generate matrix normal distribution with upper triangular cholesky factor
# m: mean
# RM: upper triangular cholesky factor of M
# RSigma: upper triangular cholesky factor of Sigma
rmn_chol <- function(m, RM, RSigma){
  p <- nrow(m)
  S <- ncol(m)
  vec <- rnorm(p*S)
  Z <- matrix(vec, nrow = p, ncol = S)
  theta <- m + crossprod(RM, Z) %*% RSigma
  return(theta)
}

# Generate matrix t distribution
# m: mean
rmt <- function(m, M, nu, D){
  Sigma <- rinvwishart(n = 1, nu = nu, Omega = D)[,,1]
  RM <- chol(M)
  RSigma <- chol(Sigma)
  theta <- rmn_chol(m, RM, RSigma)
  return(theta)
}

# packages to compare with
library(matrixsampling)

# generate random pd matrix
gen_pd_matrix <- function(dim){
  tempS <- matrix(rnorm(dim * dim), nrow  = dim, ncol  = dim)
  S <- tempS %*% t(tempS)
  return(S)
}

# generate random matrix
gen_ran_matrix <- function(nrow, ncol){
  S <- matrix(rnorm(nrow * ncol), nrow  = nrow, ncol  = ncol)
  return(S)
}

# Compare Matrix Normal with package
p <- 50
S <- 60
m <- gen_ran_matrix(p, S)
M <- gen_pd_matrix(p)
Sigma <- gen_pd_matrix(S)
RM <- chol(M)
RSigma <- chol(Sigma)
y1 <- rmn_chol(m, RM, RSigma)
y2 <- rmatrixnormal(1, m, M, Sigma)[,,1]
qqplot(y1, y2)
abline(0,1)

# Compare Matrix t with package
nu <- S + 2
D <- gen_pd_matrix(S)
y3 <- rmt(m, M, nu+p-1, D)
y4 <- rmatrixt(1, nu = nu , M = m, U = M, V = D)[,,1]
qqplot(sort(y3), sort(y4))
abline(0,1)

nsam <- 20
y5 <- c()
y6 <- c()
for (i in 1:nsam) {
  y5 <- c(y5, sort(rmt(m, M, nu+p-1, D)))
  y6 <- c(y6, sort(rmatrixt(1, nu = nu , M = m, U = M, V = D)[,,1]))
}
qqplot(sort(y5), sort(y6))
abline(0,1)
