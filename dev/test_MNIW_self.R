# Check if pacman package is installed
if (!requireNamespace("pacman", quietly = TRUE)) {
  # If not installed, install pacman
  install.packages("pacman")
}

library(pacman)
p_load(matrixsampling)
p_load(magrittr)
p_load(microbenchmark)
p_load(RcppEigen)
p_load(Rcpp)

# library(matrixsampling)
# library(magrittr)
# library(microbenchmark)
# library(RcppEigen)
# library(Rcpp)

source("naiive_MNIW.R")
sourceCpp("../src/MNIW.cpp")

set.seed(1234)
n <- 10
q <- 5
p <- 5
v <- q + 1
tempS <- matrix(rnorm(q*q), nrow  = q, ncol  = q)
S <- tempS %*% t(tempS)
C <- matrix(rnorm(p*q), nrow  = p, ncol  = q)
tempVb <- matrix(rnorm(p*p), nrow  = p, ncol  = p)
Vb <- tempVb %*% t(tempVb)
X <- matrix(rnorm(n*p), nrow  = n, ncol  = p)
sigmapre <- rinvwishart(1, nu = v, Omega = S, epsilon = 0, checkSymmetry = TRUE)
sigma <- matrix(sigmapre, ncol = q, nrow = q)
Bpre <- rmatrixnormal(1, M = C, U = Vb, V = sigma, checkSymmetry = TRUE, keep = TRUE)
B <- matrix(Bpre, nrow = p, ncol = q)
H <- toeplitz(n:1)
zero <- matrix(0, n, q)
Epre <- rmatrixnormal(1, zero, H, sigma)
E <- matrix(Epre, nrow = n, ncol = q)
Y <- X %*% B + E

res1 <- MNIW_R(X = X, Y = Y, H = H, v = v, S = S, C = C, Vb = Vb)
res2 <- MNIW_cpp(X = X, Y = Y, H = H, v = v, S = S, C = C, Vb = Vb)

# assessment
cat("Degrees of freedom difference:", res1$vnew - res2$vnew ,"\n")
cat("S Froebenius norm difference:", norm(res1$Snew - res2$Snew, "F"),"\n")
print(res1$Snew - res2$Snew)

# Ditto theta
cat("C Froebenius norm difference:", norm(res1$Cnew - res2$Cnew, "F"), "\n")
cat("Vb Froebenius norm difference:", norm(res1$Vbnew - res2$Vbnew, "F"), "\n")


microbenchmark(MNIW_R(X = X, Y = Y, H = H, v = v, S = S, C = C, Vb = Vb),
               MNIW_cpp(X = X, Y = Y, H = H, v = v, S = S, C = C, Vb = Vb))
