# Check if pacman package is installed
if (!requireNamespace("pacman", quietly = TRUE)) {
  # If not installed, install pacman
  install.packages("pacman")
}

library(pacman)
p_load(matrixsampling)
p_load(invgamma)
p_load(magrittr)
p_load(microbenchmark)
p_load(RcppEigen)
p_load(Rcpp)

# library(matrixsampling)
# library(magrittr)
# library(microbenchmark)
# library(RcppEigen)
# library(Rcpp)

# setwd("C:/Users/xiang/OneDrive/Desktop/MNIW-eigen/test")
sourceCpp("../src/eigen.cpp")
source("../src/MNIW-Dan.R") # pure R code
source("../src/sampler.R")
# 
# dyn.load("MNIW.dll")
# source("../src/MNIW.R")

set.seed(1234)
n <- 30
q <- 25
p <- 25
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

res1 <- MNIW_cpp(X = X, Y = Y, H = H, v = v, S = S, C = C, Vb = Vb)

## Dan---
# compute Sigma | Y
nu0 <- v
Psi0 <- S
V <- H
m0 <- C
M0 <- Vb
Sigma_givY <- Sigma_post(nu0, Psi0, Y, X, V, m0, M0)
ESigma <- Sigma_givY$Psi/(Sigma_givY$nu - q - 1)

# compute Theta | Y, Sigma
Theta_givYSigma <- Theta_post(m0, M0, ESigma, Y, X, V)

res3 <- list()
res3$vnew <- Sigma_givY$nu
res3$Snew <- Sigma_givY$Psi
res3$Cnew <- Theta_givYSigma$m
res3$Vbnew <- Theta_givYSigma$M

# assessment
cat("Degrees of freedom difference:", res1$vnew - res3$vnew ,"\n")
print(res1$Snew - res3$Snew)
cat("S Froebenius norm difference:", norm(res1$Snew - res3$Snew, "F"),"\n")
# Ditto theta
print(res1$Snew - res3$Snew)
cat("C Froebenius norm difference:", norm(res1$Cnew - res3$Cnew, "F"), "\n")
print(res1$Snew - res3$Snew)
cat("Vb Froebenius norm difference:", norm(res1$Vbnew - res3$Vbnew, "F"), "\n")

## sampling ---
n_sam <- 20
samp <- MNIW_sampler(n_sam = n_sam, X = X, H = H, 
              v = res1$vnew, S = res1$Snew, C = res1$Cnew, Vb = res1$Vbnew, seed = 123)

## Daniel
Psiinv_post <- chol(Sigma_givY$Psi)
Psiinv_post <- chol2inv(Psiinv_post)

# Sample inverse wishart
Sigmainv_samples <- NULL

set.seed(123)
for (i in seq(1, n_sam)) {
  Sigmainv_samples <- rbind(Sigmainv_samples, c(rWishart(1, Sigma_givY$nu, Psiinv_post)))
}


# Sample rmn NITER times for Theta and put it into a data frame
Theta_samples <- NULL
u_M <- chol(Theta_givYSigma$M)
for (i in seq(1, n_sam)) {
  Theta_samples <- rbind(Theta_samples, c(rmn(Theta_givYSigma$m, Theta_givYSigma$M, ESigma)))
}


# a <- as.vector(samp$B[,,1:10])
# b <- as.vector(samp$B[,,11:20])
# qqplot(a,b)
# abline(0,1)
# 
# a <- as.vector(Sigmainv_samples)
# b <- as.vector(samp$sigma)
# qqplot(a,b)
# abline(0,1)
# Sigmainv_samples <- Sigmainv_samples


jpeg("../out/qqplot_B.jpg", width = 1200, height = 1200)
par(mfrow=c(3,3))
par(cex=1.5) # Scales the size of points globally
par(cex.axis=0.8, cex.lab=1, cex.main=1) # Adjust font sizes for axis text, labels, and main title
for (i in 1:9) {
  B1 <- as.vector(samp$B[,,i])
  B2 <- as.vector(Theta_samples[i,])
  qqplot(B1, B2)
  abline(0,1)
}
dev.off()

jpeg("../out/qqplot_sigma.jpg", width = 1200, height = 1200)
par(mfrow=c(3,3))
par(cex=1.5) # Scales the size of points globally
par(cex.axis=0.8, cex.lab=1, cex.main=1) # Adjust font sizes for axis text, labels, and main title
for (i in 1:9) {
  sigma1 <- sort(as.vector(samp$sigma[,,i]))
  sigma2 <- sort(as.vector(Sigmainv_samples[i,]))
  qqplot(sigma1, sigma2)
  abline(0,1)
}
dev.off()
# 
# par(mfrow=c(1,1))
# par(cex=1, cex.axis=1, cex.lab=1, cex.main=1)


i <- 15
sigma1 <- as.vector(samp$sigma[,,i])
sigma2 <- as.vector(Sigmainv_samples[i,])
qqplot(sigma1, sigma2)
abline(0,1)
