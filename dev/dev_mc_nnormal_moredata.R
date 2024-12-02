source("../R/gen_data.R")
library(matrixsampling)
library(matrixNormal)
library(matrixStats)
library(magrittr)
seed <- 1234
set.seed(seed)
n <- 1000
sigma2 <- 5
beta <- 10
y <- rnorm(n = n, beta, sigma2)
ybar <- mean(y)
n0 <- 2
mu <- 1
m_post <- (n0 * mu + n * ybar) / (n0 + n)
var_post <- sigma2 / (n0 + n)
L <- 1000
beta_post <- rnorm(n = L, m_post, sd = sqrt(var_post))
qt <- m_post
Qt <- sigma2 *(1 + n / (n0 + n))

## single----
# Analytical
Y_new <- y[1]
dens_ana <- dnorm(x = Y_new, mean = qt, sd = sqrt(Qt))
log_dens_ana <- dnorm(x = Y_new, mean = qt, sd = sqrt(Qt), log = T)

# MC
nsim <- L
theta_post <- beta_post
for (i in 1:nsim) {
  if(i == 1){
    log_dens_ls <- c()
    dens_sum <- 0
  }
  
  M_temp <- theta_post[i]
  dens_temp <- dnorm(x = Y_new, mean = M_temp, sd = sqrt(Vt))
  dens_sum <- dens_sum + dens_temp
  log_dens_temp <- dnorm(x = Y_new, mean = M_temp, sd = sqrt(Vt), log = T)
  log_dens_ls <- c(log_dens_ls, log_dens_temp)
  
  
  if(i == nsim){
    dens_mean_mc = dens_sum / nsim
    log_dens_mc <- logSumExp(log_dens_ls) - log(length(log_dens_ls))
  }
}

dens_ana
dens_mean_mc
log_dens_ana
log_dens_mc

# print(paste("Analytical Density:", dens_ana))
# print(paste("MC Density:", dens_mean_mc))
# print(paste("Diff:", dens_mean_mc - dens_ana))
# print(paste("log-Analytical Density:", log_dens_ana))
# print(paste("log-MC Density:", log_dens_mc))
# print(paste("Diff:", log_dens_mc - log_dens_ana))