source("../R/gen_data.R")
library(matrixsampling)
library(matrixNormal)
library(matrixStats)
library(magrittr)
seed <- 1234
set.seed(seed)
n <- 1
n0 <- 1/2
mu0 <- 2
y <- 1.5
qt <- (n0 * mu0 + n * y) / (n + n0)
Qt <- 1/(n0 + n) + 1/n
st <- qt
St <- 1/(n0 + n)
Vt <- 1/n
## single----
# Analytical
Y_new <- y
dens_ana <- dnorm(x = Y_new, mean = qt, sd = sqrt(Qt))
log_dens_ana <- dnorm(x = Y_new, mean = qt, sd = sqrt(Qt), log = T)

# MC
nsim <- 1000
theta_post <- rnorm(n = nsim, mean = st, sd = sqrt(St))
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