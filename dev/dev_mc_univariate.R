source("../R/gen_data.R")
library(matrixsampling)
library(matrixNormal)
library(matrixStats)
library(magrittr)
seed <- 1234
set.seed(seed)
N <- 1
S <- 1
p <- 1

Ft <- gen_ran_matrix(nrow = N, ncol = p)
st <- gen_ran_matrix(nrow = p, ncol = S)
St <- gen_pd_matrix(dim = p)
Vt <- diag(N)
IS <- diag(S)

qt <- Ft %*% st
Qt <- Vt + Ft %*% St %*% t(Ft)


## single----
# Analytical
Y_new <- rnorm(n = 1, mean = qt, sd = sqrt(Qt))
dens_ana <- dnorm(x = Y_new, mean = qt, sd = sqrt(Qt))
log_dens_ana <- dnorm(x = Y_new, mean = qt, sd = sqrt(Qt), log = T)

# MC
nsim <- 10000
theta_post <- rmatrixnormal(n = nsim, M = st, U = St, V = IS)
for (i in 1:nsim) {
  if(i == 1){
    log_dens_ls <- c()
    dens_sum <- 0
  }
  
  M_temp <- Ft %*% theta_post[,,i]
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