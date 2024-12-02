library(MNIW)
library(matrixsampling)
library(microbenchmark)
library(Rcpp)
source("../R/gen_data.R")
source("../R/BS.R")
source("../R/pds.R")
source("../R/MNIW_sampler.R")
source("naiive_BS.R")
source("test_utils.R")
source("naiive_FF.R")
source("naiive_FFBS.R")
sourceCpp("../src/FF.cpp")
sourceCpp("../src/BS.cpp")

seed <- 1234
set.seed(seed)
path_data <- file.path("..", "data")
path_fig <- file.path("..", "figures")

# read in pde data----
n_input <- 50
dt_pde <- list()
for (i in 1:n_input) {
  dt_pde[[i]] <- as.matrix(read.csv(file = paste0(path_data, "/pde_solution_", i, ".csv", sep = "")))
}
pde_para <- as.matrix(read.csv(file = paste0(path_data, "/pde_para.csv", sep = "")))[1:n_input,]
nT <- dim(dt_pde[[1]])[1]
Nx <- 10
Ny <- 10
ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))

# Get covariance matrix
dist_para <- as.matrix(stats::dist(pde_para, method = "euclidean", diag = T, upper = T))
phi_para <- 0.4 * 3 / max(dist_para)
V_para <- gen_expsq_kernel(loc = pde_para, phi = phi_para, sigma2 = 1, tau2 = 0.01) # exponential kernel

dist_sp <- as.matrix(stats::dist(ind_sp, method = "euclidean", diag = T, upper = T))
phi_sp <- 0.4 * 3 / max(dist_sp)
sigma_sp <- gen_gp_kernel(loc = ind_sp, phi = phi_sp, sigma2 = 1.1, tau2 = 0.01)

# transform the pde data
for (i in 1:nT) {
  temp <- c()
  for (j in 1:n_input) {
    temp <- rbind(temp, dt_pde[[j]][i,])
  }
  write.csv(temp,
            file=file.path(path_data, paste0("pde_reform_", i, ".csv", sep = "")),
            row.names=FALSE)
}

# read in transformed data
dt_pde_reform <- list()
for (i in 1:nT) {
  dt_pde_reform[[i]] <- as.matrix(read.csv(file = paste0(path_data, "/pde_reform_", i, ".csv", sep = "")))
}
names(dt_pde_reform) <- paste0("T", seq(1:nT))

# Initialize para
Y = dt_pde_reform
D0 = sigma_sp
F_ls <- cbind(rep(1, n_input), rowSums(pde_para))
V_ls <- V_para
N <- n_input
S <- Nx * Ny
p <- dim(F_ls)[2]
nsam <- 10
G_ls <- diag(p)
W_ls <- diag(p)
n0 <- p + 3
m0 <- matrix(1, nrow = p, ncol = S)
M0 <- diag(p)
delta <- 1
# Sys.time()
# dt <- gen_ffbs_data(N = N, S = S, p = p, nT = nT)


## FFBS----
Sys.time()
res_ffbs <- FFBS_R_naiive(nsam = nsam, Y = Y, F_ls = F_ls, G_ls = G_ls, 
                      W_ls = W_ls, V_ls = V_ls, 
                      m0 = m0, M0 = M0, 
                      n0 = n0, D0 = D0, 
                      nT = nT, delta = 1.0)
Sys.time()

## Prediction----
input <- pde_para
input_new <- pde_para
n_input_new <- dim(input_new)[1]
F_new_ls <- cbind(rep(1, n_input_new), rowSums(input_new))


# check prediction on new inputs----
FFBS_predict_MC_R_naiive_old <- function(nsam, Y, res_ffbs, F_ls, F_new_ls, input, input_new,
                                         nT, delta = 1.0){
  out <- list()
  post_mean_ls <- list()
  post_V_ls <- list()
  # get spatial kernel
  input_full <- rbind(input, input_new)
  dist_input_full <- as.matrix(stats::dist(input_full, method = "euclidean", diag = T, upper = T))
  phi_input_full <- 0.4 * 3 / max(dist_input_full)
  V_full <- gen_expsq_kernel(loc = input_full, phi = phi_input_full, sigma2 = 1, tau2 = 0.1) # exponential kernel
  Vt <- V_full[1:N,1:N]
  Jt <- V_full[1:N,(N+1):dim(V_full)[2]]
  Vt_new <- V_full[(N+1):dim(V_full)[1],(N+1):dim(V_full)[2]]
  # prepare matrix
  # Vtinv2 <- solve(Vt)
  Vt_chol <- chol(Vt)
  Vtinv <- chol2inv(Vt_chol)
  tJtVtinv <- crossprod(Jt, Vtinv)
  
  # initialize Ft, Ft_new
  Ft <- F_ls
  Ft_new <- F_new_ls
  RSigma <- array(dim = dim(res_ffbs$Sigma))
  
  for (i in 1:nT) {
    # get Ft's
    if(is.list(F_ls) && length(F_ls) != 1){
      Ft <- F_ls[[i]]
    } 
    if(is.list(F_new_ls) && length(F_new_ls) != 1){
      Ft_new <- F_new_ls[[i]]
    } 
    Yt <- Y[[i]]
    post_mean <- array(dim = c(dim(Yt), nsam))
    post_V <- array(dim = c(dim(Vt_new), nsam))
    
    for(j in 1:nsam){
      # Get Thetat
      Thetat <- res_ffbs[[i]][,,j]
      if(i == 1){
        RSigma[,,j] <- chol(res_ffbs$Sigma[,,j])    # cholesky of Sigma
      }
      RSigma_j <- RSigma[,,j]
      
      # post_mean <- Ft_new %*% Thetat + t(Jt) %*% solve(Vt) %*% (Yt - Ft %*% Thetat)
      post_mean[,,j] <- Ft_new %*% Thetat + tJtVtinv %*% (Yt - Ft %*% Thetat)
      # post_var <- Vt_new - t(Jt) %*% solve(Vt) %*% Jt
      post_V[,,j] <- Vt_new - tJtVtinv %*% Jt
      
    }
    post_mean_ls[[i]] <- post_mean
    post_V_ls[[i]] <- post_V
  }
  names(post_mean_ls) <- paste("T", seq(1:nT), sep = "")
  names(post_V_ls) <- paste("T", seq(1:nT), sep = "")
  out$post_mean <- post_mean_ls
  out$post_V <- post_V_ls
  return(out)
}

# check deterministic interpolation----
res_pre <- FFBS_predict_MC_R_naiive(nsam = nsam, Y = Y, res_ffbs = res_ffbs, 
                                    input = input, input_new = input_new,
                                    F_ls = F_ls, F_new_ls = F_new_ls, 
                                    nT = nT, delta = 1.0)

tol <- 10^(-6)
# check 0 covariance matrix using F norm
for (i in 1:nT) {
  if(i == 1){
    count_V <- 0
    count_mean <- 0
    # check variance matrix
    norm_V <- norm(res_pre$post_V, "F")
    if(norm_V > tol){
      count_V <- count_V + 1
      print(paste("Variance at Time", i, "Sample", j, "has norm", norm_V))
    }
  }
  for (j in 1:nsam) {
    # check mean matrix
    norm_mean <- norm(res_pre$post_mean[[i]][,,j] - Y[[i]], "F")
    if(norm_mean > tol){
      count_mean <- count_mean + 1
      print(paste("Mean at Time", i, "Sample", j, "has norm", norm_mean))
    }
  }
  if(i == nT && j == nsam){
    if(count_V == 0 && count_mean == 0){
      print("All tests passed. It exactly interpolates.")
    } else{
      print(paste(count_V, "norms of variance are >=", tol))
      print(paste(count_mean, "norms of mean are >=", tol))
    }

  }
}

