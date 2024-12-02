require(matrixsampling)
#' Generate random matrix
#' @param nrow Number of rows
#' @param col Number of columns
#' @returns A random matrix
#' @export
gen_ran_matrix <- function(nrow, ncol){
  S <- matrix(rnorm(nrow * ncol), nrow  = nrow, ncol  = ncol)
  S <- S / abs(max(S))
  return(S)
}

#' Generate random positive definite matrix
#' @param dim Number of rows/columns
#' @returns A random pd matrix
#' @export
gen_pd_matrix <- function(dim){
  tempS <- matrix(rnorm(dim * dim), nrow  = dim, ncol  = dim)
  S <- crossprod(tempS)
  S <- S / abs(max(S))
  return(S)
}

#' Generate multivariate normal
#' @param m Mean matrix
#' @param RM Upper-triangular Cholesky factor of row covariance matrix
#' @param RSigma Upper-triangular Cholesky factor of column covariance matrix
#' @returns A random multivariate normal matrix
#' @export
rmn_chol <- function(m, RM, RSigma){
  p <- nrow(m)
  S <- ncol(m)
  vec <- rnorm(p * S)
  Z <- matrix(vec, nrow = p, ncol = S)
  theta <- m + crossprod(RM, Z) %*% RSigma
  return(theta)
}


#' Generate multivariate normal
#' @param nsam Number of samples
#' @param m Mean matrix
#' @param RM Upper-triangular Cholesky factor of row covariance matrix
#' @param RSigma Upper-triangular Cholesky factor of column covariance matrix
#' @returns A random multivariate normal matrix
#' @export
rmn_chol_more <- function(nsam, m, RM, RSigma){
  p <- nrow(m)
  S <- ncol(m)
  out <- array(dim = c(p, S, nsam))
  for (i in 1:nsam) {
    vec <- rnorm(p * S)
    Z <- matrix(vec, nrow = p, ncol = S)
    theta <- m + crossprod(RM, Z) %*% RSigma
    out[,,i] <- theta
  }
  
  return(out)
}


#' Generate FFBS data in workspace
#' @param N Number of rows of Y
#' @param S Number of columns of Y
#' @param p Number of rows of Theta
#' @param nT Number of time stamps
#' @returns A list of .csv files: Y's m0 M0 G0 F0 Sigma RSigma(chol(Sigma)) V0 RVO(chol(V0)) W0
#' @export
gen_ffbs_data <- function(N, S, p, nT){
  out <- list()
  para <- list()
  Y_out <- list()
  m0 <- gen_ran_matrix(p, S)
  M0 <- gen_pd_matrix(p)
  n0 <- S + 2
  D0 <- gen_pd_matrix(S)
  Sigmapre <- rinvwishart(1, nu = n0, Omega = D0, epsilon = 0, checkSymmetry = TRUE)
  Sigma <- matrix(Sigmapre, ncol = S, nrow = S)
  G0 <- diag(1, nrow = p)
  # generate F0
  F0 <- matrix(nrow = N, ncol = p)
  if(p > 1){
    F0[,1] <- 1
    for (i in 2:p) {
      F0[,i] <- rnorm(n = N)
    }
  } else{
    stop("p is less than 2")
  }
  # generate V0 as Gaussian kernel
  loc <- F0 # location
  dist <- as.matrix(stats::dist(loc, method = "euclidean", diag = T, upper = T))
  phi <- 0.4 * 3 / max(dist)
  sigma2 <- 1
  tau2 <- 0.01
  # V0 <- exp(-phi * dist)
  V0 <- gen_gp_kernel(loc = F0, phi = phi, sigma2 = sigma2, tau2 = tau2)
  # V0 <- gen_pd_matrix(N)
  W0 <- gen_pd_matrix(p)
  RM0 <- chol(M0)
  RSigma <- chol(Sigma)
  RV0 <- chol(V0)
  zero_NS <- matrix(0, nrow = N, ncol = S)
  zero_pS <- matrix(0, nrow = p, ncol = S)
  
  # generate and save Y's
  for(i in 1:nT){
    if(i == 1){
      Theta = rmn_chol(m0, RM0, RSigma)
      E <- rmn_chol(zero_NS, RV0, RSigma)
      Y <- F0 %*% Theta + E
    }
    else{
      Gamma <- rmn_chol(zero_pS, W0, RSigma)
      Theta <- G0 %*% Theta + Gamma
      E <- rmn_chol(zero_NS, RV0, RSigma)
      Y <- F0 %*% Theta + E
    }
    Yl <- list(Y)
    names(Yl) <- paste("Y", i, sep = "")
    Y_out <- append(Y_out, Yl)
  }
  
  
  # save other matrices
  para$loc <- loc
  para$dist <- dist
  para$m0 <- m0
  para$M0 <- M0
  para$G0 <- G0
  para$F0 <- F0
  para$Sigma <- Sigma
  para$RSigma <- RSigma
  para$V0 <- V0
  para$RV0 <- RV0
  para$W0 <- W0
  para$n0 <- n0
  para$D0 <- D0
  
  out <- list(Y = Y_out, para = para)
  return(out)
}


#' Generate FFBS related data
#' @param N Number of rows of Y
#' @param S Number of columns of Y
#' @param p Number of rows of Theta
#' @param nT Number of time stamps
#' @param path Path for saving generated data
#' @returns A list of .csv files: Y's m0 M0 G0 F0 Sigma RSigma(chol(Sigma)) V0 RVO(chol(V0)) W0
#' @export
gen_ffbs_csv <- function(N, S, p, nT, path){
  m0 <- gen_ran_matrix(p, S)
  M0 <- gen_pd_matrix(p)
  Sigma <- gen_pd_matrix(S)
  G0 <- diag(1, nrow = p)
  # generate F0
  F0 <- matrix(nrow = N, ncol = p)
  if(p > 1){
    F0[,1] <- 1
    for (i in 2:p) {
      F0[,i] <- rnorm(n = N)
    }
  } else{
    stop("p is less than 2")
  }
  # generate V0 as Gaussian kernel
  loc <- F0 # location
  dist <- as.matrix(stats::dist(loc, method = "euclidean", diag = T, upper = T))
  phi <- 0.4 * 3 / max(dist)
  sigma2 <- 1
  tau2 <- 0.01
  # V0 <- exp(-phi * dist)
  V0 <- gen_gp_kernel(loc = F0, phi = phi, sigma2 = sigma2, tau2 = tau2)
  # V0 <- gen_pd_matrix(N)
  W0 <- gen_pd_matrix(p)
  RM0 <- chol(M0)
  RSigma <- chol(Sigma)
  RV0 <- chol(V0)
  zero_NS <- matrix(0, nrow = N, ncol = S)
  zero_pS <- matrix(0, nrow = p, ncol = S)
  
  # save big matrices
  write.csv(x = m0, file = paste0(path, "m0.csv"), row.names=FALSE)
  write.csv(x = M0, file = paste0(path, "MM0.csv"), row.names=FALSE)
  write.csv(x = G0, file = paste0(path, "G0.csv"), row.names=FALSE)
  write.csv(x = F0, file = paste0(path, "F0.csv"), row.names=FALSE)
  write.csv(x = Sigma, file = paste0(path, "Sigma.csv"), row.names=FALSE)
  write.csv(x = RSigma, file = paste0(path, "RSigma.csv"), row.names=FALSE)
  write.csv(x = V0, file = paste0(path, "V0.csv"), row.names=FALSE)
  write.csv(x = RV0, file = paste0(path, "RV0.csv"), row.names=FALSE)
  write.csv(x = W0, file = paste0(path, "W0.csv"), row.names=FALSE)
  
  # generate and save Y's
  for(i in 1:nT){
    if(i == 1){
      Theta = rmn_chol(m0, RM0, RSigma)
      E <- rmn_chol(zero_NS, RV0, RSigma)
      Y <- F0 %*% Theta + E
    }
    else{
      Gamma <- rmn_chol(zero_pS, W0, RSigma)
      Theta <- G0 %*% Theta + Gamma
      E <- rmn_chol(zero_NS, RV0, RSigma)
      Y <- F0 %*% Theta + E
    }
    write.csv(x = Y, file = paste0(path, "Y", i, ".csv"), row.names=FALSE)
  }
  
}


## calculate Gaussian Process kernel
gen_gp_kernel <- function(loc, phi, sigma2, tau2){
  dist_k <- as.matrix(stats::dist(loc, method = "euclidean", diag = T, upper = T))
  n <- dim(dist_k)[1]
  var_gp <- sigma2 * exp(-phi * dist_k) + tau2 * (dist_k == 0)
  return(var_gp)
}

## calculate Gaussian Process kernel
gen_expsq_kernel <- function(loc, phi, sigma2 = 1, tau2 = 0){
  dist_k <- as.matrix(stats::dist(loc, method = "euclidean", diag = T, upper = T))
  n <- dim(dist_k)[1]
  var_gp <- sigma2 * exp(-phi * dist_k^2) + tau2 * (dist_k == 0)
  return(var_gp)
}

gen_exp_kernel <- function(loc, phi, sigma2 = 1, tau2 = 0){
  dist_k <- as.matrix(stats::dist(loc, method = "euclidean", diag = T, upper = T))
  n <- dim(dist_k)[1]
  var_gp <- sigma2 * exp(-phi * dist_k) + tau2 * (dist_k == 0)
  return(var_gp)
}


# create F_ls using auto-regressive
gen_F_ls_AR1 <- function(Y, nT){
  F_ls <- list()
  for (i in 1:nT) {
    if(i == 1){
      Ft <- Y[[1]]
      F_ls[[i]] <- Ft
    } else{
      Ft <- Y[[i - 1]]
      F_ls[[i]] <- Ft
    }
  }
  return(F_ls)
}


gen_F_ls_AR2 <- function(Y, nT){
  F_ls <- list()
  for (i in 1:nT) {
    if(i == 1){
      # Ft <- cbind(1/2 * (Y[[2]] + Y[[3]]), 
      #             1/2 * (Y[[2]] + Y[[4]]))
      # Ft <- cbind(gen_ran_matrix(nrow = nrow(Y[[1]]), ncol = ncol(Y[[1]])),
      #             gen_ran_matrix(nrow = nrow(Y[[1]]), ncol = ncol(Y[[1]])))
      Ft <- cbind(Y[[1]],
                  0.01 * gen_ran_matrix(nrow = nrow(Y[[1]]), ncol = ncol(Y[[1]])))
      F_ls[[i]] <- Ft
    } else if(i == 2){
      Ft <- cbind(Y[[1]],
                  1/2 * (Y[[1]] + Y[[3]]))
      # Ft <- cbind(1/2 * (Y[[1]] + Y[[3]]), 
      #             gen_ran_matrix(nrow = nrow(Y[[1]]), ncol = ncol(Y[[1]])))
      F_ls[[i]] <- Ft
    } else{
      Ft <- cbind(Y[[i - 2]], Y[[i - 1]])
      # Ft <- cbind(1/2 * (Y[[i]] + Y[[i - 1]]),
      #             Y[[i - 1]])
      F_ls[[i]] <- Ft
    }
  }
  return(F_ls)
}

gen_F_ls_AR1_EP <- function(Y, nT, n_b, ind){
  out <- list()
  F_ls_full <- gen_F_ls_AR1(Y = Y, nT = nT)
  for (i in 1:nT) {
    temp_F <- F_ls_full[[i]]
    for (j in 1:n_b) {
      out[[(i-1)*n_b + j]] <- temp_F[ind[j,2]:ind[j,3], ind[j,4]:ind[j,5]]
    }
  }
  # print("Ft using AR1")
  return(out)
}

gen_F_ls_AR2_EP <- function(Y, nT, n_b, ind){
  out <- list()
  S <- ncol(Y[[1]])
  F_ls_full <- gen_F_ls_AR2(Y = Y, nT = nT)
  for (i in 1:nT) {
    temp_F <- F_ls_full[[i]]
    for (j in 1:n_b) {
      out[[(i-1)*n_b + j]] <- temp_F[ind[j,2]:ind[j,3], 
                                     c(ind[j,4]:ind[j,5], (S + (ind[j,4]:ind[j,5])))]
    }
  }
  # print("Ft using AR2")
  return(out)
}


# recover form episode season model to regular data dimension
recover_from_EP_exact <- function(dat_EP, nT_ori, nT){
  if(length(dat_EP) != nT){
    stop("dimension of dat_EP doesn't match nT")
  }
  
  n_block <- nT / nT_ori
  out <- list()
  
  # cbind n_block of EP blocks to recover season
  for (i in 1:nT_ori) {
    out[[i]] <- as.matrix(as.data.frame(dat_EP[(1+(i-1)*n_block):(i*n_block)]))
  }
  
  names(out) <- paste0("T", seq(1:nT_ori))
  return(out)
}

# recover for list of arrays
recover_from_EP_MC <- function(dat_EP, nT_ori, nT, nsam){
  if(length(dat_EP) != nT){
    stop("dimension of dat_EP doesn't match nT")
  }
  
  n_block <- nT / nT_ori
  dim_block <- dim(dat_EP[[1]])
  dim_no_block <- dim_block
  dim_no_block[2] <- dim_block[2] * n_block
  n_input <- dim(dat_EP)
  out <- list()
  
  # transform list to array to simplify modification
  dat_ary <- array(dim = c(dim_block, nT))
  for (i in 1:nT) {
    dat_ary[,,,i] <- dat_EP[[i]]
  }
  
  # cbind n_block of EP blocks to recover season
  for (i in 1:nT_ori) {
    out[[i]] <- array(dim = dim_no_block)
    for (j in 1:nsam) {
      out[[i]][,,j] <- as.matrix(as.data.frame(dat_ary[,,j,(1+(i-1)*n_block):(i*n_block)]))
    }
  }
  
  names(out) <- paste0("T", seq(1:nT_ori))
  return(out)
}
