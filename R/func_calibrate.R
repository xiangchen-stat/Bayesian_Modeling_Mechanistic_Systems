scale_uniform <- function(x, low, high){
  if (all(eta >= eta_limit_low) && all(eta <= eta_limit_high)) {
    res <- (x - low) / (high - low)
  } else{
    stop("eta is out of range.")
  }
  res[is.nan(res)] <- 0
  return(res)
}

scale_back_uniform <- function(x, low, high){
  res <- x * (high - low) + low
  return(res)
}

expit <- function(x){
  return((1 + exp(1)^(-x))^(-1))
}

inv_chol <- function(X){
  RX <- chol(X)
  invX <- chol2inv(RX)
  return(invX)
}

cal_jacobian_logit_uniform <- function(eta, eta_limit_low, eta_limit_high, log = FALSE){
  ind_nonzero = eta != 0
  res <- abs(prod(1/(eta[ind_nonzero] - eta_limit_low[ind_nonzero]) + 1/(eta_limit_high[ind_nonzero] - eta[ind_nonzero])))
  if (log == TRUE) {
    res <- log(res)
  }
  return(res)
}

gen_calibrate_data <- function(ycal_mat, para_gen_cal, U_gen){
  res <- list()
  tau2 <- c()
  t_cal <- dim(ycal_mat)[1]
  s_cal <- dim(ycal_mat)[2]
  u <- matrix(nrow = dim(ycal_mat)[1], ncol = dim(ycal_mat)[2])
  z <- matrix(nrow = dim(ycal_mat)[1], ncol = dim(ycal_mat)[2])
  for (t in 1:t_cal) {
    # print(paste(t,"/",t_cal))
    y_t <- ycal_mat[t,]
    alpha <- para_gen_cal$b^t * para_gen_cal$n0
    beta <- para_gen_cal$b^t * para_gen_cal$d0
    tau2_t <- LaplacesDemon::rinvgamma(n = 1, shape = alpha, scale = beta)
    tau2 <- c(tau2, tau2_t)
    e_u <- mniw::rmNorm(n = 1, mu = rep(0, s_cal), Sigma = tau2_t * U_gen)
    if (t == 1) {
      u_t <- para_gen_cal$u0 + e_u
      u[t,] <- u_t
    } else{
      u_t <- u_t + e_u
      u[t,] <- u_t
    }
    e_z <- rnorm(n = s_cal, mean = 0, sd = sqrt(tau2_t))
    z[t,] <- y_t + u_t + e_z
  }
  res$u <- u
  res$z <- z
  res$tau2 <- tau2
  return(res)
}

gen_calibrate_data_uncorr <- function(ycal_mat, para_gen_cal, U_gen){
  res <- list()
  tau2 <- c()
  t_cal <- dim(ycal_mat)[1]
  s_cal <- dim(ycal_mat)[2]
  z <- matrix(nrow = dim(ycal_mat)[1], ncol = dim(ycal_mat)[2])
  for (t in 1:t_cal) {
    y_t <- ycal_mat[t,]
    alpha <- para_gen_cal$b^t * para_gen_cal$n0
    beta <- para_gen_cal$b^t * para_gen_cal$d0
    tau2_t <- LaplacesDemon::rinvgamma(n = 1, shape = alpha, scale = beta)
    tau2 <- c(tau2, tau2_t)
    # e_u <- mniw::rmNorm(n = 1, mu = rep(0, s_cal), Sigma = tau2_t * U_gen)
    e_u <- rnorm(n = s_cal, mean = 0, sd = sqrt(tau2_t * U_gen[1,1]))
    if (t == 1) {
      u_t <- para_gen_cal$u0 + e_u
    } else{
      u_t <- u_t + e_u
    }
    e_z <- rnorm(n = s_cal, mean = 0, sd = sqrt(tau2_t))
    z[t,] <- y_t + u_t + e_z
  }
  res$z <- z
  res$tau2 <- tau2
  return(res)
}


gen_prior_u_tau2 <- function(n0, d0, b, m0, M0, U, nT_cal){
  res <- list()
  s_cal <- length(m0)
  tau2 <- c()
  u <- array(dim = c(nT_cal, s_cal))
  zero_s <- rep(0, s_cal)
  # t = 0
  tau2_t <- LaplacesDemon::rinvgamma(n = 1, shape = n0, scale = d0)
  tau2_0 <- tau2_t
  u_t <- mniw::rmNorm(n = 1, mu = m0, Sigma = tau2_t * M0)
  u0 <- u_t
  # t = 1:T  
  for (t in 1:nT_cal) {
    alpha <- b^t * n0
    beta <- b^t * d0
    tau2_t <- LaplacesDemon::rinvgamma(n = 1, shape = alpha, scale = beta)
    tau2 <- c(tau2, tau2_t)
    e_u <- mniw::rmNorm(n = 1, mu = zero_s, Sigma = tau2_t * U)
    u_t <- u_t + e_u
    u[t,] <- u_t
  }
  
  res$u0 <- u0
  res$u <- u
  res$tau2_0 <- tau2_0
  res$tau2 <- tau2
  return(res)
}

# gen y | Y_1:T
library(deSolve)
library(ReacTran)

SIR <- function(t, y, parms) {
  Nx <- parms$Nx
  Ny <- parms$Ny
  N <- parms$N
  eta1 <- parms$eta[1]
  eta2 <- parms$eta[2]
  alpha1 <- parms$eta[3]
  alpha2 <- parms$eta[4]
  alpha3 <- parms$eta[5]
  Gridx <- parms$Gridx
  Gridy <- parms$Gridy
  
  S <- matrix(nrow = Nx, ncol = Ny, data = y[1:(Nx*Ny)])
  I <- matrix(nrow = Nx, ncol = Ny, data = y[(Nx*Ny+1) : (2*Nx*Ny)])
  R <- matrix(nrow = Nx, ncol = Ny, data = y[(2*Nx*Ny+1) : (3*Nx*Ny)])
  #NOTE: tran.2D as called here is equivalent to computing the Laplacian of S, times alpha1;
  #	with dC being the result.
  dS <- -eta1*S*I/N +
    ReacTran::tran.2D(C = S, D.x = alpha1, D.y = alpha1,
                      dx = Gridx, dy= Gridy,
                      C.x.up = 0,
                      C.y.up = 0,
                      C.x.down = 0,
                      C.y.down = 0)$dC
  dI <- eta1*S*I/N - eta2*I +
    ReacTran::tran.2D(C = I, D.x = alpha2, D.y = alpha2,
                      dx = Gridx, dy = Gridy,
                      C.x.up = 0,
                      C.y.up = 0,
                      C.x.down = 0,
                      C.y.down = 0)$dC
  dR <- eta2*I + 
    ReacTran::tran.2D(C = R, D.x = alpha3, D.y = alpha3,
                      dx = Gridx, dy = Gridy,
                      C.x.up = 0,
                      C.y.up = 0,
                      C.x.down = 0,
                      C.y.down = 0)$dC
  list(c(dS, dI, dR))
}


gen_pde <- function(eta, Nx, Ny, N, nT_ori){
  # NOTE: Nx is length of x-dimension, or rows, Ny is the length of the y-dimension, or columns.
  # set up parameters
  Nx <- Nx
  Ny <- Ny
  N = N
  times <- seq(0, nT_ori - 1)
  
  # generate grid
  locs <- expand.grid(y=1:Ny, x=1:Nx)
  locs <- locs[,c("x", "y")]
  locs = as.numeric(rownames(locs))
  S = length(locs)
  Gridx <- ReacTran::setup.grid.1D(x.up = 0, x.down = Nx, N = Nx)
  Gridy <- ReacTran::setup.grid.1D(x.up = 0, x.down = Ny, N = Ny)
  
  # set up initial values for SIR
  S0 <- matrix(nrow = Nx, ncol = Ny, data = N)
  I0 <- matrix(nrow = Nx, ncol = Ny, data = 0)
  I0[Nx/2,Ny/2] <- 1 # set midpoint to have 50 infected people
  I0[round(Nx/10), round(Ny/10)] <- N * 0.05
  R0 <- matrix(nrow = Nx, ncol = Ny, data = 0)
  yini <- c(S0, I0, R0)
  
  # pass parms to SIR within ode.2D
  parms <- list(eta = eta, Gridx = Gridx, Gridy = Gridy, Nx = Nx, Ny = Ny, N = N,
                y = yini, t = times)
  
  # Generate SIR data from the SIR ODE.
  out <- deSolve::ode.2D(y = yini, times = times, func = SIR, parms = parms,
                         nspec = 3, dimens = c(Nx, Ny), 
                         lrw = 200000000, names=c("X1", "X2", "X3"))[,-1]
  lgout <- log(out + 1)
  res <- lgout[,(1 + Nx * Ny):(2 * Nx * Ny)]
  # if(length(subset) != 1 || !is.na(subset)){
  #   res <- res[,subset]
  # }
  return(res)
}


gen_Jt <- function(input, input_new, phi_para, gp_sigma2 = 1.1, gp_tau2 = 10^(-4)){
  # deal with one input prediction
  if(is.vector(input_new)){
    input_new <- t(input_new)
  }
  N <- dim(input)[1]
  
  # get spatial kernel
  input_full <- rbind(input, input_new)
  dist_input_full <- as.matrix(stats::dist(input_full, method = "euclidean", diag = T, upper = T))
  V_full <- gen_exp_kernel(loc = input_full, phi = phi_para, sigma2 = gp_sigma2, tau2 = gp_tau2) # exponential kernel
  Jt <- V_full[1:N, (N+1):dim(V_full)[2]]
  return(Jt)
}

# TODO sigma2_cal using mean instead of all Monte Carlo samples
update_y_eta <- function(eta, nsam, Nx, Ny, N_people, nT_ori, nT, nT_cal, N_sp, N_sp_train, AR_choice, bncol, 
                     res_ffbs, n_b_cal, ind_block_cal, t_VinvY_FTheta_ls, Vinv, R_cal_train, input, input_cal, phi_para, 
                     gp_sigma2, gp_tau2,
                     ind_train_cal){
  res <- list()
  # 1. predict y_eta using PDE
  yt_eta_tt <- gen_pde(eta = eta, Nx = Nx, Ny = Ny, N = N_people, 
                       nT_ori = nT_ori)
  yt_eta_tt_ls <- list()
  for (i in 1:dim(yt_eta_tt)[1]) {
    yt_eta_tt_ls[[i]] <- t(yt_eta_tt[i,])
  }
  
  # 2. calculate post mean and variance
  # 2.1 generate F_eta
  if(AR_choice == 1){
    F_eta_ls <- gen_F_ls_AR1_EP(Y = yt_eta_tt_ls, nT = nT_ori, n_b = n_b_cal, ind = ind_block_cal)
  } else if(AR_choice == 2){
    F_eta_ls <- gen_F_ls_AR2_EP(Y = yt_eta_tt_ls, nT = nT_ori, n_b = n_b_cal, ind = ind_block_cal)
  }
  
  F_eta_Theta_EP <- list()
  for (t in 1:nT) {
    F_eta_Theta_EP[[t]] <- array(dim = c(1, bncol, nsam))
    F_eta_ls_t <- F_eta_ls[[t]]
    for (i in 1:nsam) {
      F_eta_Theta_EP[[t]][,,i] <- F_eta_ls_t %*% res_ffbs[[t]][,,i]
    }
  }
  F_eta_Theta <- recover_from_EP_MC(dat_EP = F_eta_Theta_EP, nT_ori = nT_ori, nT = nT, nsam = nsam)
  
  # 2.2 generate J_eta
  J_eta <- gen_Jt(input = input, input_new = input_cal, phi_para = phi_para, 
                  gp_sigma2 = gp_sigma2, gp_tau2 = gp_tau2)
  # 2.3 calculate mu
  mu_arry <- array(dim = c(nT_cal, N_sp, nsam))
  for (t in 1:nT_cal) {
    t_VinvY_FTheta_t <- t_VinvY_FTheta_ls[[t]]
    for (i in 1:nsam) {
      mu_arry[t,,i] <- F_eta_Theta[[t]][,,i] + as.vector(t_VinvY_FTheta_t[,,i] %*% J_eta)
    }
  }
  
  mu_arry_train <- array(dim = c(nT_cal, N_sp_train, nsam))
  for (t in 1:nT_cal) {
    mu_arry_train[t,,] <- mu_arry[t,ind_train_cal,]
  }
  # 2.4 calculate Sigma constant
  Sigma_eta_c <- 1 - crossprod(J_eta, Vinv) %*% J_eta
  # TODO 
  sigma2_cal <- as.numeric(mean(res_ffbs$Sigma) * Sigma_eta_c)
  Sigma_eta_train <- sigma2_cal * R_cal_train
  # 3. sampling
  y_arry <- array(dim = c(nT_cal, N_sp_train, nsam))
  for (t in 1:nT_cal) {
    for (i in 1:nsam) {
      y_arry[t,,i] <- mniw::rmNorm(n = 1, mu = mu_arry_train[t,,i], 
                                    Sigma = Sigma_eta_train)
    }
  }
  
  res$y <- y_arry
  res$mu <- mu_arry_train
  res$sigma2_cal <- sigma2_cal
  res$Sigma <- Sigma_eta_train
  return(res)
}

# TODO sigma2_cal using mean instead of all Monte Carlo samples
update_y_eta_one <- function(eta, ind_sam, Nx, Ny, N_people, nT_ori, nT, nT_cal, N_sp, N_sp_train, AR_choice, bncol, 
                      res_ffbs, n_b_cal, ind_block_cal, t_VinvY_FTheta_ls, Vinv, R_cal_train, input, input_cal, phi_para, 
                      gp_sigma2, gp_tau2, ind_train_cal){
  res <- list()
  nsam <- 1 # generate only one sample
  # 1. predict y_eta using PDE
  yt_eta_tt <- gen_pde(eta = eta, Nx = Nx, Ny = Ny, N = N_people, 
                       nT_ori = nT_ori)
  yt_eta_tt_ls <- list()
  for (i in 1:dim(yt_eta_tt)[1]) {
    yt_eta_tt_ls[[i]] <- t(yt_eta_tt[i,])
  }
  
  # 2. calculate post mean and variance
  # 2.1 generate F_eta
  if(AR_choice == 1){
    F_eta_ls <- gen_F_ls_AR1_EP(Y = yt_eta_tt_ls, nT = nT_ori, n_b = n_b_cal, ind = ind_block_cal)
  } else if(AR_choice == 2){
    F_eta_ls <- gen_F_ls_AR2_EP(Y = yt_eta_tt_ls, nT = nT_ori, n_b = n_b_cal, ind = ind_block_cal)
  }
  
  F_eta_Theta_EP <- list()
  for (t in 1:nT) {
    F_eta_ls_t <- F_eta_ls[[t]]
    for (i in ind_sam) {
      F_eta_Theta_EP[[t]] <- F_eta_ls_t %*% res_ffbs[[t]][,,i]
    }
  }
  F_eta_Theta <- recover_from_EP_exact(dat_EP = F_eta_Theta_EP, nT_ori = nT_ori, nT = nT)
  
  # 2.2 generate J_eta
  J_eta <- gen_Jt(input = input, input_new = input_cal, phi_para = phi_para, 
                  gp_sigma2 = gp_sigma2, gp_tau2 = gp_tau2)
  # 2.3 calculate mu
  mu_arry <- array(dim = c(nT_cal, N_sp))
  for (t in 1:nT_cal) {
    t_VinvY_FTheta_t <- t_VinvY_FTheta_ls[[t]]
    for (i in ind_sam) {
      mu_arry[t,] <- F_eta_Theta[[t]] + as.vector(t_VinvY_FTheta_t[,,i] %*% J_eta)
    }
  }
  
  mu_arry_train <- array(dim = c(nT_cal, N_sp_train))
  for (t in 1:nT_cal) {
    mu_arry_train[t,] <- mu_arry[t,ind_train_cal]
  }
  # 2.4 calculate Sigma constant
  Sigma_eta_c <- 1 - crossprod(J_eta, Vinv) %*% J_eta
  # TODO 
  sigma2_cal <- as.numeric(res_ffbs$Sigma[ind_sam] * Sigma_eta_c)
  Sigma_eta_train <- sigma2_cal * R_cal_train
  # 3. sampling
  y_arry <- array(dim = c(nT_cal, N_sp_train))
  for (t in 1:nT_cal) {
    for (i in ind_sam) {
      y_arry[t,] <- mniw::rmNorm(n = 1, mu = mu_arry_train[t,],
                                  Sigma = Sigma_eta_train)
    }
  }

  res$y <- y_arry
  res$mu <- mu_arry_train
  res$Sigma <- Sigma_eta_train
  return(res)
}

update_muSigma_eta_one <- function(eta, ind_sam, Nx, Ny, N_people, nT_ori, nT, nT_cal, N_sp, N_sp_train, AR_choice, bncol, 
                         res_ffbs, n_b_cal, ind_block_cal, t_VinvY_FTheta_ls, Vinv, R_cal_train, input, input_cal, phi_para, 
                         gp_sigma2, gp_tau2, ind_train_cal){
  res <- list()
  nsam <- 1 # generate only one sample
  # 1. predict y_eta using PDE
  yt_eta_tt <- gen_pde(eta = eta, Nx = Nx, Ny = Ny, N = N_people, 
                       nT = nT_ori)
  yt_eta_tt_ls <- list()
  for (i in 1:dim(yt_eta_tt)[1]) {
    yt_eta_tt_ls[[i]] <- t(yt_eta_tt[i,])
  }
  
  # 2. calculate post mean and variance
  # 2.1 generate F_eta
  if(AR_choice == 1){
    F_eta_ls <- gen_F_ls_AR1_EP(Y = yt_eta_tt_ls, nT = nT_ori, n_b = n_b_cal, ind = ind_block_cal)
  } else if(AR_choice == 2){
    F_eta_ls <- gen_F_ls_AR2_EP(Y = yt_eta_tt_ls, nT = nT_ori, n_b = n_b_cal, ind = ind_block_cal)
  }
  
  F_eta_Theta_EP <- list()
  for (t in 1:nT) {
    F_eta_ls_t <- F_eta_ls[[t]]
    for (i in ind_sam) {
      F_eta_Theta_EP[[t]] <- F_eta_ls_t %*% res_ffbs[[t]][,,i]
    }
  }
  F_eta_Theta <- recover_from_EP_exact(dat_EP = F_eta_Theta_EP, nT_ori = nT_ori, nT = nT)
  
  # 2.2 generate J_eta
  J_eta <- gen_Jt(input = input, input_new = input_cal, phi_para = phi_para, 
                  gp_sigma2 = gp_sigma2, gp_tau2 = gp_tau2)
  # 2.3 calculate mu
  mu_arry <- array(dim = c(nT_cal, N_sp))
  for (t in 1:nT_cal) {
    t_VinvY_FTheta_t <- t_VinvY_FTheta_ls[[t]]
    for (i in ind_sam) {
      mu_arry[t,] <- F_eta_Theta[[t]] + as.vector(t_VinvY_FTheta_t[,,i] %*% J_eta)
    }
  }
  
  mu_arry_train <- array(dim = c(nT_cal, N_sp_train))
  for (t in 1:nT_cal) {
    mu_arry_train[t,] <- mu_arry[t,ind_train_cal]
  }
  # 2.4 calculate Sigma constant
  Sigma_eta_c <- 1 - crossprod(J_eta, Vinv) %*% J_eta
  # TODO 
  sigma2_cal <- as.numeric(mean(res_ffbs$Sigma) * Sigma_eta_c)
  Sigma_eta_train <- sigma2_cal * R_cal_train
  # # 3. sampling
  # y_arry <- array(dim = c(nT_cal, N_sp_train))
  # for (t in 1:nT_cal) {
  #   for (i in ind_sam) {
  #     y_arry[t,] <- mniw::rmNorm(n = 1, mu = mu_arry_train[t,], 
  #                                 Sigma = Sigma_eta_train)
  #   }
  # }
  # 
  # res$y <- y_arry
  res$mu <- mu_arry_train
  res$Sigma <- Sigma_eta_train
  return(res)
}

sample_y_eta_one <- function(ind_sam, nT_cal, N_sp_train, mu_arry_train, Sigma_eta_train){
  res <- list()
  # 3. sampling
  y_arry <- array(dim = c(nT_cal, N_sp_train))
  for (t in 1:nT_cal) {
    for (i in ind_sam) {
      y_arry[t,] <- mniw::rmNorm(n = 1, mu = mu_arry_train[t,], 
                                 Sigma = Sigma_eta_train)
    }
  }
  
  res = y_arry
  return(res)
}

cal_Bt_bt <- function(mut, invSigma, tau2t, zt_ut, I_Stilde){
  res <- list()
  invBt <- invSigma + 1/tau2t * I_Stilde
  Bt <- inv_chol(invBt)
  bt <- invSigma %*% mut + 1/tau2t * zt_ut 
  Btbt <- Bt %*% bt
  res$Bt <- Bt
  res$Btbt <- Btbt
  return(res)
}