# set up emulation ----
PC <- "mac"
switch (PC,
        "x1c" = setwd("D:/Documents/UCLA/0-Administrative/GSR/Bayesian-Modeling-Mechanistic-Systems/dev"),
        "mac" = setwd("/Users/xiangchen/Documents/Bayesian-Modeling-Mechanistic-Systems/dev"),
        "xps" = setwd("C:/Users/wilson/Desktop/X1C/MNIW/dev")
)
if(T){
  source("emulate_EP_sigma2R.R")
}

{
  source("../R/func_calibrate.R")
  source("../R/func_plot.R")
  # inherit value from emulation
  loc_cal <- ind_sp # use the same location as emulation
  nT_cal <- nT_ori
  N_sp_cal <- N_sp # choose no of locations
  no_input_new <- n_test # choose which new input using in calibration
  input_cal <- input_new[no_input_new,]
  eta_true <- input_cal
  n_eta <- 5
  # eta_low = c(2, 0.2, rep(0, 3)) # e1 [2, 4] e2: [0.2, .4] # a1-a3 [0, 0.2]
  # eta_high = c(4, 0.4, rep(0.2, 3))
  eta_low = c(6, 2.2, 0, 1.4, 0) # e1 [6, 7] e2: [2.2, 2.4] a2 [1.4, 1.6] # a1, a3 [0, 0]
  eta_high = c(7, 2.4, 0, 1.6, 0)
  ycal_mat <- gen_pde(eta = input_cal, Nx = Nx, Ny = Ny, N = N_people, 
                      nT_ori = nT_ori)
  
  # # Compare with original generated data
  # y_emu <- c()
  # for (i in 1:length(dt_pde_test)) {
  #   y_emu <- rbind(y_emu, t(dt_pde_test[[i]][no_input_new,]))
  # }
  # summary(as.vector(ycal_mat - y_emu))
  # 
  # y_emu2 <- dat_pde[[n_train + no_input_new]]
  # input_cal - pde_para[n_train + no_input_new,]
  # summary(as.vector(ycal_mat - y_emu2))
  
  # Using predicted emulation data as true y
  # ycal_mat <- c()
  # ycal_ls <- list()
  # set y as time by location matrix
  # for (i in 1:length(res_pre_exact)) {
  #   ycal_mat <- rbind(ycal_mat, t(res_pre_exact[[i]][no_input_new,]))
  #   # ycal_ls[[i]] <- t(res_pre_exact[[i]][no_input_new,])
  # }
  
  
  # simulate field data z ----
  para_gen_cal <- list(n0 = 30, d0 = 2, b = 1, u0 = rep(0, N_sp_cal))
  
  dist_sp_cal <- as.matrix(stats::dist(loc_cal, method = "euclidean", diag = T, upper = T))
  phi_sp_gen_cal <- 3 / (gp_tune * max(dist_sp_cal))
  U_gen <- gen_gp_kernel(loc = loc_cal, phi = phi_sp_gen_cal, sigma2 = gp_sigma2, tau2 = gp_tau2)
  para_gen_cal$rho <- phi_sp_gen_cal
  
  df_gen <- gen_calibrate_data(ycal_mat = ycal_mat, para_gen_cal = para_gen_cal, 
                               U_gen = U_gen)
  tau2_gen <- df_gen$tau2
  z <- df_gen$z
  
  # plot heat map of z and u
  # dat <- z
  tstamp <- as.integer(seq(1, nT_cal, length.out = 9))
  # max_y = max(z)
  max_y <- max(as.vector(unlist(dt_pde_test)))
  cal_heat <- plot_panel_heatmap_9_cal(dat = z, tstamp = tstamp, max_y = max_y,
                                       loc_cal = loc_cal, Nx = Nx, Ny = Ny)
  u_heat <- plot_panel_heatmap_9_cal(dat = df_gen$u, tstamp = tstamp, max_y = max_y,
                                     loc_cal = loc_cal, Nx = Nx, Ny = Ny)
  
  # split to train and test
  prop_train_cal <- 1
  N_sp_train <- round(N_sp_cal * prop_train_cal)
  N_sp_test <- N_sp_cal - N_sp_train
  ind_temp <- sample(1:N_sp_cal)
  ind_train_cal <- ind_temp[1:N_sp_train]
  ind_test_cal <- ind_temp[(N_sp_train+1):N_sp_cal]
  z_train <- z[,ind_train_cal]
  z_test <- z[,ind_test_cal]
  
  loc_cal_train <- loc_cal[ind_train_cal,]
  loc_cal_test <- loc_cal[ind_test_cal,]
  
  # prepare y | Y_{1:T}
  # values_posty <- list()
  # V-1(Y-F*Theta)
  Y_FTheta_EP <- list()
  for (t in 1:nT) {
    Y_FTheta_EP[[t]] <- array(dim = c(dim(Y[[1]]), nsam))
    for (i in 1:nsam) {
      Y_FTheta_EP[[t]][,,i] <- Y[[t]] - F_ls_train[[t]] %*% res_ffbs[[t]][,,i]
    }
  }
  Y_FTheta <- recover_from_EP_MC(dat_EP = Y_FTheta_EP, nT_ori = nT_ori, nT = nT, nsam = nsam)
  
  
  Vinv <- solve(V_ls)
  # VinvY_FTheta_ls <- list()
  t_VinvY_FTheta_ls <- list()
  for (t in 1:nT_ori) {
    # VinvY_FTheta_ls[[t]] <- array(dim = dim(Y_FTheta[[1]]))
    t_VinvY_FTheta_ls[[t]] <- array(dim = c(dim(Y_FTheta[[1]])[2], 
                                            dim(Y_FTheta[[1]])[1],
                                            dim(Y_FTheta[[1]])[3]))
    for (i in 1:nsam) {
      # VinvY_FTheta_ls[[t]][,,i] <- Vinv %*% Y_FTheta[[t]][,,i]
      # t_VinvY_FTheta_ls[[t]][,,i] <- t(VinvY_FTheta_ls[[t]][,,i])
      t_VinvY_FTheta_ls[[t]][,,i] <- crossprod(Y_FTheta[[t]][,,i], Vinv)
    }
  }
  R_cal_train <- gen_gp_kernel(loc = loc_cal_train, phi = phi_sp, sigma2 = gp_sigma2, tau2 = gp_tau2)
  ind_block_cal <- generate.grid.rowsnake(fnrow = 1, fncol = N_sp_cal, bnrow = 1, bncol = bncol)
  n_b_cal <- N_sp_cal / bncol
  values_post_y <- list(t_VinvY_FTheta_ls = t_VinvY_FTheta_ls,
                        R_cal_train = R_cal_train,
                        ind_block_cal = ind_block_cal,
                        n_b_cal = n_b_cal,
                        Vinv = Vinv)
}

# prior ----
{
  # library(mcmc)
  # library(mcmcse)
  # generate prior
  prior_mcmc <- list(
    # eta = input_cal,
    # eta = scale_back_uniform(scale_uniform(input_cal, eta_low, eta_high) * 2,
    #                          low = eta_low, high = eta_high),
    eta = c(6.3, 2.25, 0.0, 1.45, 0.0), # 6.5, 2.3, 1.5
    rho = para_gen_cal$rho,
    n0 = para_gen_cal$n0,
    d0 = para_gen_cal$d0,
    b = para_gen_cal$b,
    # e_metrop = 0.5 * abs(LaplacesDemon::logit(scale_uniform(input_cal, eta_low, eta_high))),
    e_metrop = 0.3 * abs(eta_high - eta_low),
    m0 = rep(0, N_sp_train),
    M0 = diag(N_sp_train))
  prior_mcmc$e_metrop[is.infinite(prior_mcmc$e_metrop)] <- 0
  # TODO: NOT adjust eta equaling to 0
  # prior_mcmc$e_metrop[c(3,5)] <- 0
  
  dist_sp_train <- as.matrix(stats::dist(loc_cal_train, method = "euclidean", diag = T, upper = T))
  # phi_sp_train <- 3 / (gp_tune * max(dist_sp_train))
  U_cal <- gen_gp_kernel(loc = loc_cal_train, phi = prior_mcmc$rho, sigma2 = gp_sigma2, tau2 = gp_tau2)
  prior_mcmc$U <- U_cal
  prior_mcmc$Uinv <- solve(U_cal)
  u_tau2_prior <- gen_prior_u_tau2(n0 = prior_mcmc$n0,
                                   d0 = prior_mcmc$d0,
                                   b = prior_mcmc$b,
                                   m0 = prior_mcmc$m0,
                                   M0 = prior_mcmc$M0,
                                   U = prior_mcmc$U,
                                   nT_cal = nT_cal)
  prior_mcmc$tau2_0 <- u_tau2_prior$tau2_0
  prior_mcmc$tau2 <- u_tau2_prior$tau2
  prior_mcmc$u0 <- u_tau2_prior$u0
  prior_mcmc$u <- u_tau2_prior$u
  
  # TODO should update_y_eta_one
  y_prior <- update_y_eta(eta = prior_mcmc$eta, nsam = nsam, 
                          Nx = Nx, Ny = Ny, N_people = N_people, nT_ori = nT_ori, nT = nT, 
                          nT_cal = nT_cal, N_sp = N_sp, N_sp_train = N_sp_train,
                          AR_choice = AR_choice, bncol = bncol, res_ffbs = res_ffbs, 
                          n_b_cal = values_post_y$n_b_cal,
                          ind_block_cal = values_post_y$ind_block_cal, 
                          t_VinvY_FTheta_ls = values_post_y$t_VinvY_FTheta_ls,
                          Vinv = values_post_y$Vinv, R_cal_train = values_post_y$R_cal_train,
                          input = input, input_cal = input_cal, phi_para = phi_para, 
                          gp_sigma2 = gp_sigma2, gp_tau2 = gp_tau2,
                          ind_train_cal = ind_train_cal)
  prior_mcmc$y <- y_prior$y
  prior_mcmc$mu <- y_prior$mu
  prior_mcmc$Sigma <- y_prior$Sigma
  I_Stilde <- diag(N_sp_train)
}

# # debug
# hist(z_train)
# hist(y_prior$y)
# hist(y_prior$mu)
# hist(y_prior$Sigma)
# plot_panel_heatmap_9_cal(dat = z_train, tstamp = tstamp, max_y = max_y,
#                          loc_cal = loc_cal, Nx = Nx, Ny = Ny)
# plot_panel_heatmap_9_cal(dat = apply(y_prior$y, c(1, 2), mean), tstamp = tstamp, max_y = max_y,
#                          loc_cal = loc_cal, Nx = Nx, Ny = Ny)
# plot_panel_heatmap_9_cal(dat = apply(y_prior$mu, c(1, 2), mean), tstamp = tstamp, max_y = max_y,
#                          loc_cal = loc_cal, Nx = Nx, Ny = Ny)
# prior_mcmc$mu <- (df_gen$z - df_gen$u) * 0.9
# prior_mcmc$y <- prior_mcmc$mu + rnorm(n = length(prior_mcmc$mu), mean = 0, sd = tau2_0)
# mean(res_ffbs$Sigma)
# y_prior$sigma2_cal



# Metropolis----
# prior_mcmc$e_metrop <- abs(0.1 * LaplacesDemon::logit(scale_uniform(input_cal, eta_low, eta_high)))
{
  n_iter <- nsam
  n_burn <- 0
  mcmc_temp <- array(dim = c(n_iter, 8))
  colnames(mcmc_temp) <- c("p", "log_p", "y_new", "y_old", "z_new", "z_old", "jac_old", "jac_new")
  # mcmc_temp <- array(dim = c(n_iter, 6))
  # colnames(mcmc_temp) <- c("p", "log_p", "y_new", "y_old", "jac_old", "jac_new")
  mcmc_accept <- array(dim = c(n_iter))
  mcmc_eta <- array(dim = c(n_iter, n_eta))
  mcmc_tau2 <- array(dim = c(n_iter, nT_cal))
  mcmc_u <- array(dim = c(n_iter, nT_cal, N_sp_train))
  mcmc_y <- array(dim = c(n_iter, nT_cal, N_sp_train))
  n <- 1; i <- 1
}

for (n in 1:(n_burn + n_iter)) {
  if (n %% min(1000, round((n_burn + n_iter)*0.05)) == 0) {
    print(paste("iteration", n, "/", (n_burn + n_iter)))
    print(Sys.time())
  }
  ## 0. initialize values
  if (n == 1){
    n_tilde <- N_sp_train
    eta <- prior_mcmc$eta
    rho <- prior_mcmc$rho
    U <- prior_mcmc$U
    Uinv <- prior_mcmc$Uinv
    n0 <- prior_mcmc$n0
    d0 <- prior_mcmc$d0
    b <- prior_mcmc$b
    m0 <- prior_mcmc$m0
    M0 <- prior_mcmc$M0
    e_metrop <- prior_mcmc$e_metrop
    tau2_0 <- prior_mcmc$tau2_0
    tau2 <- prior_mcmc$tau2 
    u0 <- prior_mcmc$u0 
    u <- prior_mcmc$u 
    Sigma <- prior_mcmc$Sigma
    invSigma <- inv_chol(Sigma)
    ind_accept <- FALSE
    y <- prior_mcmc$y[,,n]
    mu <- prior_mcmc$mu[,,n]
    # y <- prior_mcmc$y
    # mu <- prior_mcmc$mu
  }
  
  
  # 1. update tau2_t and u_t FF
  A_arry <- array(dim = c(nT_cal, n_tilde, n_tilde))
  Q_arry <- array(dim = dim(A_arry))
  m_arry <- array(dim = c(nT_cal, length(m0)))
  M_arry <- array(dim = c(nT_cal, dim(M0)))
  for (t in 1:nT_cal) {
    zt_yt <- z_train[t,] - y[t,]
    # 1.1 update tau2_t
    if(t == 1){
      alpha_tau2 <- b^t * n0 + n_tilde
      q1 <- u[t,] - u0
      Q1 <- crossprod(q1, Uinv) %*% q1
      q2 <- zt_yt - u[t,]
      Q2 <- crossprod(q2, q2)
      beta_tau2 <- b^t * d0 + 1/2 * (Q1 + Q2)
    } else{
      alpha_tau2 <- b^t * n0 + n_tilde
      q1 <- u[t,] - u[(t - 1),]
      Q1 <- crossprod(q1, Uinv) %*% q1
      q2 <- zt_yt - u[t,]
      Q2 <- crossprod(q2, q2)
      beta_tau2 <- b^t * d0 + 1/2 * (Q1 + Q2)
    }
    tau2[t] <- LaplacesDemon::rinvgamma(n = 1, shape = alpha_tau2, scale = beta_tau2)

    # 1.2 updtate ut FF
    if (t == 1) {
      A_arry[t,,] <- tau2_0 * M0 + tau2[t] * U
      Q_arry[t,,] <- A_arry[t,,] + tau2[t] * I_Stilde
      invQt <- inv_chol(Q_arry[t,,])
      AtQtinv <- A_arry[t,,] %*% invQt
      m_arry[t,] <- m0 + AtQtinv %*% 
        (z_train[t,] - y[t,] - m0)
      M_arry[t,,] <- 1/tau2[t] * (A_arry[t,,] - AtQtinv %*% A_arry[t,,])
    } else{
      A_arry[t,,] <- tau2[t-1] * M_arry[(t-1),,] + tau2[t] * U
      Q_arry[t,,] <- A_arry[t,,] + tau2[t] * I_Stilde
      invQt <- inv_chol(Q_arry[t,,])
      AtQtinv <- A_arry[t,,] %*% invQt
      m_arry[t,] <- m_arry[(t-1),] + AtQtinv %*% 
        (z_train[t,] - y[t,] - m_arry[(t-1),])
      M_arry[t,,] <- 1/tau2[t] * (A_arry[t,,] - AtQtinv %*% A_arry[t,,])
    }
  }
  # debug
  print(tau2)
  
  # 2. ut BS
  h_arry <- array(dim = c(nT_cal, length(m0)))
  H_arry <- array(dim = c(nT_cal, dim(M0)))
  # when t = nT_cal
  h_arry[nT_cal,] <- m_arry[nT_cal,]
  H_arry[nT_cal,,] <- tau2[nT_cal] * M_arry[nT_cal,,]
  u[nT_cal,] <- mniw::rmNorm(n = 1, mu = h_arry[nT_cal,], Sigma = H_arry[nT_cal,,])
  for (t in (nT_cal-1):0) {
    if (t == 0) {
      Atp1inv <- inv_chol(A_arry[(t+1),,])
      tau2tMtAtp1inv <- tau2_0 * M0 %*% Atp1inv 
      h0 <- as.vector(m0 + tau2tMtAtp1inv %*% (h_arry[(t+1),] - m0))
      H0 <- tau2_0 * M0 - tau2tMtAtp1inv %*% 
        tcrossprod((A_arry[(t+1),,] - H_arry[(t+1),,]), tau2tMtAtp1inv)
      u0 <- mniw::rmNorm(n = 1, mu = h0, Sigma = H0)
    } else{
      Atp1inv <- inv_chol(A_arry[(t+1),,])
      tau2tMtAtp1inv <- tau2[t] * M_arry[t,,] %*% Atp1inv 
      h_arry[t,] <- as.vector(m_arry[t,] + tau2tMtAtp1inv %*% (h_arry[(t+1),] - m_arry[t,]))
      H_arry[t,,] <- tau2[t] * M_arry[t,,] - tau2tMtAtp1inv %*% 
        tcrossprod((A_arry[(t+1),,] - H_arry[(t+1),,]), tau2tMtAtp1inv)
      u[t,] <- mniw::rmNorm(n = 1, mu = h_arry[t,], Sigma = H_arry[t,,])
    }
  }
  # # debug
  # hist(as.vector(u))
  
  ## 3. update yt
  zt_ut <- array(dim = dim(z_train))
  Bt_arry <- array(dim = c(nT_cal, n_tilde, n_tilde))
  Btbt_arry <- array(dim = c(nT_cal, n_tilde))
  for (t in 1:nT_cal) {
    zt_ut[t,] <- z_train[t,] - u[t,]
    # update yt, Bt, Btbt
    Btbt_s <- cal_Bt_bt(mut = mu[t,], invSigma = invSigma, tau2t = tau2[t], 
                        zt_ut = zt_ut[t,], I_Stilde = I_Stilde)
    Bt_arry[t,,] <- Btbt_s$Bt
    Btbt_arry[t,] <- Btbt_s$Btbt
    y[t,] <- mniw::rmNorm(n = 1, mu = Btbt_arry[t,], Sigma = Bt_arry[t,,])
  }
  # # debug
  # hist(as.vector(y))
  
  ## 4. Metroplis for eta----
  ind_accept <- FALSE
  rand_walk <- mniw::rmNorm(n = 1, mu = rep(0, n_eta), Sigma = e_metrop * diag(n_eta))
  eta_trans <- scale_uniform(eta, low = eta_low, high = eta_high)
  temp_walk <- expit(LaplacesDemon::logit(eta_trans) + rand_walk)
  eta_new <- scale_back_uniform(temp_walk, eta_low, eta_high)
  
  # generate new mu and Sigma for y | Y_1:T
  muSigma_new <- update_muSigma_eta_one(eta = eta_new, ind_sam = n,
                                        Nx = Nx, Ny = Ny, N_people = N_people, nT_ori = nT_ori, nT = nT, 
                                        nT_cal = nT_cal, N_sp = N_sp, N_sp_train = N_sp_train,
                                        AR_choice = AR_choice, bncol = bncol, res_ffbs = res_ffbs, 
                                        n_b_cal = values_post_y$n_b_cal,
                                        ind_block_cal = values_post_y$ind_block_cal, 
                                        t_VinvY_FTheta_ls = values_post_y$t_VinvY_FTheta_ls,
                                        Vinv = values_post_y$Vinv, R_cal_train = values_post_y$R_cal_train,
                                        input = input, input_cal = input_cal, phi_para = phi_para, 
                                        gp_sigma2 = gp_sigma2, gp_tau2 = gp_tau2,
                                        ind_train_cal = ind_train_cal)
  mu_new <- muSigma_new$mu
  Sigma_new <- muSigma_new$Sigma
  invSigma_new <- inv_chol(Sigma_new)
  
  # get parameter for y | ~
  Bt_arry_new <- array(dim = c(nT_cal, n_tilde, n_tilde))
  Btbt_arry_new <- array(dim = c(nT_cal, n_tilde))
  y_new <- array(dim = dim(y))
  for (t in 1:nT_cal) {
    Btbt_s_new <- cal_Bt_bt(mut = mu_new[t,], invSigma = invSigma_new, tau2t = tau2[t], 
                            zt_ut = zt_ut[t,], I_Stilde = I_Stilde)
    Bt_arry_new[t,,] <- Btbt_s_new$Bt
    Btbt_arry_new[t,] <- Btbt_s_new$Btbt
    y_new[t,] <- mniw::rmNorm(n = 1, mu = Btbt_arry_new[t,], Sigma = Bt_arry_new[t,,])
  }
  # # debug
  # hist(as.vector(y_new))
  
  # likelihood of y
  loglik_y_old <- 0
  loglik_y_new <- 0
  for (t in 1:nT_cal) {
    # y_old
    yt_loglik_old <- mniw::dmNorm(x = y[t,], mu = Btbt_arry[t,], Sigma = Bt_arry[t,,], log = TRUE)
    loglik_y_old <- loglik_y_old + yt_loglik_old
    
    # y_new
    yt_loglik_new <- mniw::dmNorm(x = y_new[t,], mu = Btbt_arry_new[t,], Sigma = Bt_arry_new[t,,], log = TRUE)
    loglik_y_new <- loglik_y_new + yt_loglik_new
  }
  
  # likelihood of z
  loglik_z_old <- 0
  loglik_z_new <- 0
  for (t in 1:nT_cal) {
    for (s in 1:n_tilde) {
      # z_old
      zts_loglik_old <- dnorm(x = z_train[t, s], mean = y[t, s] + u[t, s], sd = sqrt(tau2[t]), log = TRUE)
      loglik_z_old <- loglik_z_old + zts_loglik_old
      # z_new
      zts_loglik_new <- dnorm(x = z_train[t, s], mean = y_new[t, s] + u[t, s], sd = sqrt(tau2[t]), log = TRUE)
      loglik_z_new <- loglik_z_new + zts_loglik_new
    }
  }
  
  # calculate Jacobian matrix: e1 [2, 4] e2: [0.2, .4]; a1-a3 [0, 0.2]
  log_jacobian_old <- cal_jacobian_logit_uniform(eta = eta, eta_low = eta_low, eta_high = eta_high, log = TRUE)
  log_jacobian_new <- cal_jacobian_logit_uniform(eta = eta_new, eta_low = eta_low, eta_high = eta_high, log = TRUE)
  if(is.infinite(log_jacobian_old) || is.infinite(log_jacobian_new)){
    log_jacobian_old <- 0
    log_jacobian_new <- 0
  }
  # prior?
  p_old <- 1
  p_new <- 1
  # calculate probability
  log_p_raw <- (loglik_y_new + loglik_z_new + log(p_new) + log_jacobian_old) -
    (loglik_y_old + loglik_z_old + log(p_old) + log_jacobian_new)
  
  # No jacobian
  # log_p_raw <- (loglik_y_new + loglik_z_new) - (loglik_y_old + loglik_z_old)
  
  # No jacobian and y
  # log_p_raw <- (loglik_z_new) - (loglik_z_old)
  
  # y and jacobian
  # log_p_raw <- (loglik_y_new + log(p_new) + log_jacobian_old) -
  #   (loglik_y_old + log(p_old) + log_jacobian_new)
  
  # only y
  # log_p_raw <- (loglik_y_new + log(p_new)) -
  #   (loglik_y_old + log(p_old))
  
  p_raw <- exp(log_p_raw)
  p_update <- min(1, p_raw)
  
  if(runif(1) <= p_update){
    ind_accept <- TRUE
    eta <- eta_new
    mu <- mu_new
    Sigma <- Sigma_new 
    invSigma <- invSigma_new
    # Bt_arry <- Bt_arry_new 
    # Btbt_arry <- Btbt_arry_new 
    y <- y_new
  }
  
  # save
  if (n > n_burn) {
    mcmc_accept[n - n_burn] <- ind_accept
    mcmc_eta[(n - n_burn),] <- eta
    mcmc_tau2[(n - n_burn),] <- tau2
    mcmc_y[(n - n_burn),,] <- y
    mcmc_u[(n - n_burn),,] <- u
    
    mcmc_temp[(n - n_burn),] <- c(p_raw, log_p_raw, 
                                  loglik_y_new, loglik_y_old, 
                                  loglik_z_new, loglik_z_old,
                                  log_jacobian_old, log_jacobian_new)
    # update U if we update rho
    # U <- gp_sigma2 * exp(-rho * dist_sp_train) + gp_tau2 * (dist_sp_train == 0)
  }
}

# check----
{
  print(input_cal)
  print("eta mean:")
  print(colMeans(mcmc_eta, na.rm = T))
  print("eta final:")
  print(mcmc_eta[n_iter,])
  print("accept rate:")
  print(mean(mcmc_accept, na.rm = T))
  print("y_new - y_old")
  hist(mcmc_temp[,3] - mcmc_temp[,4])
  print("j_old - j_new")
  hist(mcmc_temp[,7] - mcmc_temp[,8])
  print("z_new - z_old")
  hist(mcmc_temp[,5] - mcmc_temp[,6])
  mean((mcmc_temp[,5] - mcmc_temp[,6])>0)
  # plot(x = 1:n_iter, y = mcmc_temp[,3] - mcmc_temp[,4], type ="l")
  # plot(x = 1:n_iter, y = mcmc_temp[,5] - mcmc_temp[,6], type ="l")
  # plot(x = 1:n_iter, y = mcmc_temp[,7] - mcmc_temp[,8], type ="l")
  # View(mcmc_temp)
}

# plot----
{
  start_iter = 1
  
  df_eta <- data.frame(
    Iteration = 1:n_iter,
    eta_1 = mcmc_eta[1:n_iter,1],
    eta_2 = mcmc_eta[1:n_iter,2],
    eta_3 = mcmc_eta[1:n_iter,3],
    eta_4 = mcmc_eta[1:n_iter,4],
    eta_5 = mcmc_eta[1:n_iter,5]
  )
  
  df_eta <- data.frame(
    Iteration = 1:(n_iter-start_iter+1),
    eta_1 = mcmc_eta[start_iter:n_iter,1],
    eta_2 = mcmc_eta[start_iter:n_iter,2],
    eta_3 = mcmc_eta[start_iter:n_iter,3],
    eta_4 = mcmc_eta[start_iter:n_iter,4],
    eta_5 = mcmc_eta[start_iter:n_iter,5]
  )
  
  # Reshape the data for ggplot2
  df_eta_long <- df_eta %>%
    pivot_longer(cols = starts_with("eta"), 
                 names_to = "eta", 
                 values_to = "Value")
  
  perc_eta <- data.frame(
    Eta = c("eta_1", "eta_2", "eta_3", "eta_4", "eta_5"),
    `lower` = c(quantile(df_eta$eta_1, probs = 0.025, na.rm = T),
                quantile(df_eta$eta_2, probs = 0.025, na.rm = T),
                quantile(df_eta$eta_3, probs = 0.025, na.rm = T),
                quantile(df_eta$eta_4, probs = 0.025, na.rm = T),
                quantile(df_eta$eta_5, probs = 0.025, na.rm = T)),
    `middle` = c(quantile(df_eta$eta_1, probs = 0.5, na.rm = T),
                 quantile(df_eta$eta_2, probs = 0.5, na.rm = T),
                 quantile(df_eta$eta_3, probs = 0.5, na.rm = T),
                 quantile(df_eta$eta_4, probs = 0.5, na.rm = T),
                 quantile(df_eta$eta_5, probs = 0.5, na.rm = T)),
    `upper` = c(quantile(df_eta$eta_1, probs = 0.975, na.rm = T),
                quantile(df_eta$eta_2, probs = 0.975, na.rm = T),
                quantile(df_eta$eta_3, probs = 0.975, na.rm = T),
                quantile(df_eta$eta_4, probs = 0.975, na.rm = T),
                quantile(df_eta$eta_5, probs = 0.975, na.rm = T))
  )
  
  # Create the plot using ggplot2
  plot_eta <- ggplot(df_eta_long, aes(x = Iteration, y = Value, color = eta)) +
    geom_line() +
    geom_hline(aes(yintercept = input_cal[1]), colour = "black") +
    geom_hline(aes(yintercept = input_cal[1]), colour = "black") +
    geom_hline(aes(yintercept = input_cal[2]), colour = "black") +
    geom_hline(aes(yintercept = input_cal[3]), colour = "black") +
    geom_hline(aes(yintercept = input_cal[4]), colour = "black") +
    geom_hline(aes(yintercept = input_cal[5]), colour = "black") +
    labs(title = "MCMC eta Traces",
         x = "Iteration",
         y = "eta value",
         color = "eta") +
    theme_bw()
  
  
  ggsave(filename = "plot_eta.png",
         path = path_fig,
         plot = plot_eta,
         device = "png",
         width = 20,
         height = 15,
         units = "cm",
         dpi = 300
  )
  
  # Create individual histograms for each eta
  {
    eta_num <- 1
    plot_eta_1 <- ggplot(df_eta, aes(x = eta_1)) +
      geom_histogram(fill = "orange", color = "black") +
      geom_vline(xintercept = input_cal[eta_num], linewidth = 1.25, colour = "red") +
      geom_vline(xintercept = perc_eta[eta_num, "lower"], linetype = "dashed", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "middle"], linetype = "solid", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "upper"], linetype = "dashed", colour = "black") +
      lims(x = c(eta_low[eta_num]-0.1, eta_high[eta_num]+0.1)) +
      labs(title = "Histogram for eta_1", x = "eta_1", y = "Frequency") +
      theme_bw()
    
    plot_eta_1
    
    
    eta_num <- 2
    plot_eta_2 <- ggplot(df_eta, aes(x = eta_2)) +
      geom_histogram(fill = "lightgreen", color = "black") +
      geom_vline(xintercept = input_cal[eta_num], linewidth = 1.25, colour = "red") +
      geom_vline(xintercept = perc_eta[eta_num, "lower"], linetype = "dashed", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "middle"], linetype = "solid", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "upper"], linetype = "dashed", colour = "black") +
      lims(x = c(eta_low[eta_num]-0.01, eta_high[eta_num]+0.01)) +
      labs(title = paste("Histogram for eta_", eta_num, sep = ""), x = paste("eta_", eta_num, sep = ""), y = "Frequency") +
      theme_bw()
    plot_eta_2
    
    eta_num <- 3
    plot_eta_3 <- ggplot(df_eta, aes(x = eta_3)) +
      geom_histogram(fill = "purple", color = "black") +
      geom_vline(xintercept = input_cal[eta_num], linewidth = 1.25, colour = "red") +
      geom_vline(xintercept = perc_eta[eta_num, "lower"], linetype = "dashed", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "middle"], linetype = "solid", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "upper"], linetype = "dashed", colour = "black") +
      lims(x = c(eta_low[eta_num]-0.01, eta_high[eta_num]+0.01)) +
      labs(title = paste("Histogram for eta_", eta_num, sep = ""), x = paste("eta_", eta_num, sep = ""), y = "Frequency") +
      theme_bw()
    
    eta_num <- 4
    plot_eta_4 <- ggplot(df_eta, aes(x = eta_4)) +
      geom_histogram(fill = "lightyellow", color = "black") +
      geom_vline(xintercept = input_cal[eta_num], linewidth = 1.25, colour = "red") +
      geom_vline(xintercept = perc_eta[eta_num, "lower"], linetype = "dashed", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "middle"], linetype = "solid", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "upper"], linetype = "dashed", colour = "black") +
      lims(x = c(eta_low[eta_num]-0.01, eta_high[eta_num]+0.01)) +
      labs(title = paste("Histogram for eta_", eta_num, sep = ""), x = paste("eta_", eta_num, sep = ""), y = "Frequency") +
      theme_bw()
    
    eta_num <- 5
    plot_eta_5 <- ggplot(df_eta, aes(x = eta_5)) +
      geom_histogram(fill = "orange", color = "black") +
      geom_vline(xintercept = input_cal[eta_num], linewidth = 1.25, colour = "red") +
      geom_vline(xintercept = perc_eta[eta_num, "lower"], linetype = "dashed", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "middle"], linetype = "solid", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "upper"], linetype = "dashed", colour = "black") +
      lims(x = c(eta_low[eta_num]-0.01, eta_high[eta_num]+0.01)) +
      labs(title = paste("Histogram for eta_", eta_num, sep = ""), x = paste("eta_", eta_num, sep = ""), y = "Frequency") +
      theme_bw()
    
    # Arrange the histograms in a panel
    plot_eta_panel <- ggarrange(
      plot_eta_1, plot_eta_2, plot_eta_4,
      # plot_eta_3, plot_eta_4, plot_eta_5,
      ncol = 3, nrow = 1
    )
    
    # plot_eta_panel <- ggarrange(
    #   plot_eta_1, plot_eta_2, 
    #   plot_eta_3, plot_eta_4, plot_eta_5,
    #   ncol = 3, nrow = 2
    # )
    
    ggsave(filename = "plot_eta_panel.png",
           path = path_fig,
           plot = plot_eta_panel,
           device = "png",
           width = 40,
           height = 15,
           units = "cm",
           dpi = 300
    )
  }
  
  # tau2----
  df_tau2 <- data.frame(
    Iteration = 1:n_iter,
    tau2_1 = mcmc_tau2[1:n_iter,1],
    tau2_2 = mcmc_tau2[1:n_iter,2],
    tau2_3 = mcmc_tau2[1:n_iter,3],
    tau2_4 = mcmc_tau2[1:n_iter,4],
    tau2_5 = mcmc_tau2[1:n_iter,5],
    tau2_6 = mcmc_tau2[1:n_iter,6]
  )
  
  df_tau2 <- data.frame(
    Iteration = 1:(n_iter-start_iter+1),
    tau2_1 = mcmc_tau2[start_iter:n_iter,1],
    tau2_2 = mcmc_tau2[start_iter:n_iter,2],
    tau2_3 = mcmc_tau2[start_iter:n_iter,3],
    tau2_4 = mcmc_tau2[start_iter:n_iter,4],
    tau2_5 = mcmc_tau2[start_iter:n_iter,5],
    tau2_6 = mcmc_tau2[start_iter:n_iter,6]
  )
  
  perc_tau2 <- data.frame(
    tau2 = c("tau2_1", "tau2_2", "tau2_3", "tau2_4", "tau2_5", "tau2_6"),
    `lower` = c(quantile(df_tau2$tau2_1, probs = 0.025, na.rm = T),
                quantile(df_tau2$tau2_2, probs = 0.025, na.rm = T),
                quantile(df_tau2$tau2_3, probs = 0.025, na.rm = T),
                quantile(df_tau2$tau2_4, probs = 0.025, na.rm = T),
                quantile(df_tau2$tau2_5, probs = 0.025, na.rm = T),
                quantile(df_tau2$tau2_6, probs = 0.025, na.rm = T)),
    `middle` = c(quantile(df_tau2$tau2_1, probs = 0.5, na.rm = T),
                 quantile(df_tau2$tau2_2, probs = 0.5, na.rm = T),
                 quantile(df_tau2$tau2_3, probs = 0.5, na.rm = T),
                 quantile(df_tau2$tau2_4, probs = 0.5, na.rm = T),
                 quantile(df_tau2$tau2_5, probs = 0.5, na.rm = T),
                 quantile(df_tau2$tau2_6, probs = 0.5, na.rm = T)),
    `upper` = c(quantile(df_tau2$tau2_1, probs = 0.975, na.rm = T),
                quantile(df_tau2$tau2_2, probs = 0.975, na.rm = T),
                quantile(df_tau2$tau2_3, probs = 0.975, na.rm = T),
                quantile(df_tau2$tau2_4, probs = 0.975, na.rm = T),
                quantile(df_tau2$tau2_5, probs = 0.975, na.rm = T),
                quantile(df_tau2$tau2_6, probs = 0.975, na.rm = T))
  )
  
  # Reshape the data for ggplot2
  df_tau2_long <- df_tau2 %>%
    pivot_longer(cols = starts_with("tau2"), 
                 names_to = "tau2", 
                 values_to = "Value")
  
  # Create the plot using ggplot2
  plot_tau2 <- ggplot(df_tau2_long, aes(x = Iteration, y = Value, color = tau2)) +
    geom_line() +
    geom_hline(aes(yintercept = tau2_gen[1]), colour = "black") +
    geom_hline(aes(yintercept = tau2_gen[2]), colour = "black") +
    geom_hline(aes(yintercept = tau2_gen[3]), colour = "black") +
    geom_hline(aes(yintercept = tau2_gen[4]), colour = "black") +
    geom_hline(aes(yintercept = tau2_gen[5]), colour = "black") +
    geom_hline(aes(yintercept = tau2_gen[6]), colour = "black") +
    labs(title = "MCMC tau2 Traces",
         x = "Iteration",
         y = "tau2 value",
         color = "tau2") +
    theme_bw()
  
  
  ggsave(filename = "plot_tau2.png",
         path = path_fig,
         plot = plot_tau2,
         device = "png",
         width = 20,
         height = 15,
         units = "cm",
         dpi = 300
  )
  
  print(input_cal)
  print(colMeans(mcmc_eta, na.rm = T))
  print(mcmc_eta[n_iter,])
  print("accept rate:")
  print(mean(mcmc_accept, na.rm = T))
  # mean((input_cal - mcmc_eta[n_iter,])^2, , na.rm = T)
  
  # table ----
  tbl_eta <- cbind(cbind(input_cal, perc_eta$middle, perc_eta$lower, perc_eta$upper))
  colnames(tbl_eta) <- c("true", "median", "2.5%", "97.5%")
  rownames(tbl_eta) <- c("eta_1", "eta_2", "eta_3", "eta_4", "eta_5")
  tbl_eta
  print(round(tbl_eta[c(1,2,4),], 3))
  
  tbl_tau2 <- cbind(cbind(tau2_gen, perc_tau2$middle, perc_tau2$lower, perc_tau2$upper))
  colnames(tbl_tau2) <- c("true", "median", "2.5%", "97.5%")
  rownames(tbl_tau2) <- c("tau2_1", "tau2_2", "tau2_3", "tau2_4", "tau2_5", "tau2_6")
  tbl_tau2
  print(tbl_tau2[c(1,2,4),], 3)
}

# quantile heat map----
eta_low <- c(quantile(mcmc_eta[,1], 0.025), quantile(mcmc_eta[,2], 0.025), 
             quantile(mcmc_eta[,3], 0.025), quantile(mcmc_eta[,4], 0.025), 
             quantile(mcmc_eta[,5], 0.025))
eta_high <- c(quantile(mcmc_eta[,1], 0.975), quantile(mcmc_eta[,2], 0.975), 
              quantile(mcmc_eta[,3], 0.975), quantile(mcmc_eta[,4], 0.975), 
              quantile(mcmc_eta[,5], 0.975))

y_true <- gen_pde(eta = input_cal, Nx = Nx, Ny = Ny, N = N_people, nT_ori = nT_ori)
y_low <- gen_pde(eta = eta_low, Nx = Nx, Ny = Ny, N = N_people, nT_ori = nT_ori)
y_high <- gen_pde(eta = eta_high, Nx = Nx, Ny = Ny, N = N_people, nT_ori = nT_ori)

plot_heat_true <- plot_panel_heatmap_9_cal(dat = y_true, tstamp = tstamp, max_y = max_y,
                                           loc_cal = loc_cal, Nx = Nx, Ny = Ny)
plot_heat_low <- plot_panel_heatmap_9_cal(dat = y_low, tstamp = tstamp, max_y = max_y,
                                          loc_cal = loc_cal, Nx = Nx, Ny = Ny)
plot_heat_high <- plot_panel_heatmap_9_cal(dat = y_high, tstamp = tstamp, max_y = max_y,
                                           loc_cal = loc_cal, Nx = Nx, Ny = Ny)

ggsave(filename = "plot_heat_true.png",
       path = path_fig,
       plot = plot_heat_true,
       device = "png",
       width = 40,
       height = 15,
       units = "cm",
       dpi = 300
)

ggsave(filename = "plot_heat_low.png",
       path = path_fig,
       plot = plot_heat_low,
       device = "png",
       width = 40,
       height = 15,
       units = "cm",
       dpi = 300
)

ggsave(filename = "plot_heat_high.png",
       path = path_fig,
       plot = plot_heat_high,
       device = "png",
       width = 40,
       height = 15,
       units = "cm",
       dpi = 300
)

# Arrange the histograms in a panel
plot_heat_panel <- ggarrange(
  plot_heat_true[[5]], plot_heat_true[[7]], plot_heat_true[[9]], 
  heat_low, heat_high,
  ncol = 3, nrow = 3
)


ggsave(filename = "plot_heat_panel.png",
       path = path_fig,
       plot = plot_heat_panel,
       device = "png",
       width = 40,
       height = 15,
       units = "cm",
       dpi = 300
)

# # y
# # plot heat map of z
# dat <- z
# tstamp <- as.integer(seq(1, nT_cal, length.out = 9))
# # max_y = max(z)
# max_y <- max(as.vector(unlist(dt_pde_test)))
# cal_heat <- plot_panel_heatmap_9_cal(dat = z, tstamp = tstamp, max_y = max_y,
#                                      loc_cal = loc_cal, Nx = Nx, Ny = Ny)







# plot(x = 1:n_iter, y = mcmc_eta[1:n_iter,1])
# plot(x = 1:n_iter, y = mcmc_eta[1:n_iter,2])
# plot(x = 1:n_iter, y = mcmc_eta[1:n_iter,3])
# plot(x = 1:n_iter, y = mcmc_eta[1:n_iter,4])
# plot(x = 1:n_iter, y = mcmc_eta[1:n_iter,5])
# hist(mcmc_eta[,n_iter,5])
# res1 <- mcmc_eta
# save(mcmc_eta, file = paste(path_data, "/mcmc_eta_50.RData", sep = ""))

# F_eta_Theta_train <- list()
# for (i in 1:length(F_eta_Theta)) {
#   F_eta_Theta_train[[i]] <- F_eta_Theta[[i]][,ind_train_cal,] # TODO 
# }


# try to adapt from emulation prediction function
# res_pre_MC_EP <- FFBS_predict_MC(nsam = nsam, Y = Y, res_ffbs = res_ffbs,
#                                  input = input, input_new = input_new[1,],
#                                  F_ls = F_ls, F_new_ls = F_new_ls,
#                                  nT = nT, gp_tune = gp_tune, gp_sigma2 = gp_sigma2,
#                                  gp_tau2 = gp_tau2, delta = 1.0)
# temp <- F_new_ls
# Ftt <- list()
# for (i in 1:nT) {
#   Ftt[[i]] <- t(F_new_ls[[i]][1,])
# }
# F_new_ls <- Ftt
