# set up emulation ----
{
# preparation----
# need change: n_input, Nx, Ny, nsam
if (!library("here", logical.return = T)) {
  install.packages("here", dependencies = TRUE)
  library("here", character.only = TRUE)
}
setwd(here())
source(here("R", "init_lib.r"))
seed <- 1234
set.seed(seed)
path <- "calibration"
path_data <- here("data", path)
path_fig <- here("figures", path)
if (!dir.exists(path_fig)) {
  dir.create(path_fig)
}

ind_old_data <- F
AR_choice <- 2
nsam <- 3 * 10^4
# nsam <- 3
n_input <- 50
nT <- 26
nT_ori <- nT
Nx <- 6
Ny <- 6
N_people <- 10000
prop_train <- 0.8
gp_tune <- 0.5
gp_sigma2 <- 1.1
gp_tau2 <- 10^(-4)
N_sp <- Nx * Ny
n_train <- round(n_input * prop_train)
n_test <- n_input - n_train
bnrow <- n_train #
bnrow_test <- n_test #
bncol <- Ny * 2 #
pde_para <- as.matrix(read.csv(file = paste0(path_data, "/pde_para.csv", sep = "")))[1:n_input,]
pde_para_train <- pde_para[1:n_train,]
pde_para_test <- pde_para[(n_train + 1):n_input,]
ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))

if(ind_old_data == TRUE){
  # read old data
  load(file.path(paste(path_data, "/dat_pde.RData", sep = "")))
  load(file.path(paste(path_data, "/dt_pde_train.RData", sep = "")))
  load(file.path(paste(path_data, "/dt_pde_test.RData", sep = "")))
} else{
  # manipulate new data
  dat_pde <- list()
  for (i in 1:n_input) {
    dat_pde[[i]] <- as.matrix(read.csv(file = paste0(path_data, "/pde_solution_", i, ".csv", sep = "")))
    if(i %% 10 == 0){
      print(paste("read in pde_solution", i, "/", n_input))
      print(Sys.time())
    }
  }
  save(dat_pde, file = file.path(paste(path_data, "/dat_pde.RData", sep = "")))
  
  # transform the pde train data
  # read in transformed train data
  dt_pde_train <- list()
  for (i in 1:nT) {
    temp <- c()
    for (j in 1:n_train) {
      temp <- rbind(temp, dat_pde[[j]][i,])
    }
    dt_pde_train[[i]] <- temp
    write.csv(temp,
              file=file.path(path_data, paste0("pde_train_", i, ".csv", sep = "")),
              row.names=FALSE)
  }
  names(dt_pde_train) <- paste0("T", seq(1:nT))
  save(dt_pde_train, file = file.path(paste(path_data, "/dt_pde_train.RData", sep = "")))
  
  # transform the pde test data
  # read in transformed test data
  dt_pde_test <- list()
  for (i in 1:nT) {
    temp <- c()
    for (j in (n_train + 1):n_input) {
      temp <- rbind(temp, dat_pde[[j]][i,])
    }
    dt_pde_test[[i]] <- temp
    write.csv(temp,
              file=file.path(path_data, paste0("pde_test_", i, ".csv", sep = "")),
              row.names=FALSE)
  }
  names(dt_pde_test) <- paste0("T", seq(1:nT))
  save(dt_pde_test, file = file.path(paste(path_data, "/dt_pde_test.RData", sep = "")))
}


# Initialize para for FFBS
# Get covariance matrix
dist_para <- as.matrix(stats::dist(pde_para_train, method = "euclidean", diag = T, upper = T))
phi_para <- 3 / (gp_tune * max(dist_para))
V_para <- gen_exp_kernel(loc = pde_para_train, phi = phi_para, sigma2 = gp_sigma2, tau2 = gp_tau2) # exponential kernel

ind_sp_EP <- ind_sp[1:bncol, ]
dist_sp <- as.matrix(stats::dist(ind_sp_EP, method = "euclidean", diag = T, upper = T))
phi_sp <- 3 / (gp_tune * max(dist_sp))
R <- gen_gp_kernel(loc = ind_sp_EP, phi = phi_sp, sigma2 = gp_sigma2, tau2 = gp_tau2)
# R <- sigma_sp[1:bncol, 1:bncol]

# generate season-episode data
V_para_full <- V_para #
Y_full = dt_pde_train #
Y_test_full <- dt_pde_test #

fnrow_train <- n_train #
fnrow_test <- n_test #
fncol <- Nx * Ny
N <- bnrow
S <- bncol
n_block <- fnrow_train / bnrow * fncol / bncol
nT_block <- n_block * nT
# generate index
ind <- generate.grid.rowsnake(fnrow = fnrow_train, fncol = fncol, bnrow = bnrow, bncol = bncol)
ind_test <- generate.grid.rowsnake(fnrow = fnrow_test, fncol = fncol, bnrow = bnrow_test, bncol = bncol)
n_b <- n_block
n_b_test <- dim(ind_test)[1]
Y <- list()
# Y_test <- list()
# V_para <- list()

# set F
if(AR_choice == 1){
  F_ls_train <- gen_F_ls_AR1_EP(nT = nT_ori, n_b = n_b, ind = ind, Y = Y_full)
  F_ls_test <- gen_F_ls_AR1_EP(nT = nT_ori, n_b = n_b, ind = ind_test, Y = Y_test_full)
  print("using AR1")
} else if(AR_choice == 2){
  F_ls_train <- gen_F_ls_AR2_EP(nT = nT_ori, n_b = n_b, ind = ind, Y = Y_full)
  F_ls_test <- gen_F_ls_AR2_EP(nT = nT_ori, n_b = n_b, ind = ind_test, Y = Y_test_full)
  print("using AR2")
}


# set Y
for (i in 1:nT) {
  temp_Y <- Y_full[[i]]
  for (j in 1:n_b) {
    Y[[(i-1)*n_b + j]] <- temp_Y[ind[j,2]:ind[j,3], ind[j,4]:ind[j,5]]
  }
}


F_ls <- F_ls_train
V_ls <- V_para
N <- bnrow
S <- bncol
nT_ori <- nT
nT <- nT_block
p <- dim(F_ls[[1]])[2]
G_ls <- diag(p)
W_ls <- diag(p)
n0 <- p + 3
m0 <- matrix(1, nrow = p, ncol = S)
M0 <- diag(p)
D0 <- p + S #
delta <- 1
gc()
}
# print("calibration succeed")
## FFBS----
{
time_start <- Sys.time()
time_start
para_ffbs <- FFBS_sigma2R(Y = Y, F_ls = F_ls, G_ls = G_ls,
                          W_ls = W_ls, V_ls = V_ls,
                          m0 = m0, M0 = M0,
                          n0 = n0, D0 = D0,
                          nT = nT, R = R, delta = 1.0)
time_end <- Sys.time()
time_end
print(time_end - time_start)

# sampling
time_start <- Sys.time()
time_start
res_ffbs <- FFBS_sampling_sigma2R(nsam = nsam, para_ffbs = para_ffbs, 
                                  F_ls = F_ls, G_ls = G_ls,
                                  nT = nT, R = R, delta = 1)
time_end <- Sys.time()
time_end
print(time_end - time_start)
# save(res_ffbs, file = file.path(paste(path_data, "/res_ffbs_10by10_fixBSbug.RData", sep = "")))
# load(file.path(paste(path_data, "/res_ffbs_10by10_fixBSbug.RData", sep = "")))


# Prediction----
input <- pde_para_train
input_new <- pde_para_test
n_input_new <- dim(input_new)[1]
F_new_ls <- F_ls_test

res_pre_exact_EP <- FFBS_predict_exact(Y = Y, para_ffbs = para_ffbs,
                                       input = input, input_new = input_new,
                                       F_ls = F_ls, F_new_ls = F_new_ls,
                                       nT = nT, phi_para = phi_para, gp_sigma2 = gp_sigma2,
                                       gp_tau2 = gp_tau2, delta = 1.0)
res_pre_MC_EP <- FFBS_predict_MC(nsam = nsam, Y = Y, res_ffbs = res_ffbs,
                                 input = input, input_new = input_new,
                                 F_ls = F_ls, F_new_ls = F_new_ls,
                                 nT = nT, phi_para = phi_para, gp_sigma2 = gp_sigma2,
                                 gp_tau2 = gp_tau2, delta = 1.0)


### transfer back ----
# from episode season to season only
res_pre_exact <- recover_from_EP_exact(dat_EP = res_pre_exact_EP, nT_ori = nT_ori, nT = nT)
res_pre_MC <- recover_from_EP_MC(dat_EP = res_pre_MC_EP, nT_ori = nT_ori, nT = nT, nsam = nsam)
}


# calibration----
{
source(here("R/func_calibrate.R"))
source(here("R/func_plot.R"))
# inherit value from emulation
loc_cal <- ind_sp # use the same location as emulation
nT_cal <- nT_ori
N_sp_cal <- N_sp # choose no of locations
no_input_new <- n_test # choose which new input using in calibration
input_cal <- input_new[no_input_new,]
eta_true <- input_cal
n_eta <- 5
eta_limit_low = c(2, 0.2, rep(0, 3)) # e1 [2, 4] e2: [0.2, .4] # a1-a3 [0, 0.2]
eta_limit_high = c(4, 0.4, rep(0.2, 3))
# eta_limit_low = c(6, 2.2, 0, 1.4, 0) # e1 [6, 7] e2: [2.2, 2.4] a2 [1.4, 1.6] # a1, a3 [0, 0]
# eta_limit_high = c(7, 2.4, 0, 1.6, 0)
ycal_mat <- gen_pde(eta = input_cal, Nx = Nx, Ny = Ny, N = N_people, 
                     nT_ori = nT_ori)

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
# plot_heat_z <- plot_panel_heatmap_9_cal(dat = z, tstamp = tstamp, max_y = max_y,
#                                      loc_cal = loc_cal, Nx = Nx, Ny = Ny, filename = "plot_heat_z_")
# plot_heat_u <- plot_panel_heatmap_9_cal(dat = df_gen$u, tstamp = tstamp, max_y = max_y,
#                                      loc_cal = loc_cal, Nx = Nx, Ny = Ny, filename = "plot_heat_u_")

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
# generate prior
prior_mcmc <- list(
                   # eta = input_cal,
                   # eta = scale_back_uniform(scale_uniform(input_cal, eta_limit_low, eta_limit_high) * 2,
                   #                          low = eta_limit_low, high = eta_limit_high),
                   eta = c(2.7, 0.32, 0.0, 0.07, 0.0), # 3, 0.3, 0.1
                   # eta = c(6.3, 2.25, 0.0, 1.45, 0.0), # 6.5, 2.3, 1.5
                   rho = para_gen_cal$rho,
                   n0 = para_gen_cal$n0,
                   d0 = para_gen_cal$d0,
                   b = para_gen_cal$b,
                   # e_metrop = 0.5 * abs(LaplacesDemon::logit(scale_uniform(input_cal, eta_limit_low, eta_limit_high))),
                   e_metrop = 0.1 * abs(eta_limit_high - eta_limit_low),
                   m0 = rep(0, N_sp_train),
                   M0 = diag(N_sp_train))
prior_mcmc$e_metrop[is.infinite(prior_mcmc$e_metrop)] <- 0

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

y_prior <- update_y_eta_one(eta = prior_mcmc$eta, ind_sam = 1, 
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
I_T <- diag(nT_cal)


# prior_mcmc$e_metrop <- abs(0.1 * LaplacesDemon::logit(scale_uniform(input_cal, eta_limit_low, eta_limit_high)))
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
mcmc_z <- array(dim = c(n_iter, nT_cal, N_sp_train))
n <- 1; i <- 1
}

# Metropolis----
for (n in 1:(n_burn + n_iter)) {
  if (n %% max(min(1000, round((n_burn + n_iter)*0.05)), 1) == 0) {
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
    # y <- prior_mcmc$y[,,n]
    # mu <- prior_mcmc$mu[,,n]
    y <- prior_mcmc$y
    mu <- prior_mcmc$mu
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
  
  
  ## 4. Metroplis for eta----
  ind_accept <- FALSE
  rand_walk <- mniw::rmNorm(n = 1, mu = rep(0, n_eta), Sigma = e_metrop * diag(n_eta))
  eta_trans <- scale_uniform(eta, low = eta_limit_low, high = eta_limit_high)
  temp_walk <- expit(LaplacesDemon::logit(eta_trans) + rand_walk)
  eta_new <- scale_back_uniform(temp_walk, eta_limit_low, eta_limit_high)
  
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
  log_jacobian_old <- cal_jacobian_logit_uniform(eta = eta, eta_limit_low = eta_limit_low, eta_limit_high = eta_limit_high, log = TRUE)
  log_jacobian_new <- cal_jacobian_logit_uniform(eta = eta_new, eta_limit_low = eta_limit_low, eta_limit_high = eta_limit_high, log = TRUE)
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
    # debug 
    mcmc_temp[(n - n_burn),] <- c(p_raw, log_p_raw, 
                                  loglik_y_new, loglik_y_old, 
                                  loglik_z_new, loglik_z_old,
                                  log_jacobian_old, log_jacobian_new)
    
    # results
    mcmc_accept[n - n_burn] <- ind_accept
    mcmc_eta[(n - n_burn),] <- eta
    mcmc_tau2[(n - n_burn),] <- tau2
    mcmc_y[(n - n_burn),,] <- y
    mcmc_u[(n - n_burn),,] <- u
    # generate post prediction of z
    epsilon <- mniw::rmNorm(n = 1, mu = rep(0, nT_cal), Sigma = tau2 * I_T)
    mcmc_z[(n - n_burn),,] <- y + u + epsilon
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
print(mcmc_eta[n-1,])
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

# PLOT----
{
  ## 0. dataframe----
  start_iter <- 1
  # start_iter <- 1000
  # n_iter <- 11000
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
    `lower` = apply(df_eta[,-1], MARGIN = 2, FUN = quantile, probs = 0.025, na.rm = T),
    `middle` = apply(df_eta[,-1], MARGIN = 2, FUN = quantile, probs = 0.5, na.rm = T),
    `upper` = apply(df_eta[,-1], MARGIN = 2, FUN = quantile, probs = 0.975, na.rm = T)
  )
  rownames(perc_eta) <- NULL
  
  # 1. ETA----
  ## eta trace----
  plot_eta_trace <- ggplot(df_eta_long, aes(x = Iteration, y = Value, color = eta)) +
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
         color = "eta") 
  
  
  # ggsave(filename = "plot_eta_trace.png",
  #        path = path_fig,
  #        plot = plot_eta_trace,
  #        device = "png",
  #        width = 30,
  #        height = 20,
  #        units = "cm",
  #        dpi = 300
  # )
  
  ## eta histogram----
  {
    eta_num <- 1
    plot_eta_1 <- ggplot(df_eta, aes(x = eta_1)) +
      geom_histogram(fill = "orange", color = "black") +
      geom_vline(xintercept = input_cal[eta_num], linewidth = 1, colour = "red") +
      geom_vline(xintercept = perc_eta[eta_num, "lower"], linetype = "dashed", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "middle"], linetype = "solid", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "upper"], linetype = "dashed", colour = "black") +
      lims(x = c(eta_limit_low[eta_num]-0.1, eta_limit_high[eta_num]+0.1)) +
      # labs(x = paste("eta ", eta_num, sep = ""), y = NULL) +
      # labs(title = paste("Histogram for eta ", eta_num, sep = ""), x = paste("eta ", eta_num, sep = ""), y = "Frequency") +
      labs(x = NULL, y = NULL, title = paste("eta ", eta_num, sep = "")) 
    # plot_eta_1
    
    
    eta_num <- 2
    plot_eta_2 <- ggplot(df_eta, aes(x = eta_2)) +
      geom_histogram(fill = "lightgreen", color = "black") +
      geom_vline(xintercept = input_cal[eta_num], linewidth = 1, colour = "red") +
      geom_vline(xintercept = perc_eta[eta_num, "lower"], linetype = "dashed", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "middle"], linetype = "solid", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "upper"], linetype = "dashed", colour = "black") +
      lims(x = c(eta_limit_low[eta_num]-0.01, eta_limit_high[eta_num]+0.01)) +
      # labs(x = paste("eta ", eta_num, sep = ""), y = NULL) +   
      # labs(title = paste("Histogram for eta ", eta_num, sep = ""), x = paste("eta ", eta_num, sep = ""), y = "Frequency") 
      labs(x = NULL, y = NULL, title = paste("eta ", eta_num, sep = "")) 
      # plot_eta_2
      
      eta_num <- 3
    plot_eta_3 <- ggplot(df_eta, aes(x = eta_3)) +
      geom_histogram(fill = "purple", color = "black") +
      geom_vline(xintercept = input_cal[eta_num], linewidth = 1, colour = "red") +
      geom_vline(xintercept = perc_eta[eta_num, "lower"], linetype = "dashed", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "middle"], linetype = "solid", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "upper"], linetype = "dashed", colour = "black") +
      lims(x = c(eta_limit_low[eta_num]-0.01, eta_limit_high[eta_num]+0.01)) +
      # labs(x = paste("eta ", eta_num, sep = ""), y = NULL) +   
      labs(x = NULL, y = NULL, title = paste("eta ", eta_num, sep = "")) 
      # labs(title = paste("Histogram for eta ", eta_num, sep = ""), x = paste("eta ", eta_num, sep = ""), y = "Frequency") 
      
      eta_num <- 4
    plot_eta_4 <- ggplot(df_eta, aes(x = eta_4)) +
      geom_histogram(fill = "lightyellow", color = "black") +
      geom_vline(xintercept = input_cal[eta_num], linewidth = 1, colour = "red") +
      geom_vline(xintercept = perc_eta[eta_num, "lower"], linetype = "dashed", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "middle"], linetype = "solid", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "upper"], linetype = "dashed", colour = "black") +
      lims(x = c(eta_limit_low[eta_num]-0.01, eta_limit_high[eta_num]+0.01)) +
      # labs(x = paste("eta ", eta_num, sep = ""), y = NULL) +   
      labs(x = NULL, y = NULL, title = paste("eta ", eta_num, sep = "")) 
      # labs(title = paste("Histogram for eta ", eta_num, sep = ""), x = paste("eta ", eta_num, sep = ""), y = "Frequency") 
      
      eta_num <- 5
    plot_eta_5 <- ggplot(df_eta, aes(x = eta_5)) +
      geom_histogram(fill = "orange", color = "black") +
      geom_vline(xintercept = input_cal[eta_num], linewidth = 1, colour = "red") +
      geom_vline(xintercept = perc_eta[eta_num, "lower"], linetype = "dashed", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "middle"], linetype = "solid", colour = "black") +
      geom_vline(xintercept = perc_eta[eta_num, "upper"], linetype = "dashed", colour = "black") +
      lims(x = c(eta_limit_low[eta_num]-0.01, eta_limit_high[eta_num]+0.01)) +
      # labs(x = paste("eta ", eta_num, sep = ""), y = NULL) +   
      labs(x = NULL, y = NULL, title = paste("eta ", eta_num, sep = "")) 
      # labs(title = paste("Histogram for eta ", eta_num, sep = ""), x = paste("eta ", eta_num, sep = ""), y = "Frequency") 
      
      # Arrange the histograms in a panel
      plot_eta_panel <- ggarrange(
        plot_eta_1, plot_eta_2, plot_eta_4,
        # plot_eta_3, plot_eta_4, plot_eta_5,
        ncol = 3, nrow = 1,
        # labels = c("eta 1", "eta 2", "eta 4"),
        font.label = list(size = 28),
        vjust = 1,
        hjust = -0.3,
        align = "hv",
        common.legend = T,
        legend = "right"
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
           width = 45,
           height = 15,
           units = "cm",
           dpi = 300
    )
  }
  
  # 2. TAU2----
  n_tau2 <- ncol(mcmc_tau2)
  
  # Create a data frame with dynamic column names
  # df_tau2 <- data.frame(
  #   Iteration = 1:n_iter,
  #   setNames(as.data.frame(mcmc_tau2[1:n_iter, 1:n_tau2]), paste0("tau2_", 1:n_tau2))
  # )
  
  df_tau2 <- data.frame(
    Iteration = 1:(n_iter - start_iter + 1),
    setNames(as.data.frame(mcmc_tau2[start_iter:n_iter, 1:n_tau2]), paste0("tau2_", 1:n_tau2))
  )
  
  perc_tau2 <- data.frame(
    tau2 = paste0("tau2_", 1:n_tau2),
    `true` = tau2_gen,
    `lower` = apply(df_tau2[,-1], MARGIN = 2, FUN = quantile, probs = 0.025, na.rm = T),
    `middle` = apply(df_tau2[,-1], MARGIN = 2, FUN = quantile, probs = 0.5, na.rm = T),
    `upper` = apply(df_tau2[,-1], MARGIN = 2, FUN = quantile, probs = 0.975, na.rm = T),
    `time` = 0:(nT_cal-1)
  )
  rownames(perc_tau2) <- NULL
  
  # Reshape the data for ggplot2
  df_tau2_long <- df_tau2 %>%
    pivot_longer(cols = starts_with("tau2"), 
                 names_to = "tau2", 
                 values_to = "Value")
  
  plot_tau2_panel <- ggplot(df_tau2_long, aes(x = Iteration, y = Value, color = as.factor(tau2))) +
    geom_line(alpha = 0.7) +
    geom_hline(aes(yintercept = Value), data = data.frame(tau2 = unique(df_tau2_long$tau2), Value = tau2_gen), 
               colour = "red") +
    labs(title = "MCMC tau2 Traces",
         x = "Iteration",
         y = "tau2 value",
         color = "tau2") +
    facet_wrap(~ tau2, scales = "free_y") + 
    theme(legend.position = "none")
  
  # ggsave(filename = "plot_tau2_panel.png",
  #        path = path_fig,
  #        plot = plot_tau2_panel,
  #        device = "png",
  #        width = 40,
  #        height = 30,
  #        units = "cm",
  #        dpi = 300
  # )
  
  ## band plot----
  plot_band_tau2 <- ggplot(perc_tau2, aes(x = time)) +
    # Confidence interval band with a softer color and transparency
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#A6CEE3", alpha = 0.6) +
    # Middle line with a more distinctive color and line type
    geom_line(aes(y = middle), color = "#1F78B4", size = 1, linetype = "dotdash") +
    # True line in a bold red for emphasis
    geom_line(aes(y = true), color = "#E31A1C", size = 1, linetype = "solid") +
    labs(
      title = expression("Variance (" * tau^2 * ")"),
      # subtitle = "True values in red, with middle estimate and confidence interval",
      x = "Time",
      y = "Value"
    ) +
    # xlim(c(1,nT_cal))+
    ylim(c(0, max(perc_tau2[,"upper"]))) 
    # Fine-tuning theme elements for clarity
    # Add minor aesthetic tweaks
    # scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    # scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) 
    
    
    # ggsave(filename = "plot_band_tau2.png",
    #        path = path_fig,
    #        plot = plot_band_tau2,
    #        device = "png",
    #        width = 30,
    #        height = 25,
    #        units = "cm",
    #        dpi = 300
    # )
  
  # 3. Y and U----
  ## quantile heat map for y----
  eta_lower <- apply(mcmc_eta, MARGIN = 2, FUN = quantile, probs = 0.025, na.rm = T)
  eta_upper <- apply(mcmc_eta, MARGIN = 2, FUN = quantile, probs = 0.975, na.rm = T)
  
  y_true <- gen_pde(eta = eta_true, Nx = Nx, Ny = Ny, N = N_people, nT_ori = nT_ori)
  y_low <- gen_pde(eta = eta_lower, Nx = Nx, Ny = Ny, N = N_people, nT_ori = nT_ori)
  y_high <- gen_pde(eta = eta_upper, Nx = Nx, Ny = Ny, N = N_people, nT_ori = nT_ori)
  
  plot_heat_y_true <- plot_panel_heatmap_9_cal_nolab(dat = y_true, tstamp = tstamp, max_y = max_y, savei = F,
                                                     loc_cal = loc_cal, Nx = Nx, Ny = Ny, filename = "plot_heat_y_true")
  plot_heat_y_low <- plot_panel_heatmap_9_cal_nolab(dat = y_low, tstamp = tstamp, max_y = max_y, savei = F,
                                                    loc_cal = loc_cal, Nx = Nx, Ny = Ny, filename = "plot_heat_y_low")
  plot_heat_y_high <- plot_panel_heatmap_9_cal_nolab(dat = y_high, tstamp = tstamp, max_y = max_y, savei = F,
                                                     loc_cal = loc_cal, Nx = Nx, Ny = Ny, filename = "plot_heat_y_high")
  
  t_heat <- c(2, 5, 8)
  plot_heat_y_panel <- ggarrange(
    plot_heat_y_true[[t_heat[1]]], plot_heat_y_true[[t_heat[2]]], plot_heat_y_true[[t_heat[3]]], 
    plot_heat_y_low[[t_heat[1]]], plot_heat_y_low[[t_heat[2]]], plot_heat_y_low[[t_heat[3]]], 
    plot_heat_y_high[[t_heat[1]]], plot_heat_y_high[[t_heat[2]]], plot_heat_y_high[[t_heat[3]]], 
    ncol = 3, nrow = 3,
    labels = c(paste("True eta: t =", tstamp[t_heat[1]]-1),
               paste("t =", tstamp[t_heat[2]]-1),
               paste("t =", tstamp[t_heat[3]]-1),
               paste("2.5% eta: t =", tstamp[t_heat[1]]-1),
               paste("t =", tstamp[t_heat[2]]-1),
               paste("t =", tstamp[t_heat[3]]-1),
               paste("97.5% eta: t =", tstamp[t_heat[1]]-1),
               paste("t =", tstamp[t_heat[2]]-1),
               paste("t =", tstamp[t_heat[3]]-1)
    ),
    font.label = list(size = 28),
    vjust = 1.2,
    hjust = -0.3,
    align = "hv",
    common.legend = T,
    legend = "right"
  ) %>% annotate_figure(left = text_grob("Y", rot = 90, vjust = 1, size = 28),
                        bottom = text_grob("X", vjust = 0, size = 28))
  
  ggsave(filename = "plot_heat_y_panel.png",
         path = path_fig,
         plot = plot_heat_y_panel,
         device = "png",
         width = 40,
         height = 35,
         units = "cm",
         dpi = 300
  )
  
  ## quantile heat map for u----
  u_true <- df_gen$u
  u_low <- apply(mcmc_u, MARGIN = c(2,3), FUN = quantile, probs = 0.025, na.rm = T)
  u_high <- apply(mcmc_u, MARGIN = c(2,3), FUN = quantile, probs = 0.975, na.rm = T)
  
  # plot_heat_u_true <- plot_panel_heatmap_9_cal_nolab(dat = u_true, tstamp = tstamp, max_y = max_y,
  #                                                    loc_cal = loc_cal, Nx = Nx, Ny = Ny, filename = "plot_heat_u_true")
  # plot_heat_u_low <- plot_panel_heatmap_9_cal_nolab(dat = u_low, tstamp = tstamp, max_y = max_y,
  #                                                   loc_cal = loc_cal, Nx = Nx, Ny = Ny, filename = "plot_heat_u_low")
  # plot_heat_u_high <- plot_panel_heatmap_9_cal_nolab(dat = u_high, tstamp = tstamp, max_y = max_y,
  #                                                    loc_cal = loc_cal, Nx = Nx, Ny = Ny, filename = "plot_heat_u_high")
  # 
  # t_heat <- c(2, 5, 8)
  # plot_heat_u_panel <- ggarrange(
  #   plot_heat_u_true[[t_heat[1]]], plot_heat_u_true[[t_heat[2]]], plot_heat_u_true[[t_heat[3]]], 
  #   plot_heat_u_low[[t_heat[1]]], plot_heat_u_low[[t_heat[2]]], plot_heat_u_low[[t_heat[3]]], 
  #   plot_heat_u_high[[t_heat[1]]], plot_heat_u_high[[t_heat[2]]], plot_heat_u_high[[t_heat[3]]], 
  #   ncol = 3, nrow = 3,
  #   labels = c(paste("True eta: t =", tstamp[t_heat[1]]-1),
  #              paste("t =", tstamp[t_heat[2]]-1),
  #              paste("t =", tstamp[t_heat[3]]-1),
  #              paste("2.5% eta: t =", tstamp[t_heat[1]]-1),
  #              paste("t =", tstamp[t_heat[2]]-1),
  #              paste("t =", tstamp[t_heat[3]]-1),
  #              paste("97.5% eta: t =", tstamp[t_heat[1]]-1),
  #              paste("t =", tstamp[t_heat[2]]-1),
  #              paste("t =", tstamp[t_heat[3]]-1)
  #   ),
  #   font.label = list(size = 28),
  #   vjust = 1.2,
  #   hjust = -0.3,
  #   align = "hv",
  #   common.legend = T,
  #   legend = "right"
  # ) %>% annotate_figure(left = text_grob("Y", rot = 90, vjust = 1, size = 28),
  #                       bottom = text_grob("X", vjust = 0, size = 28))
  # 
  # ggsave(filename = "plot_heat_u_panel.png",
  #        path = path_fig,
  #        plot = plot_heat_u_panel,
  #        device = "png",
  #        width = 40,
  #        height = 35,
  #        units = "cm",
  #        dpi = 300
  # )
  
  # 4. Z and U----
  ## z dif plot----
  z_temp <- apply(mcmc_z, c(2, 3), mean, na.rm = T)
  z_mean <- z_temp[,order(ind_train_cal)] # indices are randomized, transfer back
  # dif_temp <- z[,ind_train_cal] - z_mean
  dif_temp <- z - z_mean

  ## scatter plot z ----
  mcmc_z_lower <- apply(mcmc_z, MARGIN = c(2, 3), FUN = quantile, probs = 0.025, na.rm = T)
  mcmc_z_med <- apply(mcmc_z, MARGIN = c(2, 3), FUN = quantile, probs = 0.5, na.rm = T)
  mcmc_z_upper <- apply(mcmc_z, MARGIN = c(2, 3), FUN = quantile, probs = 0.975, na.rm = T)
  
  z_reorder <- z[,ind_train_cal]
  # hist(as.vector(z_reorder - mcmc_z_med))
  df_stat_z <- data.frame(true = matrix(data = z_reorder, ncol = 1, byrow = T),
                          lower = matrix(data = mcmc_z_lower, ncol = 1, byrow = T),
                          med = matrix(data = mcmc_z_med, ncol = 1, byrow = T),
                          upper = matrix(data = mcmc_z_upper, ncol = 1, byrow = T))
  # coverage
  # sum(mcmc_z_lower <= z_reorder & z_reorder <= mcmc_z_upper) / length(z_reorder)
  
  
  alpha_error <- 0.5
  error_width <- (max(z) - min(z)) / 30 # denominator is for aesthetic 
  
  plot_scatter_z <- df_stat_z %>% ggplot(aes(x = true, y = med)) + 
    # geom_linerange(aes(ymin = lower, ymax = upper)) 
    geom_pointrange(aes(ymin = lower, ymax = upper), size =.3, alpha = alpha_error)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width, alpha = alpha_error) + 
    geom_abline(col = "red", linewidth = 0.8) + 
    labs(x = "True field value", y = "Posterior prediction")
  # plot_scatter_z
  
  # ggsave(filename = paste("plot_error_", as.numeric(Sys.time()), ".png", sep = ""),
  # ggsave(filename = paste("plot_error_z.png", sep = ""),
  #        path = path_fig,
  #        plot = plot_scatter_z,
  #        device = "png",
  #        width = 30,
  #        height = 25,
  #        units = "cm",
  #        dpi = 300
  # ) 
  
  
  ## scatter u----
  mcmc_u_lower <- apply(mcmc_u, MARGIN = c(2, 3), FUN = quantile, probs = 0.025, na.rm = T)
  mcmc_u_med <- apply(mcmc_u, MARGIN = c(2, 3), FUN = quantile, probs = 0.5, na.rm = T)
  mcmc_u_upper <- apply(mcmc_u, MARGIN = c(2, 3), FUN = quantile, probs = 0.975, na.rm = T)
  
  u_reorder <- u
  hist(as.vector(u_reorder - mcmc_u_med))
  df_stat_u <- data.frame(true = matrix(data = u_reorder, ncol = 1, byrow = T),
                          lower = matrix(data = mcmc_u_lower, ncol = 1, byrow = T),
                          med = matrix(data = mcmc_u_med, ncol = 1, byrow = T),
                          upper = matrix(data = mcmc_u_upper, ncol = 1, byrow = T))
  # coverage
  # sum(mcmc_u_lower <= u_reorder & u_reorder <= mcmc_u_upper) / length(u_reorder)
  
  alpha_error <- 0.5
  error_width <- (max(u_true) - min(u_true)) / 30 # denominator is for aesthetic 
  
  plot_scatter_u <- df_stat_u %>% ggplot(aes(x = true, y = med)) + 
    # geom_linerange(aes(ymin = lower, ymax = upper)) 
    geom_pointrange(aes(ymin = lower, ymax = upper), size =.3, alpha = alpha_error)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width, alpha = alpha_error) + 
    geom_abline(col = "red", linewidth = 0.8) + 
    labs(x = "True spatial bias", y = "Posterior prediction")
  
  # ggsave(filename = paste("plot_error_", as.numeric(Sys.time()), ".png", sep = ""),
  # ggsave(filename = paste("plot_error_u.png", sep = ""),
  #        path = path_fig,
  #        plot = plot_scatter_u,
  #        device = "png",
  #        width = 30,
  #        height = 25,
  #        units = "cm",
  #        dpi = 300
  # ) 
  
  
  ## band plot z ----
  plot_band_predict_z <- df_stat_z %>%
    ggplot(aes(x = true, y = med)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper, x = true), fill = "#A6CEE3", alpha = 1) +
    geom_point(aes(y = med), color = "#1F78B4", size = 0.7, alpha = 1) +
    geom_abline(color = "#E31A1C", linewidth = 1, linetype = "solid") + 
    labs(title = expression("Field Value (" * z * ")"),
         x = "True field value", 
         y = "Posterior prediction"
         )
  
  # ggsave(filename = paste("plot_band_predict_z.png", sep = ""),
  #        path = path_fig,
  #        plot = plot_band_predict_z,
  #        device = "png",
  #        width = 25,
  #        height = 25,
  #        units = "cm",
  #        dpi = 300
  # ) 
  
  ## band plot u----
  plot_band_predict_u <- df_stat_u %>%
    ggplot(aes(x = true, y = med)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper, x = true), fill = "#A6CEE3", alpha = 1) +
    geom_point(aes(y = med), color = "#1F78B4", size = 0.7, alpha = 1) +
    geom_abline(color = "#E31A1C", linewidth = 1, linetype = "solid") + 
    labs(title = expression("Spatial Bias (" * u * ")"),
         x = "True spatial bias", 
         y = "Posterior prediction"
         )
  
  # ggsave(filename = paste("plot_band_predict_u.png", sep = ""),
  #        path = path_fig,
  #        plot = plot_band_predict_u,
  #        device = "png",
  #        width = 25,
  #        height = 25,
  #        units = "cm",
  #        dpi = 300
  # ) 
  
  # 5. panel TAU2 U Z----
  plot_panel_tau_u_z <- ggarrange(
    plot_band_tau2, plot_band_predict_u, plot_band_predict_z,
    ncol = 3, nrow = 1,
    # labels = c("Variance (tau^2)", "Spatial Bias (u)", "Field Value (z)"),
    font.label = list(size = 24),
    # vjust = 1,
    # hjust = -0.3,
    align = "hv",
    common.legend = T,
    legend = "right"
  )
  
  ggsave(filename = "plot_panel_tau_u_z.png",
         path = path_fig,
         plot = plot_panel_tau_u_z,
         device = "png",
         width = 45,
         height = 15,
         units = "cm",
         dpi = 300
  )
  
  
  # table ----
  tbl_eta <- cbind(cbind(input_cal, perc_eta$middle, perc_eta$lower, perc_eta$upper))
  colnames(tbl_eta) <- c("true", "median", "2.5%", "97.5%")
  rownames(tbl_eta) <- c("eta_1", "eta_2", "eta_3", "eta_4", "eta_5")
  tbl_eta
  print(round(tbl_eta[c(1,2,4),], 3))
  cover_eta <- sum(tbl_eta[,3] <= tbl_eta[,1] & tbl_eta[,1] <= tbl_eta[,4]) / length(tbl_eta[,1])
  
  
  tbl_tau2 <- cbind(cbind(tau2_gen, perc_tau2$middle, perc_tau2$lower, perc_tau2$upper))
  colnames(tbl_tau2) <- c("true", "median", "2.5%", "97.5%")
  rownames(tbl_tau2) <- paste0("tau2_", 1:n_tau2)
  tbl_tau2
  print(round(tbl_tau2, 3))
  cover_tau2 <- sum(tbl_tau2[,3] <= tbl_tau2[,1] & tbl_tau2[,1] <= tbl_tau2[,4]) / length(tbl_tau2[,1])
  
  print(paste("Coverage eta:", cover_eta))
  print(paste("Coverage tau2:", cover_tau2))
}


