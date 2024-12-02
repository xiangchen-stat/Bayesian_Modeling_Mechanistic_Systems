{
# preparation----
# need change: n_input, Nx, Ny, nsam
PC <- "mac"
switch (PC,
        "x1c" = setwd("D:/Documents/UCLA/0-Administrative/GSR/Bayesian-Modeling-Mechanistic-Systems/dev"),
        "mac" = setwd("/Users/xiangchen/Documents/Bayesian-Modeling-Mechanistic-Systems/dev"),
        "xps" = setwd("C:/Users/wilson/Desktop/X1C/MNIW/dev")
)
source("init_lib.r")
seed <- 1234
set.seed(seed)
path <- "Dan_pde_nopermute_train_25x25_26_20_1"
path_data <- file.path("..", "data", path)
path_fig <- file.path("..", "figures", path)
if (!dir.exists(path_fig)) {
  dir.create(path_fig)
}
# Dan_pde_nopermute_train_6x6_6_25_1
# Dan_pde_nopermute_train_10x10_11_25_1
# Dan_pde_nopermute_train_10x10_51_20_1
# Dan_pde_nopermute_train_100x100_51_20_1
# Dan_pde_nopermute_train_20x20_21_25_1
# Dan_pde_nopermute_train_50x50_21_25_1
# D:/Documents/UCLA/0-Administrative/GSR/MNIW/dev
# C:/Users/wilson/Desktop/X1C/MNIW/dev

ind_old_data <- F
AR_choice <- 2
nsam <- 200
n_input <- 20
nT <- 26
nT_ori <- nT
Nx <- 25
Ny <- 25
N_people <- 10000
prop_train <- 0.5
gp_tune <- 0.5
gp_sigma2 <- 1.1
gp_tau2 <- 10^(-4)
N_sp <- Nx * Ny
n_train <- round(n_input * prop_train)
n_test <- n_input - n_train
bnrow <- n_train #
bnrow_test <- n_test #
bncol <- Ny * 5 #
pde_para <- as.matrix(read.csv(file = paste0(path_data, "/pde_para.csv", sep = "")))[1:n_input,]
pde_para_train <- pde_para[1:n_train,]
pde_para_test <- pde_para[(n_train + 1):n_input,]
ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))
ind_sp_EP <- ind_sp

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
# dist_sp <- as.matrix(stats::dist(ind_sp, method = "euclidean", diag = T, upper = T))
# phi_sp <- 3 / (gp_tune * max(dist_sp))
# sigma_sp <- gen_gp_kernel(loc = ind_sp, phi = phi_sp, sigma2 = gp_sigma2, tau2 = gp_tau2)


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
D0 <- diag(S)
delta <- 1
gc()
}

## FFBS----
{
time_start <- Sys.time()
time_start
para_ffbs <- FFBS(Y = Y, F_ls = F_ls, G_ls = G_ls,
                      W_ls = W_ls, V_ls = V_ls,
                      m0 = m0, M0 = M0,
                      n0 = n0, D0 = D0,
                      nT = nT, delta = 1.0)
time_end <- Sys.time()
time_end
print(time_end - time_start)

# sampling
time_start <- Sys.time()
time_start
res_ffbs <- FFBS_sampling(nsam = nsam, para_ffbs = para_ffbs, 
                               F_ls = F_ls, G_ls = G_ls,
                               nT = nT, delta = 1)
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

### WAIC----
if (F) {
time_start <- Sys.time()
loglik_train <- matrix(nrow = nsam, ncol = nT)
for (t in 1:nT_block) {
  if (t %% round(nT_block*0.1) == 0) {
    print(paste("WAIC", t, "/", nT))
    print(Sys.time())
  }
  for (s in 1:nsam) {
    loglik_train[s, t] <- mniw::dMNorm(X = Y[[t]], Lambda = F_ls[[t]] %*% res_ffbs[[t]][,,s],
                                       SigmaR = V_ls, SigmaC = res_ffbs$Sigma[,,s], log=TRUE)
  }
}
# save(loglik_train, file = here("data", "loglik_train_MNIW_500.RData"))
waic_train <- waic.matrix(loglik_train)
print(waic_train)
print(paste("lppd is", round(waic_train$waic / -2 + waic_train$p_waic)))
time_end <- Sys.time()
print(time_end - time_start)


## GPD----
# G
G_MNIW <- 0
for (t in 1:nT) {
  temp <- norm(x = Y[[t]] - F_ls[[t]] %*% para_ffbs$bs$st[,,t], type = "F")^2
  G_MNIW <- G_MNIW + temp
}
print(G_MNIW)

# P
P_MNIW <- 0
P_scalar_MNIW <- LaplacesDemon::tr(para_ffbs$ff[[nT]]$Dt) / (para_ffbs$ff[[nT]]$nt - 2)
for (i in 1:nT) {
  if (i == 1) {
    trace_St <- 0
  }
  trace_St <- trace_St + LaplacesDemon::tr(para_ffbs$bs$St[,,i])
}
P_MNIW <- P_scalar_MNIW * trace_St
print(P_MNIW)

# D
D_MNIW <- G_MNIW + P_MNIW
print(D_MNIW)
}
}


# Plot file ----
source("plot_emulation.R")
