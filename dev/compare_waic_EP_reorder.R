## START----
{
# need change: n_input, Nx, Ny, nsam
setwd("C:/Users/wilson/Desktop/X1C/MNIW_20240724/dev")
source("init_lib.r")
source("func_lppd_ana.r")
seed <- 1234
set.seed(seed)
path <- "Dan_pde_nopermute_train_10x10_51_20_1"
path_data <- file.path("..", "data", path)
path_fig <- file.path("..", "figures", path)
# Dan_pde_nopermute_train_10x10_51_20_1
# Dan_pde_nopermute_train_100x100_51_20_1

# read in pde data
ind_old_data <- T
AR_choice <- 2
nsam <- 10
n_input <- 20
nT <- 51
Nx <- 10
Ny <- 10
N_people <- 10000
prop_train <- 0.5
N_sp <- Nx * Ny
n_train <- round(n_input * prop_train)
n_test <- n_input - n_train
bnrow <- n_train #
bnrow_test <- n_test #
bncol <- Ny * 1 #
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
phi_para <- 3 / (0.5 * max(dist_para))
V_para <- gen_exp_kernel(loc = pde_para_train, phi = phi_para, sigma2 = 1.1, tau2 = 10^(-4)) # exponential kernel

ind_sp_epi <- ind_sp[1:bncol, ]
dist_sp <- as.matrix(stats::dist(ind_sp_epi, method = "euclidean", diag = T, upper = T))
phi_sp <- 3 / (0.5 * max(dist_sp))
R <- gen_gp_kernel(loc = ind_sp_epi, phi = phi_sp, sigma2 = 1.1, tau2 = 10^(-4))

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
Y_test <- list()


# set F
if(AR_choice == 1){
  F_ls_train <- gen_F_ls_AR1_EP(nT = nT, n_b = n_b, ind = ind, Y = Y_full)
  F_ls_test <- gen_F_ls_AR1_EP(nT = nT, n_b = n_b, ind = ind_test, Y = Y_test_full)
  print("using AR1")
} else if(AR_choice == 2){
  F_ls_train <- gen_F_ls_AR2_EP(nT = nT, n_b = n_b, ind = ind, Y = Y_full)
  F_ls_test <- gen_F_ls_AR2_EP(nT = nT, n_b = n_b, ind = ind_test, Y = Y_test_full)
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
IS <- diag(S) #
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

lppd_ana <- matrix(nrow = 3, ncol = nT)
rownames(lppd_ana) <- c("MNIW", "MNIG", "I")

tbl_time <- matrix(nrow = 3, ncol = 4)
rownames(tbl_time) <- c("MNIW", "MNIG", "I")
colnames(tbl_time) <- c("FFBS", "Sampling", "Model", "WAIC")
}

# save(Y, file = file.path(paste(path_data, "/Y.RData", sep = "")))
# save(F_ls, file = file.path(paste(path_data, "/F_ls_ar1.RData", sep = "")))
# save(F_ls_test, file = file.path(paste(path_data, "/F_ls_test_ar1.RData", sep = "")))



# MNIW----
{
# FFBS
time_start <- Sys.time()
time_start
para_ffbs_mniw <- FFBS(Y = Y, F_ls = F_ls, G_ls = G_ls,
                      W_ls = W_ls, V_ls = V_ls,
                      m0 = m0, M0 = M0,
                      n0 = n0, D0 = D0,
                      nT = nT, delta = 1.0)
time_end <- Sys.time()
time_end
print(time_end - time_start)
tbl_time[1, 1] <- time_end - time_start

# sampling
time_start <- Sys.time()
time_start
res_ffbs_mniw <- FFBS_sampling(nsam = nsam, para_ffbs = para_ffbs_mniw, 
                               F_ls = F_ls, G_ls = G_ls,
                               nT = nT, delta = 1)
time_end <- Sys.time()
time_end
print(time_end - time_start)
tbl_time[1, 2] <- time_end - time_start

### WAIC----
## loglik
time_start <- Sys.time()
loglik_mniw <- matrix(nrow = nsam, ncol = nT)
for (t in 1:nT) {
  if (t %% round(nT_block*0.1) == 0) {
    print(paste("WAIC", t, "/", nT))
    print(Sys.time())
  }
  for (s in 1:nsam) {
    loglik_mniw[s, t] <- mniw::dMNorm(X = Y[[t]], Lambda = F_ls[[t]] %*% res_ffbs_mniw[[t]][,,s],
                                      SigmaR = V_ls, SigmaC = res_ffbs_mniw$Sigma[,,s], log=TRUE)
  }
}
time_end <- Sys.time()
print(time_end - time_start)
tbl_time[1, 4] <- time_end - time_start

waic_mniw <- waic.matrix(loglik_mniw)
print(waic_mniw)
lppd_mniw <- sum(apply(X = loglik_mniw, MARGIN = 2, FUN = matrixStats::logSumExp) - log(nsam))
pwaic_mniw <- sum(apply(X = loglik_mniw, MARGIN = 2, FUN = var))
waic_man_mniw <- -2 * (lppd_mniw - pwaic_mniw)
for (t in 1:nT) {
  # if (t %% round(0.1*nT) == 0) {
  #   print(paste("WAIC ana", t, "/", nT))
  #   print(Sys.time())
  # }
  lppd_ana[1, t] <- lppd_IW_1t(Yt = Y[[t]], Ft = F_ls[[t]], Vt = V_ls, st = para_ffbs_mniw$bs$st[,,t], St = para_ffbs_mniw$bs$St[,,t],
                               nt = para_ffbs_mniw$ff[[nT]]$nt, Dt = para_ffbs_mniw$ff[[nT]]$Dt)
  # lppd_ana[2, t] <- lppd_IG_1t(Yt = Y[[t]], Ft = F_ls[[t]], Vt = V_ls, st = para_ffbs_mnig$bs$st[,,t], St = make_pds(para_ffbs_mnig$bs$St[,,t]),
  #                              nt = para_ffbs_mnig$ff[[t]]$nt, Dt = para_ffbs_mnig$ff[[t]]$Dt, R = R)
  # lppd_ana[3, t] <- lppd_id_1t(Yt = Y[[t]], Ft = F_ls[[t]], Vt = V_ls, st = para_ffbs_i$bs$st[,,t], St = make_pds(para_ffbs_i$bs$St[,,t]))
}
save(waic_mniw, file = file.path(paste(path_data, "/waic_mniw.RData", sep = "")))
save(loglik_mniw, file = file.path(paste(path_data, "/loglik_mniw.RData", sep = "")))

## GPD----
# G
G_mniw <- 0
for (t in 1:nT) {
  temp <- norm(x = Y[[t]] - F_ls[[t]] %*% para_ffbs_mniw$bs$st[,,t], type = "F")^2
  G_mniw <- G_mniw + temp
}

# P
P_mniw <- 0
P_scalar_mniw <- LaplacesDemon::tr(para_ffbs_mniw$ff[[nT]]$Dt) / (para_ffbs_mniw$ff[[nT]]$nt - 2)
for (i in 1:nT) {
  if (i == 1) {
    trace_St <- 0
  }
  trace_St <- trace_St + LaplacesDemon::tr(para_ffbs_mniw$bs$St[,,i])
}
P_mniw <- P_scalar_mniw * trace_St

# D
D_mniw <- G_mniw + P_mniw


rm(para_ffbs_mniw)
rm(res_ffbs_mniw)
gc()
}

# ff <- para_ffbs_mniw$ff
# bs <- para_ffbs_mniw$bs
# save(ff, file = file.path(paste(path_data, "/ff.RData", sep = "")))
# save(bs, file = file.path(paste(path_data, "/bs.RData", sep = "")))



# MNIG----
{
D0 <- p + S
time_start <- Sys.time()
time_start
para_ffbs_mnig <- FFBS_sigma2R(Y = Y, F_ls = F_ls, G_ls = G_ls,
                                   W_ls = W_ls, V_ls = V_ls,
                                   m0 = m0, M0 = M0,
                                   n0 = n0, D0 = D0,
                                   nT = nT, R = R, delta = 1.0)
time_end <- Sys.time()
time_end
print(time_end - time_start)
tbl_time[2, 1] <- time_end - time_start

# sampling
time_start <- Sys.time()
time_start
res_ffbs_mnig <- FFBS_sampling_sigma2R(nsam = nsam, para_ffbs = para_ffbs_mnig, 
                                  F_ls = F_ls, G_ls = G_ls,
                                  nT = nT, R = R, delta = 1)
time_end <- Sys.time()
time_end
print(time_end - time_start)
tbl_time[2, 2] <- time_end - time_start

### WAIC----
time_start <- Sys.time()
loglik_mnig <- matrix(nrow = nsam, ncol = nT)
for (t in 1:nT) {
  if (t %% round(nT_block*0.1) == 0) {
    print(paste("WAIC", t, "/", nT))
    print(Sys.time())
  }
  for (s in 1:nsam) {
    loglik_mnig[s, t] <- mniw::dMNorm(X = Y[[t]], Lambda = F_ls[[t]] %*% res_ffbs_mnig[[t]][,,s],
                                      SigmaR =  V_ls, SigmaC = res_ffbs_mnig$Sigma[s] * R, log = TRUE)
  }
}
time_end <- Sys.time()
tbl_time[2, 4] <- time_end - time_start

waic_mnig <- waic.matrix(loglik_mnig)
print(waic_mnig)
lppd_mnig <- sum(apply(X = loglik_mnig, MARGIN = 2, FUN = matrixStats::logSumExp) - log(nsam))
pwaic_mnig <- sum(apply(X = loglik_mnig, MARGIN = 2, FUN = var))
waic_man_mnig <- -2 * (lppd_mnig - pwaic_mnig)
for (t in 1:nT) {
  # if (t %% round(0.1*nT) == 0) {
  #   print(paste("WAIC ana", t, "/", nT))
  #   print(Sys.time())
  # }
  # lppd_ana[1, t] <- lppd_IW_1t(Yt = Y[[t]], Ft = F_ls[[t]], Vt = V_ls, st = para_ffbs_mniw$bs$st[,,t], St = make_pds(para_ffbs_mniw$bs$St[,,t]),
  #                              nt = para_ffbs_mniw$ff[[t]]$nt, Dt = make_pds(para_ffbs_mniw$ff[[t]]$Dt))
  lppd_ana[2, t] <- lppd_IG_1t(Yt = Y[[t]], Ft = F_ls[[t]], Vt = V_ls, st = para_ffbs_mnig$bs$st[,,t], St = para_ffbs_mnig$bs$St[,,t],
                               nt = para_ffbs_mnig$ff[[nT]]$nt, Dt = para_ffbs_mnig$ff[[nT]]$Dt, R = R)
  # lppd_ana[3, t] <- lppd_id_1t(Yt = Y[[t]], Ft = F_ls[[t]], Vt = V_ls, st = para_ffbs_i$bs$st[,,t], St = make_pds(para_ffbs_i$bs$St[,,t]))
}
save(waic_mnig, file = file.path(paste(path_data, "/waic_mnig.RData", sep = "")))
save(loglik_mnig, file = file.path(paste(path_data, "/loglik_mnig.RData", sep = "")))


## GPD----
# G
G_mnig <- 0
for (t in 1:nT) {
  temp <- norm(x = Y[[t]] - F_ls[[t]] %*% para_ffbs_mnig$bs$st[,,t], type = "F")^2
  G_mnig <- G_mnig + temp
}

# P
P_mnig <- 0
P_scalar_mnig <- LaplacesDemon::tr(R) * para_ffbs_mnig$ff[[nT]]$Dt / (para_ffbs_mnig$ff[[nT]]$nt - 1)
for (i in 1:nT) {
  if (i == 1) {
    trace_St <- 0
  }
  trace_St <- trace_St + LaplacesDemon::tr(para_ffbs_mnig$bs$St[,,i])
}
P_mnig <- P_scalar_mnig * trace_St

# D
D_mnig <- G_mnig + P_mnig


rm(para_ffbs_mnig)
rm(res_ffbs_mnig)
gc()
}


# I----
{
time_start <- Sys.time()
time_start
para_ffbs_i <- FFBS_I(Y = Y, F_ls = F_ls, G_ls = G_ls,
                             W_ls = W_ls, V_ls = V_ls,
                             m0 = m0, M0 = M0,
                             nT = nT, delta = 1.0)
time_end <- Sys.time()
time_end
print(time_end - time_start)
tbl_time[3, 1] <- time_end - time_start

# sampling
time_start <- Sys.time()
time_start
res_ffbs_i <- FFBS_sampling_I(nsam = nsam, para_ffbs = para_ffbs_i, 
                            F_ls = F_ls, G_ls = G_ls,
                            nT = nT, delta = 1)
time_end <- Sys.time()
time_end
print(time_end - time_start)
tbl_time[3, 2] <- time_end - time_start

### WAIC----
time_start <- Sys.time()
loglik_i <- matrix(nrow = nsam, ncol = nT)
for (t in 1:nT) {
  if (t %% round(nT_block*0.1) == 0) {
    print(paste("WAIC", t, "/", nT))
    print(Sys.time())
  }
  for (s in 1:nsam) {
    loglik_i[s, t] <- mniw::dMNorm(X = Y[[t]], Lambda = F_ls[[t]] %*% res_ffbs_i[[t]][,,s],
                                   SigmaR = V_ls, SigmaC = IS, log = TRUE)
  }
}
time_end <- Sys.time()
tbl_time[3, 4] <- time_end - time_start

waic_i <- waic.matrix(loglik_i)
print(waic_i)
lppd_i <- sum(apply(X = loglik_i, MARGIN = 2, FUN = matrixStats::logSumExp) - log(nsam))
pwaic_i <- sum(apply(X = loglik_i, MARGIN = 2, FUN = var))
waic_man_i <- -2 * (lppd_i - pwaic_i)
for (t in 1:nT) {
  # if (t %% round(0.1*nT) == 0) {
  #   print(paste("WAIC ana", t, "/", nT))
  #   print(Sys.time())
  # }
  # lppd_ana[1, t] <- lppd_IW_1t(Yt = Y[[t]], Ft = F_ls[[t]], Vt = V_ls, st = para_ffbs_mniw$bs$st[,,t], St = make_pds(para_ffbs_mniw$bs$St[,,t]),
  #                              nt = para_ffbs_mniw$ff[[t]]$nt, Dt = make_pds(para_ffbs_mniw$ff[[t]]$Dt))
  # lppd_ana[2, t] <- lppd_IG_1t(Yt = Y[[t]], Ft = F_ls[[t]], Vt = V_ls, st = para_ffbs_mnig$bs$st[,,t], St = make_pds(para_ffbs_mnig$bs$St[,,t]),
  #                              nt = para_ffbs_mnig$ff[[t]]$nt, Dt = para_ffbs_mnig$ff[[t]]$Dt, R = R)
  lppd_ana[3, t] <- lppd_id_1t(Yt = Y[[t]], Ft = F_ls[[t]], Vt = V_ls, st = para_ffbs_i$bs$st[,,t], St = para_ffbs_i$bs$St[,,t])
}
save(waic_i, file = file.path(paste(path_data, "/waic_i.RData", sep = "")))
save(loglik_i, file = file.path(paste(path_data, "/loglik_i.RData", sep = "")))


## GPD----
# G
G_i <- 0
for (t in 1:nT) {
  temp <- norm(x = Y[[t]] - F_ls[[t]] %*% para_ffbs_i$bs$st[,,t], type = "F")^2
  G_i <- G_i + temp
}

# P
P_i <- 0
P_scalar_i <- S
for (i in 1:nT) {
  if (i == 1) {
    trace_St <- 0
  }
  trace_St <- trace_St + LaplacesDemon::tr(para_ffbs_i$bs$St[,,i])
}
P_i <- P_scalar_i * trace_St

# D
D_i <- G_i + P_i
 


rm(para_ffbs_i)
rm(res_ffbs_i)
gc()
}



# Summary tables
{
tbl_time[, 3] <- tbl_time[, 1] + tbl_time[, 2]   

tbl_man <- cbind(c(lppd_mniw, lppd_mnig, lppd_i), c(pwaic_mniw, pwaic_mnig, pwaic_i))
colnames(tbl_man) <- c("lppd", "pwaic")
rownames(tbl_man) <- c("MNIW", "MNIG", "I")

tbl_waic <- cbind(c(waic_mniw$waic, waic_mnig$waic, waic_i$waic), c(waic_man_mniw, waic_man_mnig, waic_man_i))
colnames(tbl_waic) <- c("waic_ana", "waic_man")
rownames(tbl_waic) <- c("MNIW", "MNIG", "I")

tbl_lppd <- cbind(c(lppd_mniw, lppd_mnig, lppd_i), rowSums(lppd_ana))
colnames(tbl_lppd) <- c("lppd", "lppd_ana")
rownames(tbl_lppd) <- c("MNIW", "MNIG", "I")

tbl_gpd <- cbind(c(G_mniw, G_mnig, G_i), c(P_mniw, P_mnig, P_i), c(D_mniw, D_mnig, D_i))
colnames(tbl_gpd) <- c("G", "P", "D")
rownames(tbl_gpd) <- c("MNIW", "MNIG", "I")

tbl_out <- matrix(nrow = 3, ncol = 10)
rownames(tbl_out) <- c("MNIW", "MNIG", "I")
colnames(tbl_out) <- c("lppd", "lppd_ana", "WAIC", "G", "P", "D", "Time", "t_FFBS", "t_Sampling","t_WAIC")
tbl_out[, c(1:2)] <- tbl_lppd[, c(1:2)]
tbl_out[, 3] <- tbl_waic[, 1]
tbl_out[, 4:6] <- tbl_gpd
tbl_out[, 7] <- tbl_time[, 3]
tbl_out[, 8:10] <- tbl_time[, c(1, 2, 4)]

print(round(tbl_out, 1))
write.csv(tbl_out, file = paste(path_fig, "/tbl_out_", round(as.numeric(Sys.time())), ".csv", sep = ""))
}


# lppd waic over time----
{
theme_set(theme_minimal(base_size = 22))
mat_waic <- rbind(waic_mniw$pointwise[,3], waic_mnig$pointwise[,3], waic_i$pointwise[,3])
mat_lppd <- rbind(apply(X = loglik_mniw, MARGIN = 2, FUN = matrixStats::logSumExp) - log(nsam),
                  apply(X = loglik_mnig, MARGIN = 2, FUN = matrixStats::logSumExp) - log(nsam),
                  apply(X = loglik_i, MARGIN = 2, FUN = matrixStats::logSumExp) - log(nsam))

rowSums(mat_waic)
rowSums(mat_lppd)

df_waic <- data.frame(time = 1:nT, MNIW = mat_waic[1,], 
                      MNIG = mat_waic[2,], 
                      I = mat_waic[3,]) %>% 
  pivot_longer(!time, names_to = "model", values_to = "waic")


df_lppd <- data.frame(time = 1:nT, 
                      MNIW = mat_lppd[1,], 
                      MNIG = mat_lppd[2,], 
                      I = mat_lppd[3,]) %>% 
  pivot_longer(!time, names_to = "model", values_to = "lppd")

plot_waic <- df_waic %>% 
  ggplot(aes(x = time, y = waic, group = model, colour = model)) + 
  geom_line() + 
  xlab("episode")  
  # ylim(c(-2*10^4, 8*10^4))
plot_waic

plot_lppd <- df_lppd %>%
  ggplot(aes(x = time, y = lppd, group = model, colour = model)) +
  geom_line() +
  xlab("episode")
plot_lppd

ggsave(filename = paste("plot_WL_EP_", round(as.numeric(Sys.time())), ".png", sep = ""),
       path = path_fig,
       plot = ggarrange(plot_waic, plot_lppd,
                        ncol = 2, nrow = 1,
                        labels = c("WAIC", "lppd"),
                        font.label = list(size = 28),
                        vjust = 1.2,
                        hjust = -4,
                        align = "hv",
                        common.legend = T,
                        legend = "right"

       ),
       device = "png",
       width = 40,
       height = 20,
       units = "cm",
       dpi = 100
)


ggsave(filename = paste("plot_lppd_EP_", round(as.numeric(Sys.time())), ".png", sep = ""),
       path = path_fig,
       plot = plot_lppd,
       device = "png",
       width = 40,
       height = 20,
       units = "cm",
       dpi = 100
)

ggsave(filename = paste("plot_waic_EP_", round(as.numeric(Sys.time())), ".png", sep = ""),
       path = path_fig,
       plot = plot_waic,
       device = "png",
       width = 40,
       height = 20,
       units = "cm",
       dpi = 100
)
}


