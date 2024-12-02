## START----
{
# need change: n_input, Nx, Ny, nsam
setwd("D:/Documents/UCLA/0-Administrative/GSR/MNIW/dev")
source("init_lib.r")
seed <- 1234
set.seed(seed)
path <- "Dan_pde_nopermute_train_100x100_51_20_1"
path_data <- file.path("..", "data", path)
path_fig <- file.path("..", "figures", path)

# read in pde data
ind_old_data <- T
nsam <- 10
n_input <- 20
nT <- 51
Nx <- 100
Ny <- 100
N_people <- 10000
prop_train <- 0.5
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
# V_para <- list()

# set F
F_ls_train <- gen_F_ls_AR1_EP(nT = nT, n_b = n_b, ind = ind, Y = Y_full)
F_ls_test <- gen_F_ls_AR1_EP(nT = nT, n_b = n_b, ind = ind_test, Y = Y_test_full)
print("using AR1")

# F_ls_train <- gen_F_ls_AR2_EP(nT = nT, n_b = n_b, ind = ind, Y = Y_full)
# F_ls_test <- gen_F_ls_AR2_EP(nT = nT, n_b = n_b, ind = ind_test, Y = Y_test_full)
# print("using AR2")

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


# FFBS
## MNIW
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

# sampling
time_start <- Sys.time()
time_start
res_ffbs_mniw <- FFBS_sampling(nsam = nsam, para_ffbs = para_ffbs_mniw, 
                               F_ls = F_ls, G_ls = G_ls,
                               nT = nT, delta = 1)
time_end <- Sys.time()
time_end
print(time_end - time_start)


# MNIG
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

# sampling
time_start <- Sys.time()
time_start
res_ffbs_mnig <- FFBS_sampling_sigma2R(nsam = nsam, para_ffbs = para_ffbs_mnig, 
                                  F_ls = F_ls, G_ls = G_ls,
                                  nT = nT, R = R, delta = 1)
time_end <- Sys.time()
time_end
print(time_end - time_start)

# I
time_start <- Sys.time()
time_start
para_ffbs_i <- FFBS_I(Y = Y, F_ls = F_ls, G_ls = G_ls,
                             W_ls = W_ls, V_ls = V_ls,
                             m0 = m0, M0 = M0,
                             nT = nT, delta = 1.0)
time_end <- Sys.time()
time_end
print(time_end - time_start)

# sampling
time_start <- Sys.time()
time_start
res_ffbs_i <- FFBS_sampling_I(nsam = nsam, para_ffbs = para_ffbs_i, 
                            F_ls = F_ls, G_ls = G_ls,
                            nT = nT, delta = 1)
time_end <- Sys.time()
time_end
print(time_end - time_start)
}

# WAIC----
{
## MNIW
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


## MNIG
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


## I
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

waic_mniw <- waic.matrix(loglik_mniw)
waic_mnig <- waic.matrix(loglik_mnig)
waic_i <- waic.matrix(loglik_i)

print(waic_mniw)
print(waic_mnig)
print(waic_i)
}


{
sum(waic_mniw$pointwise[,3])
sum(waic_mnig$pointwise[,3])
sum(waic_i$pointwise[,3])
# waic_mniw$waic / -2 + waic_mniw$p_waic

# ## lppd mannually: the same as loo package
lppd_mniw <- sum(apply(X = loglik_mniw, MARGIN = 2, FUN = matrixStats::logSumExp) - log(nsam))
pwaic_mniw <- sum(apply(X = loglik_mniw, MARGIN = 2, FUN = var))
waic_man_mniw <- -2 * (lppd_mniw - pwaic_mniw)

lppd_mnig <- sum(apply(X = loglik_mnig, MARGIN = 2, FUN = matrixStats::logSumExp) - log(nsam))
pwaic_mnig <- sum(apply(X = loglik_mnig, MARGIN = 2, FUN = var))
waic_man_mnig <- -2 * (lppd_mnig - pwaic_mnig)

lppd_i <- sum(apply(X = loglik_i, MARGIN = 2, FUN = matrixStats::logSumExp) - log(nsam))
pwaic_i <- sum(apply(X = loglik_i, MARGIN = 2, FUN = var))
waic_man_i <- -2 * (lppd_i - pwaic_i)

tbl_man <- cbind(c(lppd_mniw, lppd_mnig, lppd_i), c(pwaic_mniw, pwaic_mnig, pwaic_i))
colnames(tbl_man) <- c("lppd", "pwaic")
rownames(tbl_man) <- c("MNIW", "MNIG", "I")
tbl_man

tbl_waic <- cbind(c(waic_mniw$waic, waic_mnig$waic, waic_i$waic), c(waic_man_mniw, waic_man_mnig, waic_man_i))
colnames(tbl_waic) <- c("waic_ana", "waic_man")
rownames(tbl_waic) <- c("MNIW", "MNIG", "I")
tbl_waic


## lppd analytical ----
# dmatrixtscS <- function(q, m, M, nu, d, R, log=TRUE) {
#   p <- nrow(q)
#   S <- ncol(q)
#   
#   uM <- chol(M)
#   uR <- chol(R)
#   Qbig <- q - m
#   Qbig <- backsolve(uM, Qbig, transpose=TRUE)
#   Qbig <- backsolve(uR, t(Qbig), transpose=TRUE)
#   
#   logp <- -(p * S/2 + nu) * log(1 + sum(Qbig * Qbig)/(2 * d)) + lgamma(p*S/2 + nu) - lgamma(nu) - 
#     p*S/2 * log(2 * pi * d) - S/2 * LaplacesDemon::logdet(M) - p/2 * LaplacesDemon::logdet(R)
#   
#   if (!log) return(exp(logp))
#   return(logp)
# }

dMTig <- function(Y, m, M, nu, d, R, log=TRUE) {
  p <- nrow(Y)
  S <- ncol(Y)
  Rinv <- solve(R)
  Minv <- solve(M)
  q <- Y - m
  Q <- Rinv %*% t(q) %*% Minv %*% q
  temp <- 1 + 1/(2*d) * LaplacesDemon::tr(Q)
  out <- lgamma(p*S/2 + nu) - lgamma(nu) - p*S/2 * log(2 * pi * d) - 
    S/2 * LaplacesDemon::logdet(M) - p/2 * LaplacesDemon::logdet(R) - 
    (p*S/2 + nu) * log(temp)
  
  if (log == TRUE) {
    return(out)
  } else{
    return(exp(out))
  }
}


lppd_IW_1t <- function(Yt, Ft, Vt, st, St, nt, Dt) {
  # note nu = nt - S + 1; hyper T is not matrix t
  lppd <- mniw::dMT(X = Yt, Lambda = Ft %*% st, SigmaR = Ft %*% St %*% t(Ft) + Vt,
                    nu = nt - S + 1, SigmaC = Dt, log = TRUE)
  return(lppd)
}

lppd_IG_1t_dan <- function(Yt, Ft, Vt, st, St, nt, Dt, R) {
  lppd <- dmatrixtscS(q = Yt, m = Ft %*% st, M = Ft %*% St %*% t(Ft) + Vt, 
                      nu = nt, d = Dt, R = R, log = TRUE)
  return(lppd)
}


lppd_IG_1t <- function(Yt, Ft, Vt, st, St, nt, Dt, R) {
  # note nu = nt - S + 1; hyper T is not matrix t
  lppd <- dMTig(Y = Yt, m = Ft %*% st, M = Ft %*% St %*% t(Ft) + Vt, 
                      nu = nt - S + 1, d = Dt, R = R, log = TRUE)
  return(lppd)
}

lppd_id_1t <- function (Yt, Ft, Vt, st, St) {
  lppd <- mniw::dMNorm(X = Yt, Lambda = Ft %*% st,
                       SigmaR = make_pds(Ft %*% St %*% t(Ft) + Vt), SigmaC = diag(ncol(Yt)), log=TRUE)
  return(lppd)
}

# calculate
lppd_ana <- matrix(nrow = 3, ncol = nT)
rownames(lppd_ana) <- c("MNIW", "MNIG", "I")
for (t in 1:nT) {
  if (t %% round(0.1*nT) == 0) {
    print(paste("WAIC ana", t, "/", nT))
    print(Sys.time())
  }
  lppd_ana[1, t] <- lppd_IW_1t(Yt = Y[[t]], Ft = F_ls[[t]], Vt = V_ls, st = para_ffbs_mniw$bs$st[,,t], St = make_pds(para_ffbs_mniw$bs$St[,,t]),
                              nt = para_ffbs_mniw$ff[[t]]$nt, Dt = make_pds(para_ffbs_mniw$ff[[t]]$Dt))
  lppd_ana[2, t] <- lppd_IG_1t(Yt = Y[[t]], Ft = F_ls[[t]], Vt = V_ls, st = para_ffbs_mnig$bs$st[,,t], St = make_pds(para_ffbs_mnig$bs$St[,,t]),
                               nt = para_ffbs_mnig$ff[[t]]$nt, Dt = para_ffbs_mnig$ff[[t]]$Dt, R = R)
  lppd_ana[3, t] <- lppd_id_1t(Yt = Y[[t]], Ft = F_ls[[t]], Vt = V_ls, st = para_ffbs_i$bs$st[,,t], St = make_pds(para_ffbs_i$bs$St[,,t]))
}
rowSums(lppd_ana) / 10^3

tbl_lppd <- cbind(c(lppd_mniw, lppd_mnig, lppd_i), rowSums(lppd_ana))
colnames(tbl_lppd) <- c("lppd", "lppd_ana")
rownames(tbl_lppd) <- c("MNIW", "MNIG", "I")
print(tbl_lppd)
print(tbl_man)
print(tbl_waic)
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
plot_waic

plot_lppd <- df_lppd %>% 
  ggplot(aes(x = time, y = lppd, group = model, colour = model)) + 
  geom_line() + 
  xlab("episode")
plot_lppd

ggsave(filename = paste("plot_waic_EP_AR2_20_", as.numeric(Sys.time()), ".png", sep = ""),
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
}

# ## lppd MC mean----
# {
#   ## MNIW
#   time_start <- Sys.time()
#   loglik_mcmean_mniw <- matrix(nrow = 1, ncol = nT)
#   for (t in 1:nT) {
#     if (t %% 10 == 0) {
#       print(paste("WAIC", t, "/", nT))
#       print(Sys.time())
#     }
#     loglik_mcmean_mniw[1, t] <- mniw::dMNorm(X = Y[[t]], Lambda = F_ls[[t]] %*% apply(res_ffbs_mniw[[t]], c(1,2), mean),
#                                              SigmaR = V_ls, SigmaC = apply(res_ffbs_mniw$Sigma, c(1,2), mean), log = TRUE)
#   }
#   time_end <- Sys.time()
#   print(time_end - time_start)
#   
#   ## MNIG
#   loglik_mcmean_mnig <- matrix(nrow = 1, ncol = nT)
#   for (t in 1:nT) {
#     if (t %% 10 == 0) {
#       print(paste("WAIC", t, "/", nT))
#       print(Sys.time())
#     }
#     loglik_mcmean_mnig[1, t] <- mniw::dMNorm(X = Y[[t]], Lambda = F_ls[[t]] %*% apply(res_ffbs_mnig[[t]], c(1,2), mean),
#                                              SigmaR = V_ls, SigmaC = mean(res_ffbs_mnig$Sigma) * R, log = TRUE)
#   }
#   
#   
#   ## I
#   loglik_mcmean_i <- matrix(nrow = 1, ncol = nT)
#   for (t in 1:nT) {
#     if (t %% 10 == 0) {
#       print(paste("WAIC", t, "/", nT))
#       print(Sys.time())
#     }
#     loglik_mcmean_i[1, t] <- mniw::dMNorm(X = Y[[t]], Lambda = F_ls[[t]] %*% apply(res_ffbs_i[[t]], c(1,2), mean),
#                                           SigmaR = V_ls, SigmaC = IS, log = TRUE)
#   }
#   
# }
# 
# rowSums(loglik_mcmean_mniw) / 10^3
# rowSums(loglik_mcmean_mnig) / 10^3
# rowSums(loglik_mcmean_i) / 10^3
# 
# 
# # calculate waic based on mean of MC samples
# # lppd_mniw_p <- apply(X = loglik_mniw, MARGIN = 2, FUN = mean)
# # summary(lppd_mniw_p)
# # lppd_mnig <- sum(apply(X = loglik_mnig, MARGIN = 2, FUN = matrixStats::logSumExp))
# # pwaic_mnig <- sum(apply(X = loglik_mnig, MARGIN = 2, FUN = var))
# # 
# # lppd_i <- sum(apply(X = loglik_i, MARGIN = 2, FUN = matrixStats::logSumExp))
# # pwaic_i <- sum(apply(X = loglik_i, MARGIN = 2, FUN = var))
# # 
# # tbl_man <- cbind(c(lppd_mniw, lppd_mnig, lppd_i), c(pwaic_mniw, pwaic_mniw, pwaic_i))
# # colnames(tbl_man) <- c("lppd", "pwaic")
# # rownames(tbl_man) <- c("MNIW", "MNIG", "I")
# # tbl_man
# 
# 
# # check mean difference----
# ## MC
# meandif_mniw <- matrix(nrow = nsam, ncol = nT)
# meandif_mnig <- matrix(nrow = nsam, ncol = nT)
# meandif_i <- matrix(nrow = nsam, ncol = nT)
# for (t in 1:nT) {
#   if (t %% 10 == 0) {
#     print(paste("Mean difference", t, "/", nT))
#     print(Sys.time())
#   }
#   for (s in 1:nsam) {
#     meandif_mniw[s, t] <- norm(F_ls[[t]] %*% res_ffbs_mniw[[t]][,,s] - Y[[t]], type = "F")
#     meandif_mnig[s, t] <- norm(F_ls[[t]] %*% res_ffbs_mnig[[t]][,,s] - Y[[t]], type = "F")
#     meandif_i[s, t] <- norm(F_ls[[t]] %*% res_ffbs_i[[t]][,,s] - Y[[t]], type = "F")
#   }
# }
# 
# mean(meandif_mniw)
# mean(meandif_mnig)
# mean(meandif_i)
# 
# ## MC average
# meandif_mcmean_mniw <- matrix(nrow = 1, ncol = nT)
# meandif_mcmean_mnig <- matrix(nrow = 1, ncol = nT)
# meandif_mcmean_i <- matrix(nrow = 1, ncol = nT)
# for (t in 1:nT) {
#   if (t %% 10 == 0) {
#     print(paste("Mean difference", t, "/", nT))
#     print(Sys.time())
#   }
#   meandif_mcmean_mniw[1, t] <- norm(F_ls[[t]] %*% apply(res_ffbs_mniw[[t]], c(1,2), mean) - Y[[t]], type = "F")
#   meandif_mcmean_mnig[1, t] <- norm(F_ls[[t]] %*% apply(res_ffbs_mnig[[t]], c(1,2), mean) - Y[[t]], type = "F")
#   meandif_mcmean_i[1, t] <- norm(F_ls[[t]] %*% apply(res_ffbs_i[[t]], c(1,2), mean) - Y[[t]], type = "F")
#   
# }
# 
# mean(meandif_mcmean_mniw)
# mean(meandif_mcmean_mnig)
# mean(meandif_mcmean_i)
# 
# 
# # analytical
# meandif_ana_mniw <- matrix(nrow = 1, ncol = nT)
# meandif_ana_mnig <- matrix(nrow = 1, ncol = nT)
# meandif_ana_i <- matrix(nrow = 1, ncol = nT)
# for (t in 1:nT) {
#   if (t %% 10 == 0) {
#     print(paste("Mean difference", t, "/", nT))
#     print(Sys.time())
#   }
#   meandif_ana_mniw[1, t] <- norm(F_ls[[t]] %*% para_ffbs_mniw$bs$st[,,t] - Y[[t]], type = "F")
#   meandif_ana_mnig[1, t] <- norm(F_ls[[t]] %*% para_ffbs_mnig$bs$st[,,t] - Y[[t]], type = "F")
#   meandif_ana_i[1, t] <- norm(F_ls[[t]] %*% para_ffbs_i$bs$st[,,t] - Y[[t]], type = "F")
#   
# }
# 
# mean(meandif_ana_mniw)
# mean(meandif_ana_mnig)
# mean(meandif_ana_i)
# 
# save(res_ffbs_mniw, file = file.path(paste(path_data, "/res_ffbs_mniw.RData", sep = "")))
# save(res_ffbs_mnig, file = file.path(paste(path_data, "/res_ffbs_mnig.RData", sep = "")))
# save(res_ffbs_i, file = file.path(paste(path_data, "/res_ffbs_i.RData", sep = "")))

