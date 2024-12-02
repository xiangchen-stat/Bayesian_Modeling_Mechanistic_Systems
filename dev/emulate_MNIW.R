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

# read in pde data----
ind_old_data <- F
AR_choice <- 2
nsam <- 10
n_input <- 20
nT <- 26
nT_ori <- nT
Nx <- 25
Ny <- 25
N_people <- 10000
N_sp <- Nx * Ny
prop_train <- 0.5
gp_tune <- 0.5
gp_sigma2 <- 1.1
gp_tau2 <- 10^(-4)
n_train <- round(n_input * prop_train)
n_test <- n_input - n_train
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

# dist_sp <- as.matrix(stats::dist(ind_sp, method = "euclidean", diag = T, upper = T))
# phi_sp <- 3 / (0.5 * max(dist_sp))
# sigma_sp <- gen_gp_kernel(loc = ind_sp, phi = phi_sp, sigma2 = 1.1, tau2 = 10^(-4))


# generate para
Y = dt_pde_train
Y_test <- dt_pde_test

# set F
if(AR_choice == 1){
  F_ls_train <- gen_F_ls_AR1(nT = nT, Y = Y)
  F_ls_test <- gen_F_ls_AR1(nT = nT, Y = Y_test)
  print("Using AR 1")
} else if(AR_choice == 2){
  F_ls_train <- gen_F_ls_AR2(nT = nT, Y = Y)
  F_ls_test <- gen_F_ls_AR2(nT = nT, Y = Y_test)
  print("Using AR 2")
}

F_ls <- F_ls_train

V_ls <- V_para
N <- n_train
S <- Nx * Ny
D0 = diag(S)
# D0 = sigma_sp
p <- dim(F_ls[[1]])[2]
G_ls <- diag(p)
W_ls <- diag(p)
n0 <- p + 3
m0 <- matrix(1, nrow = p, ncol = S)
M0 <- diag(p)
# M0 <- gen_pd_matrix(p)
delta <- 1

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

## WAIC----
time_start <- Sys.time()
loglik_train <- matrix(nrow = nsam, ncol = nT)
for (t in 1:nT) {
  # if (t %% round(nT*0.1) == 0) {
  #   print(paste("WAIC", t, "/", nT))
  #   print(Sys.time())
  # }
  for (s in 1:nsam) {
    # loglik_train[s, t] <- dmatnorm(A = Y[[t]], M = F_ls_train[[t]] %*% res_ffbs[[t]][,,s], 
    #                                U = V_ls, V = res_ffbs$Sigma[,,s], log = TRUE)
    # loglik_train[s, t] <- dmn(y = Y[[t]], mu = F_ls_train[[t]] %*% res_ffbs[[t]][,,s], 
    #                                 V = V_ls, S = res_ffbs$Sigma[,,s], log = TRUE)
    # loglik_train[s, t] <- dmatrixnorm(X = Y[[t]], M = F_ls_train[[t]] %*% res_ffbs[[t]][,,s],
    #                                U = V_ls, V = res_ffbs$Sigma[,,s], log = TRUE)
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
}

# ttt <- matrix(data = c(0,0,0,10^5), nrow = 2, ncol = 2)
# waic.matrix(ttt)

{
# Prediction----
input <- pde_para_train
input_new <- pde_para_test
n_input_new <- dim(input_new)[1]
F_new_ls <- F_ls_test

# check with test data
res_pre_MC <- FFBS_predict_MC(nsam = nsam, Y = Y, res_ffbs = res_ffbs,
                                    input = input, input_new = input_new,
                                    F_ls = F_ls, F_new_ls = F_new_ls,
                                    nT = nT, phi_para = phi_para, gp_sigma2 = gp_sigma2,
                           gp_tau2 = gp_tau2, delta = 1.0)
res_pre_exact <- FFBS_predict_exact(Y = Y, para_ffbs = para_ffbs,
                                    input = input, input_new = input_new,
                                    F_ls = F_ls, F_new_ls = F_new_ls,
                                    nT = nT, phi_para = phi_para, gp_sigma2 = gp_sigma2,
                                    gp_tau2 = gp_tau2, delta = 1.0)

# integrate out beta and sigma
# mean_arry <- array(dim = c(dim(res_pre$T1)[1:2], nT))
# dif_arry <- array(dim = c(dim(res_pre$T1)[1:2], nT))
# for (i in 1:nT) {
#   sum_mean <- 0
#   for (j in 1:nsam) {
#     sum_mean <- sum_mean + res_pre[[i]][,,j]
#   }
#   avg_mean <- sum_mean / nsam
#   mean_arry[,,i] <- avg_mean
#   dif_arry[,,i] <- avg_mean - Y_test[[i]]
# }
postm_ls <- list()
dif_ls <- list()
for (i in 1:nT) {
  sum_mean <- 0
  for (j in 1:nsam) {
    sum_mean <- sum_mean + res_pre[[i]][,,j]
  }
  avg_mean <- sum_mean / nsam
  postm_ls[[i]] <- avg_mean
  dif_ls[[i]] <- avg_mean - Y_test[[i]]
}

# plot----
source("../R/func_plot.R")
# Set ggplot theme
theme_set(theme_minimal(base_size = 22))
col_epa <- c("#00e400", "#ffff00", "#ff7e00", "#ff0000", "#99004c", "#7e0023")
col_bgr <- c("#d5edfc", "#a5d9f6", "#7eb4e0", "#588dc8", "#579f8b", "#5bb349",
             "#5bb349", "#f3e35a", "#eda742", "#e36726", "#d64729", "#c52429",
             "#a62021", "#871b1c")

## Plot PDE results ----
# ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))
# plot_ls <- list() 
# tstamp <- as.integer(seq(1, nT, length.out = 9))
# 
# plot_panel_heatmap_9 <- function(dat, input_num, tstamp, max_y, Nx, Ny, nT){
#   plot_ls <- list()
#   # tstamp <- as.integer(seq(1, nT, length.out = 10))
#   ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))
#   ind_plot <- 1
# 
#   for (i in tstamp) {
#     if(i == 1){
#       ind_plot <- 1 # start counting
#     }
#     temp <- dat[[i]][input_num,]
#     rownames(temp) <- NULL
#     colnames(temp) <- NULL
#     dt <- data.frame(row = ind_sp$row, col = ind_sp$col, sol =  temp)%>% 
#       as.data.frame()
#     
#     p <- ggplot(dt, aes(x = col, y = row, fill = sol)) +
#       geom_raster() +
#       scale_fill_gradientn(colours = col_bgr,
#                            limits = c(0, max_y),
#                            oob = scales::squish) +
#       labs(x = "x", y = "y", fill = "Value")+
#       # scale_x_continuous(limits = c(-123.8, -114.2), expand = c(0, 0)) +
#       # scale_y_continuous(limits = c(32.15, 42.04), expand = c(0, 0)) +
#       theme(text = element_text(size=28),
#             legend.text = element_text(size = 28),
#             legend.key.size = unit(1.5, "cm"))
#     
#     plot_ls[[ind_plot]] <- p
#     ind_plot <- ind_plot + 1
#   }
#   
#   ggsave(filename = paste("plot_panel_", as.numeric(Sys.time()), ".png", sep = ""),
#          path = path_fig,
#          plot = ggarrange(plot_ls[[1]], plot_ls[[2]], plot_ls[[3]],
#                           plot_ls[[4]], plot_ls[[5]], plot_ls[[6]],
#                           plot_ls[[7]], plot_ls[[8]], plot_ls[[9]],
#                           ncol = 3, nrow = 3,
#                           # labels = c("(1) Data", "(2) Gaussian Process", "(3) Bilinear", 
#                           #            "(4) MBI", "(5) MBA", "(6) Graphical",
#                           #            "(7) Inverse Graphical", "(8) Nearest Neighbor"),
#                           font.label = list(size = 28),
#                           vjust = 1.2,
#                           # hjust = -1,
#                           align = "hv",
#                           common.legend = T,
#                           legend = "right"
#                           
#          ),
#          device = "png",
#          width = 60,
#          height = 50,
#          units = "cm",
#          dpi = 100
#   )
#   return(plot_ls)
# }
# 
# cal_errorbar <- function(X){
#   out <- data.frame(med = apply(X = X, MARGIN = 1, FUN = median),
#                     lower = apply(X = X, MARGIN = 1, FUN = quantile, prob = 0.025),
#                     upper = apply(X = X, MARGIN = 1, FUN = quantile, prob = 0.975))
#   return(out)
# }
# 
# cal_errorbar_mean <- function(X){
#   out <- data.frame(med = apply(X = X, MARGIN = 1, FUN = mean),
#                     lower = apply(X = X, MARGIN = 1, FUN = quantile, prob = 0.025),
#                     upper = apply(X = X, MARGIN = 1, FUN = quantile, prob = 0.975))
#   return(out)
# }

{
input_num <- 1
tstamp <- as.integer(seq(1, nT, length.out = 9))
# tstamp <- tstamp - 1
# tstamp[[1]] <- 1
dat <- dt_pde_test
max_y <- max(as.vector(unlist(dat))) # set max limit for all plots

pde_heat <- plot_panel_heatmap_9(dat = dat, tstamp = tstamp,
                    input_num = input_num, max_y = max_y, Nx = Nx, Ny = Ny, nT = nT)

# dat <- postm_ls # this is monte carlo
dat <- res_pre_exact # this is exact
ffbs_heat <- plot_panel_heatmap_9(dat = dat, tstamp = tstamp,
                    input_num = input_num, max_y = max_y, Nx = Nx, Ny = Ny, nT = nT)


ggsave(filename = paste("plot_panel_compare_", as.numeric(Sys.time()), ".png", sep = ""),
       path = path_fig,
       plot = ggarrange(pde_heat[[2]], pde_heat[[3]], pde_heat[[4]], 
                        ffbs_heat[[2]], ffbs_heat[[3]], ffbs_heat[[4]], 
                        ncol = 3, nrow = 2,
                        labels = c(paste("PDE: t =", tstamp[2]), 
                                   paste("PDE: t =", tstamp[3]), 
                                   paste("PDE: t =", tstamp[4]), 
                                   paste("FFBS: t =", tstamp[2]), 
                                   paste("FFBS: t =", tstamp[3]), 
                                   paste("FFBS: t =", tstamp[4])
                                   ),
                        font.label = list(size = 28),
                        vjust = 1.2,
                        # hjust = -1,
                        align = "hv",
                        common.legend = T,
                        legend = "right"
                        
       ),
       device = "png",
       width = 60,
       height = 35,
       units = "cm",
       dpi = 100
)
}

}

{
# error plot----
# for specific time and location, all inputs
time_num <- 20
sp_num <- 2
y_true <- dt_pde_test[[time_num]][,sp_num]
y_pre <- res_pre[[time_num]][,sp_num,]
y_pre_stat <- data.frame(y_true = y_true,
                         med = apply(X = y_pre, MARGIN = 1, FUN = median),
                         lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
                         upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
error_width <- (max(y_true) - min(y_true)) / 30
y_pre_stat %>% ggplot(aes(x = y_true, y = med)) + 
  # geom_linerange(aes(ymin = lower, ymax = upper)) 
  geom_pointrange(aes(ymin = lower, ymax = upper), size =.2)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width) + 
  geom_abline(col = "red")

# for specific time, input, all locations
{
time_num <- 20
sp_num <- c(1:N_sp)
input_num <- 10
y_true <- dt_pde_test[[time_num]][input_num,sp_num]
y_pre <- res_pre[[time_num]][input_num,sp_num,]
y_pre_stat <- data.frame(y_true = y_true,
                         med = apply(X = y_pre, MARGIN = 1, FUN = median),
                         lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
                         upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
y_pre_stat <- y_pre_stat / N_people
error_width <- (max(y_pre_stat["y_true"]) - min(y_pre_stat["y_true"])) / 30
y_pre_stat %>% ggplot(aes(x = y_true, y = med)) + 
  # geom_linerange(aes(ymin = lower, ymax = upper)) 
  geom_pointrange(aes(ymin = lower, ymax = upper), size =.2)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width) + 
  geom_abline(col = "red") + 
  labs(x = "PDE solution", y = "FFBS solution")
}

# for all locations
time_num <- c(3, 5, 10, 20)
for (t in 1:length(time_num)) {
  if(t == 1){
    plot_ls <- list()
  }
  y_pre_stat_full <- data.frame()
  for (i in 1:N_sp) {
    y_true <- dt_pde_test[[time_num[t]]][,i]
    y_pre <- res_pre[[time_num[t]]][,i,]
    y_pre_stat <- data.frame(y_true = y_true,
                             med = apply(X = y_pre, MARGIN = 1, FUN = median),
                             lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
                             upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
    y_pre_stat_full <- rbind(y_pre_stat_full, y_pre_stat)
  }
  
  y_pre_stat_full <- y_pre_stat_full / N_people
  error_width <- (max(y_pre_stat_full["y_true"]) - min(y_pre_stat_full["y_true"])) / 30
  y_true_max <- max(y_pre_stat_full["y_true"])
  
  p <- y_pre_stat_full %>% 
    dplyr::filter(lower > -1, upper < 1) %>% 
    ggplot(aes(x = y_true, y = med)) + 
    geom_smooth() + 
    geom_pointrange(aes(ymin = lower, ymax = upper), size =.2, alpha = 0.5)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width, alpha = 0.5) + 
    geom_abline(col = "red") + 
    # lims(y = c(-y_true_max + 0.1, y_true_max + 0.1)) +
    labs(x = "PDE solution", y = "FFBS solution")
  plot_ls[[t]] <- p
}

ggsave(filename = paste("plot_panel_errorbar4_", as.numeric(Sys.time()), ".png", sep = ""),
       path = path_fig,
       plot = ggarrange(plot_ls[[1]], plot_ls[[2]], plot_ls[[3]],
                        plot_ls[[4]], 
                        ncol = 2, nrow = 2,
                        labels = c(paste("t =", time_num[1]), 
                                   paste("t =", time_num[2]),
                                   paste("t =", time_num[3]),
                                   paste("t =", time_num[4])),
                        font.label = list(size = 28),
                        vjust = 1.2,
                        # hjust = -1,
                        align = "hv",
                        common.legend = T,
                        legend = "right"
                        
       ),
       device = "png",
       width = 42,
       height = 35,
       units = "cm",
       dpi = 100
) 



# tranform results----
# list: input; col: (space)_time_1-(space)_time_nT; row: replication
res_reformback_all <- list()
for (i in 1:n_input_new) {
  res_reformback <- c()
  for (j in 1:nT) {
    res_reformback <- rbind(res_reformback, res_pre[[j]][i,,]) # each row includes replications for a location
  }
  # write.csv(res_reformback, file = paste(path_data, "/res_reformback_", i, ".csv", sep = ""))
  res_reformback_all[[i]] <- res_reformback
}

# choose input not far away
input_full <- rbind(input, input_new)
dist_input_full <- as.matrix(stats::dist(input_full, method = "euclidean", diag = T, upper = T))
dist_to_train <- t(dist_input_full[1:N,(N+1):dim(dist_input_full)[2]])
dist_to_train_mean <- rowMeans(dist_to_train)
order_input <- order(dist_to_train_mean)
sort(dist_to_train_mean)
order_select <- order_input[1:9]

for (i in 1:length(order_select)) {
  if(i == 1){
    plot_ls <- list()
    plot_ls_pt <- list()
  }
  y_true_mat <- dat_pde[[(n_train + order_select[i])]]
  y_true <- matrix(t(y_true_mat), ncol = 1)
  dt_pre <- cal_errorbar_mean(X = res_reformback_all[[order_select[i]]])
  dt_pre <- cbind(y_true, dt_pre) / N_people
  
  error_width <- (max(dt_pre["y_true"]) - min(dt_pre["y_true"])) / 50
  p <- dt_pre %>% 
    dplyr::filter(lower >= -1, upper <= 1) %>% 
    ggplot(aes(x = y_true, y = med)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper), size =.2)+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width) + 
    geom_abline(color = "red") +
    labs(x = "PDE solution", y = "FFBS solution")
    # lims(x = c(0, 1), y = c(-1, 1))
  plot_ls[[i]] <- p
  
  p2 <- dt_pre %>% 
    dplyr::filter(lower >= -1, upper <= 1) %>% 
    ggplot(aes(x = y_true, y = med)) + 
    geom_point(aes(alpha = 0.05)) + 
    geom_smooth(se = T) +
    # geom_pointrange(aes(ymin = lower, ymax = upper), size =.2)+
    # geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width) + 
    geom_abline(color = "red") +
    labs(x = "PDE solution", y = "FFBS solution")
  # lims(x = c(0, 1), y = c(-1, 1))
  plot_ls_pt[[i]] <- p2
}

ggsave(filename = paste("plot_panel_errorbar9_", as.numeric(Sys.time()), ".png", sep = ""),
       path = path_fig,
       plot = ggarrange(plot_ls[[1]], plot_ls[[2]], plot_ls[[3]],
                        plot_ls[[4]], plot_ls[[5]], plot_ls[[6]], 
                        plot_ls[[7]], plot_ls[[8]], plot_ls[[9]],
                        ncol = 3, nrow = 3,
                        # labels = c("(1) Data", "(2) Gaussian Process", "(3) Bilinear", 
                        #            "(4) MBI", "(5) MBA", "(6) Graphical",
                        #            "(7) Inverse Graphical", "(8) Nearest Neighbor"),
                        font.label = list(size = 28),
                        vjust = 1.2,
                        # hjust = -1,
                        align = "hv",
                        common.legend = T,
                        legend = "right"
                        
       ),
       device = "png",
       width = 60,
       height = 50,
       units = "cm",
       dpi = 100
) 

ggsave(filename = paste("plot_panel_pt9_", as.numeric(Sys.time()), ".png", sep = ""),
       path = path_fig,
       plot = ggarrange(plot_ls_pt[[1]], plot_ls_pt[[2]], plot_ls_pt[[3]],
                        plot_ls_pt[[4]], plot_ls_pt[[5]], plot_ls_pt[[6]], 
                        plot_ls_pt[[7]], plot_ls_pt[[8]], plot_ls_pt[[9]],
                        ncol = 3, nrow = 3,
                        # labels = c("(1) Data", "(2) Gaussian Process", "(3) Bilinear", 
                        #            "(4) MBI", "(5) MBA", "(6) Graphical",
                        #            "(7) Inverse Graphical", "(8) Nearest Neighbor"),
                        font.label = list(size = 28),
                        vjust = 1.2,
                        # hjust = -1,
                        align = "hv",
                        common.legend = T,
                        legend = "right"
                        
       ),
       device = "png",
       width = 60,
       height = 50,
       units = "cm",
       dpi = 100
) 
}
