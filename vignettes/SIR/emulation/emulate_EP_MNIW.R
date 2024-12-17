## 1. emulation plots---------------------------
{
# preparation----
# need change: n_input, Nx, Ny, nsam
library(here)
library(lhs)
library(deSolve)
library(ReacTran)

top.dir <- here()
setwd(top.dir)
source(here("R", "init_lib.r"))
seed <- 1234
set.seed(seed)
path <- "emulation"
path_data <- here("data", path)
path_fig <- here("figures", path)

dir.create(path_data, showWarnings=FALSE)
if (!dir.exists(path_fig)) {
  dir.create(path_fig)
}

N_people <- 10000
Gridx <- setup.grid.1D(x.up = 0, x.down = Nx, N = Nx)
Gridy <- setup.grid.1D(x.up = 0, x.down = Ny, N = Ny)

SIR <- function(t, y, parms) {
    S <- matrix(nrow = Nx, ncol = Ny, data = y[1:(Nx*Ny)])
    I <- matrix(nrow = Nx, ncol = Ny, data = y[(Nx*Ny+1) : (2*Nx*Ny)])
    R <- matrix(nrow = Nx, ncol = Ny, data = y[(2*Nx*Ny+1) : (3*Nx*Ny)])
    #NOTE: tran.2D as called here is equivalent to computing the Laplacian of S, times alpha1;
    #	with dC being the result.
    dS <- -eta1*S*I/N_people +
  	  tran.2D(C = S, D.x = alpha1, D.y = alpha1,
  	  dx = Gridx, dy= Gridy,
  	  C.x.up = 0,
  	  C.y.up = 0,
  	  C.x.down = 0,
  	  C.y.down = 0)$dC
    dI <- eta1*S*I/N_people - eta2*I +
          tran.2D(C = I, D.x = alpha2, D.y = alpha2,
                  dx = Gridx, dy = Gridy,
                  C.x.up = 0,
                  C.y.up = 0,
                  C.x.down = 0,
                  C.y.down = 0)$dC
    dR <- eta2*I + 
  	  tran.2D(C = R, D.x = alpha3, D.y = alpha3,
              dx = Gridx, dy = Gridy,
              C.x.up = 0,
              C.y.up = 0,
              C.x.down = 0,
              C.y.down = 0)$dC
    list(c(dS, dI, dR))
}

# Generate pde_para.csv here for all.
sims = 20 #25
cube = maximinLHS(sims,k = 5)
e1 = qunif(cube[,1],2,4)
#set.seed(1)
e2 = qunif(cube[,2],0.2,0.4)
a1 = qunif(cube[,3],0,.2)
a2 = qunif(cube[,4],0,.2)
a3 = qunif(cube[,5],0,.2)

#theta_design = cbind(b, g, a)
theta_design = cbind(e1, e2, a1, a2, a3)

#write out infected
write.csv(theta_design, 
          file=file.path(here(), "data/emulation", paste0("pde_para.csv", sep = "")), 
          row.names=FALSE)


ind_old_data <- TRUE #FALSE
AR_choice <- 2
nsam <- 10
n_input <- 20
times <- seq(0, 50)
nT <- length(times)
nT_ori <- nT
Nx <- 100
Ny <- 100
S <- Nx * Ny

locs <- expand.grid(y=1:Ny, x=1:Nx)
locs <- locs[,c("x", "y")]
locs = as.numeric(rownames(locs))
  
prop_train <- 0.5
gp_tune <- 0.5
gp_sigma2 <- 1.1
gp_tau2 <- 10^(-4)
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

S0 <- matrix(nrow = Nx, ncol = Ny, data = N_people)
I0 <- matrix(nrow = Nx, ncol = Ny, data = 0)
I0[Nx/2,Ny/2] <- 1 # set midpoint to have 50 infected people
I0[round(Nx/10), round(Ny/10)] <- N_people * 0.05
R0 <- matrix(nrow = Nx, ncol = Ny, data = 0)
yini <- c(S0, I0, R0)

# generate actual PDE solutions here.
ysp_old = array(0, dim=c(n_input,length(times), S))
ysp = array(0, dim=c(n_input,S,length(times)))
ysp2 = array(0, dim=c(n_input,length(times),S))
FF1 = array(0, dim=c(n_input, S, length(times)-1))
FFc = array(0, dim=c(n_input, 1, length(times)-1, S))
FF2 = array(0, dim=c(n_input, 1, length(times)-1, S))

for(i in 1:n_input){
  if(i %% 5 == 0){
    print(Sys.time())
    print(paste(i, "/", n_input))
  }
  eta1 = e1[i]
  eta2 = e2[i]
  alpha1 = a1[i]
  alpha2 = a2[i]
  alpha3 = a3[i]
  out <- ode.2D(y = yini, parms = c(eta1=eta1,eta2=eta2,
				    alpha1=alpha1, alpha2=alpha2, alpha3=alpha3), 
                func = SIR,
                nspec = 3, dimens = c(Nx, Ny), times = times,
                lrw = 200000000, names=c("X1", "X2", "X3"))[,-1]

  #write out infected
  write.csv(out[,(1 + Nx * Ny):(2 * Nx * Ny)], 
            file=file.path(path_data, paste0("pde_solution_", i, ".csv", sep = "")), 
	      row.names=FALSE)

  # save the matrix of infected people in one row
  ysp_ = matrix((out[,(Nx*Ny+1) : (2*Nx*Ny)]), ncol=(Nx*Ny))
  ysp_ = log(ysp_+1)
  for(sp in 1:length(locs)){
    ysp_old[i,,sp] = (ysp_[,locs[sp]])
    FF2[i,1,,sp] = (ysp_old[i,-length(times),sp])
    ysp[i,sp,] = (ysp_[,locs[sp]])
    FF1[i,sp,] = (ysp[i,sp,-length(times)])
    ysp2[i,,sp] = (ysp_[,locs[sp]])
    FFc[i,1,,sp] = (ysp2[i,-length(times),sp])
  }
}

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
# R <- gen_gp_kernel(loc = ind_sp_EP, phi = phi_sp, sigma2 = gp_sigma2, tau2 = gp_tau2)
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

}

{
# Plot file ----
## Set up for plotting ----
source(here("R", "func_plot.R"))

## Plot PDE results ----
### Heat Map----
{
  # ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))
  # plot_ls <- list() 
  # tstamp <- as.integer(seq(1, nT, length.out = 9))
  input_num <- 1
  tstamp <- as.integer(seq(1, nT_ori, length.out = 9))
  dat <- dt_pde_test
  max_y <- max(as.vector(unlist(dat))) # set max limit for all plots
  # max_y <- max(as.vector(unlist(res_pre_exact))) # set max limit for all plots
  
  pde_heat <- plot_panel_heatmap_9(dat = dat, tstamp = tstamp, #filename = "plot_panel_pde",
                                   input_num = input_num, max_y = max_y, Nx = Nx, Ny = Ny, nT = nT)
  
  # dat <- postm_ls # this is monte carlo
  dat <- res_pre_exact # this is exact
  ffbs_heat <- plot_panel_heatmap_9(dat = dat, tstamp = tstamp, #filename = "plot_panel_ffbs", savei = F,
                                    input_num = input_num, max_y = max_y, Nx = Nx, Ny = Ny, nT = nT)
  
  # ffbs_heat19 <- plot_panel_heatmap_9(dat = dat, tstamp = 1:9,
  #                                   input_num = input_num, max_y = max_y, Nx = Nx, Ny = Ny, nT = nT)
  
  ggsave(filename = file.path(top.dir,
			      "figures/emulation",
			      paste("plot_heat_compare", ".png", sep = "")),
         path = path_fig,
         plot = ggarrange(pde_heat[[4]], pde_heat[[6]], pde_heat[[8]],
                          ffbs_heat[[4]], ffbs_heat[[6]], ffbs_heat[[8]],
                          ncol = 3, nrow = 2,
                          labels = c(paste("PDE solution: t =", tstamp[4]-1),
                                     paste("t =", tstamp[6]-1),
                                     paste("t =", tstamp[8]-1),
                                     paste("FFBS emulation: t =", tstamp[4]-1),
                                     paste("t =", tstamp[6]-1),
                                     paste("t =", tstamp[8]-1)
                          ),
                          font.label = list(size = 28),
                          vjust = 1.2,
                          hjust = -0.5,
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

{
  ### Error plot----
  #### specific time, location, all inputs (point)----
  {
    alpha_error <- 0.1
    res_pre <- res_pre_MC
    
    time_num <- 15
    sp_num <- 2
    y_true <- dt_pde_test[[time_num]][,sp_num]
    y_pre <- res_pre[[time_num]][,sp_num,]
    y_pre_stat <- data.frame(y_true = y_true,
                             med = apply(X = y_pre, MARGIN = 1, FUN = median),
                             lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
                             upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
    error_width <- (max(y_true) - min(y_true)) / 30 # denominator is for aesthetic 
    y_pre_stat %>% ggplot(aes(x = y_true, y = med)) + 
      # geom_linerange(aes(ymin = lower, ymax = upper)) 
      geom_pointrange(aes(ymin = lower, ymax = upper), size =.2)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width) + 
      geom_abline(col = "red")
  }
  
  #### specific time, input, all locations (field)----
  {
    time_num <- 20
    sp_num <- c(1:N_sp)
    # input_num <- 1
    y_true <- dt_pde_test[[time_num]][input_num,sp_num]
    y_pre <- res_pre[[time_num]][input_num,sp_num,]
    y_pre_stat <- data.frame(y_true = y_true,
                             med = apply(X = y_pre, MARGIN = 1, FUN = median),
                             lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
                             upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
    y_pre_stat <- y_pre_stat / N_people
    error_width <- (max(y_pre_stat["y_true"]) - min(y_pre_stat["y_true"])) / 30
    plot_error_1 <- y_pre_stat %>% ggplot(aes(x = y_true, y = med)) + 
      # geom_linerange(aes(ymin = lower, ymax = upper)) 
      geom_pointrange(aes(ymin = lower, ymax = upper), size =.2, alpha = alpha_error)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width, alpha = alpha_error) + 
      geom_abline(col = "red") + 
      labs(x = "PDE solution", y = "FFBS prediction")
    
    # ggsave(filename = paste("plot_error_1", ".png", sep = ""),
    #        path = path_fig,
    #        plot = plot_error_1,
    #        device = "png",
    #        width = 42,
    #        height = 35,
    #        units = "cm",
    #        dpi = 100
    # ) 
  }
  
  # panel
  time_p <- tstamp[c(4, 6, 8)]
  sp_num <- c(1:N_sp)
  plot_error_comp_ls <- list()
  plot_band_ls <- list()
  y_pre_stat_error_comp_ls <- list()
  for (t in 1:length(time_p)) {
    time_num <- time_p[t]
    y_true <- dt_pde_test[[time_num]][input_num,sp_num]
    y_pre <- res_pre[[time_num]][input_num,sp_num,]
    y_pre_stat <- data.frame(y_true = y_true,
                             med = apply(X = y_pre, MARGIN = 1, FUN = median),
                             lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
                             upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
    y_pre_stat <- y_pre_stat / N_people
    y_pre_stat_error_comp_ls[[t]] <- y_pre_stat
  }
  
  for (t in 1:length(time_p)) {
    # scatter plot
    error_width <- (max(y_pre_stat_error_comp_ls[[t]]["y_true"]) - min(y_pre_stat_error_comp_ls[[t]]["y_true"])) / 30
    plot_error_1 <- y_pre_stat_error_comp_ls[[t]] %>% ggplot(aes(x = y_true, y = med)) + 
      # geom_linerange(aes(ymin = lower, ymax = upper)) 
      geom_pointrange(aes(ymin = lower, ymax = upper), size =.2, alpha = alpha_error / 3)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width, alpha = alpha_error / 3) + 
      geom_abline(col = "red") + 
      labs(x = "PDE solution", y = "FFBS emulation")
    plot_error_comp_ls[[t]] <- plot_error_1
    
    # error band plot
    # fsize <- 34
    plot_band_predict <- y_pre_stat_error_comp_ls[[t]] %>%
      ggplot(aes(x = y_true, y = med)) + 
      geom_ribbon(aes(ymin = lower, ymax = upper, x = y_true), fill = "#A6CEE3", alpha = 1) +
      geom_point(aes(y = med), color = "#1F78B4", size = 0.7, alpha = 1) +
      geom_abline(color = "#E31A1C", linewidth = 1) + 
      labs(x = "PDE solution", y = "FFBS emulation") +
      # theme(
      #   axis.title.x = element_text(size = fsize),
      #   axis.title.y = element_text(size = fsize),
      #   axis.text = element_text(size = fsize)  # Adjust axis tick labels as needed
      # ) + 
      theme_minimal()
    plot_band_ls[[t]] <- plot_band_predict
  }
  
  # ggsave(filename = paste("plot_band_compare", ".png", sep = ""),
  #        path = path_fig,
  #        plot = ggarrange(plot_band_ls[[1]], 
  #                         plot_band_ls[[2]],
  #                         plot_band_ls[[3]],
  #                         ncol = 3, nrow = 1,
  #                         labels = c(paste("t =", tstamp[4]-1),
  #                                    paste("t =", tstamp[6]-1),
  #                                    paste("t =", tstamp[8]-1)),
  #                         font.label = list(size = 28),
  #                         vjust = 1.2,
  #                         hjust = -2,
  #                         align = "hv",
  #                         common.legend = T,
  #                         legend = "right"
  #                         
  #        ),
  #        device = "png",
  #        width = 60,
  #        height = 20,
  #        units = "cm",
  #        dpi = 300
  # )
  
  ggsave(filename = file.path(top.dir,
			      "figures/emulation",
			      paste("plot_error_compare", ".png", sep = "")),
         path = path_fig,
         plot = ggarrange(plot_error_comp_ls[[1]], 
                          plot_error_comp_ls[[2]],
                          plot_error_comp_ls[[3]],
                          ncol = 3, nrow = 1,
                          labels = c(paste("t =", tstamp[4]-1),
                                     paste("t =", tstamp[6]-1),
                                     paste("t =", tstamp[8]-1)),
                          font.label = list(size = 28),
                          vjust = 1.2,
                          hjust = -2,
                          align = "hv",
                          common.legend = T,
                          legend = "right"
                          
         ),
         device = "png",
         width = 60,
         height = 20,
         units = "cm",
         dpi = 300
  )
}

}


## 2. model comparison---------------------------
## START----
{
  # need change: n_input, Nx, Ny, nsam
  library(here)
  setwd(here())
  source(here("R", "init_lib.r"))
  source(here("R", "func_lppd_ana.r"))
  seed <- 1234
  set.seed(seed)
  # path <- "calibration" # Dan_pde_nopermute_train_6x6_26_50_1
  path <- "emulation" # Dan_pde_nopermute_train_100x100_51_20_1
  path_data <- here("data", path)
  path_fig <- here("figures", path)
  if (!dir.exists(path_fig)) {
    dir.create(path_fig)
  }
  
  # read in pde data
  ind_old_data <- F
  AR_choice <- 2
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
  # save(waic_mniw, file = file.path(paste(path_data, "/waic_mniw.RData", sep = "")))
  # save(loglik_mniw, file = file.path(paste(path_data, "/loglik_mniw.RData", sep = "")))
  
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
  # save(waic_mnig, file = file.path(paste(path_data, "/waic_mnig.RData", sep = "")))
  # save(loglik_mnig, file = file.path(paste(path_data, "/loglik_mnig.RData", sep = "")))
  
  
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
  # save(waic_i, file = file.path(paste(path_data, "/waic_i.RData", sep = "")))
  # save(loglik_i, file = file.path(paste(path_data, "/loglik_i.RData", sep = "")))
  
  
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
  
  # print(round(tbl_out, 1))
  # write.csv(tbl_out, file = paste(path_fig, "/tbl_out.csv", sep = ""))
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
  
  ggsave(filename = file.path(top.dir,
			      "figures/emulation",
			      paste("plot_WL_EP.png", sep = "")),
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

# Lotka-Volterra emulation
{
  library(matrixStats)
  library(matrixsampling)
  library(readr)
  library(fst)
  library(lhs)
  library(deSolve)
  library(ReacTran)
  library(ggplot2)
  library(dplyr)
  
  
  top.dir <- here() #gsub("\n", "", read_file("../../pkgdir.txt"))
  #top.dir <- refine.path(top.dir)
  
  source(file.path(top.dir,"R/MNIW.R"))
  source(file.path(top.dir,"R/matrixT.R"))
  
  source(file.path(top.dir,"R/generate_grid.R"))
  source(file.path(top.dir,"R/big_data_set.R"))
  
  source(file.path(top.dir,"R/FF.R"))
  source(file.path(top.dir,"R/BS.R"))
  source(file.path(top.dir,"R/FFBS.R"))
  source(file.path(top.dir,"R/predict.R"))
  source(file.path(top.dir,"R/pds.R"))
  
  # Generate data here.
  dir.create(file.path(top.dir, "data"), showWarnings=FALSE)
  
  p = 1
  # NOTE: Nx is length of x-dimension, or rows, Ny is the length of the y-dimension, or columns.
  
  # No space. Just 2 columns
  nsims_all <- 100
  nsims <- nsims_all/2
  #bnrow <- nsims
  #bncol <- 25
  AR <- 1
  T <- 20 #100 #50
  S <- 2
  
  # NOTE: generate this many latin squares. Jus in case some fail.
  nsims_test <- nsims_all * 10
  
  mu_prey <- 30
  mu_pred <- 4
  u0 <- mu_prey #matrix(nrow = Nx, ncol = Ny, data = rnorm(Nx * Ny, mu_prey, sd_prey))
  v0 <- mu_pred #matrix(nrow = Nx, ncol = Ny, data = rnorm(Nx * Ny, mu_pred, sd_pred))
  
  out.format <- "fst"
  
  #S = length(locs)
  
  zparm_ix <- 3
  logmeans <- c(0, -3, 0, -3)
  logsds <- c(0.5, 0.5, 0.5, 0.5)
  
  # output parameters
  train.out.dir <- file.path(top.dir, sprintf("data/LV_train_Zparms%d_u%dv%d_AR%d_T%d", zparm_ix, u0, v0, AR, T))
  dir.create(train.out.dir, showWarnings=FALSE)
  train.out.head <- "LV_train"
  
  source(file.path(top.dir, "R/generate_grid.R"))
  
  parameter_fname_train <- file.path(train.out.dir, "parameters_train.csv")
  
  test.out.dir <- file.path(top.dir, sprintf("data/LV_test_Zparms%d_u%dv%d_AR%d_T%d", zparm_ix, u0, v0, AR, T))
  dir.create(test.out.dir, showWarnings=FALSE)
  
  parameter_fname_test <- file.path(test.out.dir, "parameters_test.csv")
     
  e1 <- e2 <- e3 <- e4 <- NULL
  
  set.seed(12345)
      
  cube = maximinLHS(nsims_test, k = 4)
  # write this cube
  cube_df <- as.data.frame(cube)
  names(cube_df) <- c("eta1", "eta2", "eta3", "eta4")

  # write the file later.
  print("Initializing parameters.")

  # Some constants yield non-solutions
  e1 = exp(qnorm(cube[,1], mean=logmeans[1], sd=logsds[1]))
  e2 = exp(qnorm(cube[,2], mean=logmeans[2], sd=logsds[2]))
  e3 = exp(qnorm(cube[,3], mean=logmeans[3], sd=logsds[3]))
  e4 = exp(qnorm(cube[,4], mean=logmeans[4], sd=logsds[4]))

  LV <- function(t, y, parms) {
      u <- matrix(nrow = 1, ncol = 1, data = y[1])
      v <- matrix(nrow = 1, ncol = 1, data = y[2])
  
      du <- eta1 * u - eta2 * u * v
      dv <- -eta3 * v + eta4 * u * v
  
      list(c(du, dv))
  }

  yini <- c(u0, v0)

  times <- seq(0, T)

  # Generate LV data from the LV ODE.
  # If last Y file and last F file do not exist, generate the data.
  
  eta_out <- data.frame(eta1 = rep(0, nsims),
                        eta2 = rep(0, nsims),
                        eta3 = rep(0, nsims),
                        eta4 = rep(0, nsims))
  
  ysp = array(0, dim=c(nsims,S,length(times)))
  
  ix_y <- 1
  ix_f <- 0
  
  for(i in 1:nsims_test){
      eta1 = e1[i]
      eta2 = e2[i]
      eta3 = e3[i]
      eta4 = e4[i]
  
      parms <- c(eta1, eta2, eta3, eta4)
  	
      out <- tryCatch({ode(y = yini, parms = parms,
                    func = LV, times = times#,
  		      )[,-1]
                       }, error = function(e) {
                           print(e)
  		         return(NULL)
                       })
  
      if (is.null(out)) {
          print("This set of parameters could not produce data. Proceeding to next parameters.")
          next
      }
    
  	# save both predator and prey.
      ysp_ = out #matrix((out[,(Nx*Ny+1) : (2*Nx*Ny)]), ncol=(Nx*Ny))
      ysp_ = log(ysp_)
  
      # skip if bad data
      if (any(is.nan(ysp_))) { #stop("We have NaN's for a certain set of parameters!")
          print("This set of parametes produced bad data with NaN's. Proceeding to next parameters.")
          next
      }
  
      ysp[ix_y,1,] <- ysp_[,1]
      ysp[ix_y,2,] <- ysp_[,2]
  
      eta_out[ix_y,"eta1"] <- eta1
      eta_out[ix_y,"eta2"] <- eta2
      eta_out[ix_y,"eta3"] <- eta3
      eta_out[ix_y,"eta4"] <- eta4
  
      ix_y <- ix_y + 1
      if (ix_y > nsims) {
          print("Aquired all required training samples for each parameter. Continuing with data production.")
          ix_f <- i
          break
      }
  }
  
  if (ix_y <= nsims) stop("Did not produce enough good data. Please adjust parameters or take a larger sample.")
    
  # allocate new F_t as the row-wise concatenation of the past two rounds.
  # reshape ysp into space x simulation matrices, separated by time.
  F_t <- matrix(nrow = nsims, ncol = AR * S, data = 0)
 
  # reshape the data into 100x100 spatial blocks per line and stack the blocks per sim 
  for (t in times) {
      y_t <- ysp[,,t+1]

#        y_t <- permute_data(y_t, Nx, Ny, sqrt(bncol), sqrt(bncol))
      if (out.format == "csv") {
          write.csv(y_t, file=file.path(train.out.dir, sprintf("%s_Y%d.csv", train.out.head, t)),
                    row.names=FALSE, col.names=TRUE)

	} else if (out.format == "fst") {
          write.fst(as.data.frame(y_t), file.path(train.out.dir, sprintf("%s_Y%d.fst", train.out.head, t)))
	}

      if (t < AR) next

	# alternate the resulting matrices
#        for (pp in 1:p) F_t[,seq(pp, p*S, p)] <- ysp[,,t+1-pp]
      for (pp in 1:AR) F_t[,seq((pp-1)*S + 1, pp*S)] <- ysp[,,t+1-pp]
#        for (pp in 1:p) F_t[,((pp-1)*S+1):(pp*S)] <- ysp[,,t+1-pp]

      if (out.format == "csv") {
          write.csv(F_t, file=file.path(train.out.dir, sprintf("%s_F%d.csv", train.out.head, t)),
                    row.names=FALSE, col.names=TRUE)

      } else if (out.format == "fst") {
          write.fst(as.data.frame(F_t), file.path(train.out.dir, sprintf("%s_F%d.fst", train.out.head, t)))
      }
  }

  # Generate test data.
  # include the prediction data generation here.
  nsims_new <- nsims
  ix_y <- 1

  # generate and predict new data
  #print("Generating new data with existing parameters.")
  print("Generating new data with new parameters.")
  test.out.head <- "LV_test"
   
  if (file.exists(parameter_fname_test)) {
      param_df_test <- read.csv(parameter_fname_test)
  
      e1new <- param_df_test$eta1[1:nsims_new]
      e2new <- param_df_test$eta2[1:nsims_new]
      e3new <- param_df_test$eta3[1:nsims_new]
      e4new <- param_df_test$eta4[1:nsims_new]
  } 
  
  #gen_SIR_data(nsims_new, e1new, e2new,
  #	               a1new, a2new, a3new, 
  #		       Nx, Ny, times, locs, out.head = test.out.head)
  
  ysp = array(0, dim=c(nsims_new,S,length(times)))

  eta_out_new <- data.frame(eta1 = rep(0, nsims_new),
                            eta2 = rep(0, nsims_new),
                            eta3 = rep(0, nsims_new),
                            eta4 = rep(0, nsims_new))

  print("Initializing new data...")
 
#    if (!file.exists(file.path(test.out.dir, sprintf("%s_F%d.csv", test.out.head, max(times))))) {
  for(i in 1:(nsims_test - ix_f)){
      eta1 = e1[ix_f + i]
      eta2 = e2[ix_f + i]
      eta3 = e3[ix_f + i]
      eta4 = e4[ix_f + i]

      parms <- c(eta1=eta1, eta2=eta2, eta3=eta3, eta4=eta4)

	#TODO: try/catch this
      out <- tryCatch({ode(y = yini, parms = parms,
                    func = LV, times = times#,
#                      nspec = 2, dimens = NULL, 
#                      lrw = 200000000,
#		      names=c("X1", "X2")
		          )[,-1]
                       }, error = function(e) {
                           print(e)
		             return(NULL)
                       })
      if (is.null(out)) {
          print("This set of parameters could not produce data. Proceeding to next parameters.")
          next
      }
    
    #  if (i == sims) write.table(out, file=file.path(test.out.dir, "SIR_solution.csv"), 
    #	      sep=",", row.names=FALSE, col.names=FALSE)
   

      # save the matrix of infected people in one row
#            ysp_ = matrix((out[,(Nx*Ny+1) : (2*Nx*Ny)]), ncol=(Nx*Ny))
      ysp_ = out
      ysp_ = log(ysp_)
	    
      if (any(is.nan(ysp_))) { #stop("We have NaN's for a certain set of parameters!")
          print("This set of parametes produced bad data with NaN's. Proceeding to next parameters.")
          next
      }

      ysp[ix_y, 1,] <- ysp_[,1]
      ysp[ix_y, 2,] <- ysp_[,2]

      eta_out_new[ix_y, "eta1"] <- eta1
      eta_out_new[ix_y, "eta2"] <- eta2
      eta_out_new[ix_y, "eta3"] <- eta3
      eta_out_new[ix_y, "eta4"] <- eta4

      ix_y <- ix_y + 1
      if (ix_y > nsims_new) {
          print("Aquired all required test samples for each parameter. Continuing with data production.")
          break
      }
  }
  
  # reshape ysp into space x simulation matrices, separated by time.
 
  # allocate new F_t as the row-wise concatenation of the past two rounds.
  # reshape ysp into space x simulation matrices, separated by time.
#    F_t <- matrix(nrow = nsims, ncol = p * S, data = 0)
  F_t <- matrix(nrow = nsims_new, ncol = AR * S, data = 0)
 
  # reshape the data into 100x100 spatial blocks per line and stack the blocks per sim 
  for (t in times) {
      y_t <- ysp[,,t+1]
 
      if (out.format == "csv") {
          write.csv(y_t, file=file.path(test.out.dir, sprintf("%s_Y%d.csv", test.out.head, t)),
                    row.names=FALSE, col.names=TRUE)

	    } else if (out.format == "fst") {
          write.fst(as.data.frame(y_t), file.path(test.out.dir, sprintf("%s_Y%d.fst", test.out.head, t)))
      }	    
      
      if (t < AR) next

      for (pp in 1:AR) F_t[,seq((pp-1)*S + 1, pp*S)] <- ysp[,,t+1-pp]

      if (out.format == "csv") {
          write.csv(F_t, file=file.path(test.out.dir, sprintf("%s_F%d.csv", test.out.head, t)),
                    row.names=FALSE, col.names=TRUE)

      } else if (out.format == "fst") {
          write.fst(as.data.frame(F_t), file.path(test.out.dir, sprintf("%s_F%d.fst", test.out.head, t)))
      }	    
  }
    

  if (ix_y <= nsims_new) stop("Did not produce enough good test data. Please adjust parameters or take a larger sample.")

  write.csv(eta_out, file=parameter_fname_train, row.names=FALSE)
  write.csv(eta_out_new, file=parameter_fname_test, row.names=FALSE)

  print("New data initialized.")

  # setup for running the FFBS on Lotka-Volterra
  test_new_predictions <- TRUE
  
  # start time to time the run of the program from start to finish.
  t0 <- Sys.time()
  
  times <- seq(0, T)#,by=10)
  
  print(sprintf("Running results with AR %d, T = %d, u0 = %d, v0 = %d...", AR, T, u0, v0))
  # output parameters
  train.data.dir <- file.path(top.dir, sprintf("data/LV_train_Zparms%d_u%dv%d_AR%d_T%d", zparm_ix, u0, v0, AR, T))
  test.data.dir <- file.path(top.dir, sprintf("data/LV_test_Zparms%d_u%dv%d_AR%d_T%d", zparm_ix, u0, v0, AR, T))
  
  train.out.head <- "LV_train"
  
  parameter_fname <- file.path(train.data.dir, "parameters_train.csv")
     
  param_df <- read.csv(parameter_fname)
  e1 <- param_df$eta1
  e2 <- param_df$eta2
  e3 <- param_df$eta3
  e4 <- param_df$eta4
  
  nsims=nrow(param_df)
  
  # generate Vt_mat
  X2_all <- matrix(nrow=nsims, ncol=nsims, data=0)
  for (param in list(e1=e1, e2=e2, e3=e3, e4=e4)) {
      Xp <- matrix(rep(param, nsims), nrow=nsims)
      X2_all <- X2_all + (Xp - t(Xp))^2
  }
  
  X2_all <- sqrt(X2_all)
  dmax <- max(X2_all)
  deff <- 0.5 * dmax
  #beta <- 0.4 * 3/dmax
  beta <- 3/deff
  tau2 <- 1e-4 #0.04
  Vt_mat <- exp(-beta * X2_all) + diag(nsims) * tau2
  
  # check the parameter estimates
  #print(sprintf("eta = (%.3f, %.3f), alpha = (%.3f, %.3f, %.3f)", e1[sims], e2[sims], a1[sims], a2[sims], a3[sims]))
  
  theta_design = cbind(e1, e2, e3, e4)
  
  # package into a big_data_set object
  Y_filename_list <- file.path(train.data.dir, sprintf("%s_Y%d.%s", train.out.head, seq(AR,max(times)), out.format))
  F_filename_list <- file.path(train.data.dir, sprintf("%s_F%d.%s", train.out.head, seq(AR,max(times)), out.format))
  
  bnrow <- r <- Inf
  bncol <- S
  ptilde <- AR * bncol #p * min(c, S)
  
  #bncol <- c <- Ny #if (Nx == 100) 100 else if (Nx = 10) 25
  nrowblock_list <- rep(bnrow, length(times))
  ncolblock_list <- rep(bncol, length(times))
  
  train_data_set <- big_data_set(Y_filename_list,
                                 F_filename_list,
                                 rep(nsims, length(times)), 
                                 rep(bncol, length(times)), 
                                 rep(AR * bncol, length(times)), 
                                     nrowblock_list = nrowblock_list,
                                 ncolblock_list = ncolblock_list,
                                 F_ncolblock_list = AR * ncolblock_list,
                                 header = TRUE,
                                 split_col_F = TRUE,
                                 grid.traversal.mode="rowsnake",
                                 uniform_grid = TRUE)
  
  # run the FF.bigdata here.
  n0 <- AR * bncol + 3
  D0 <- diag(bncol)
  d0 <- AR * bncol + bncol
  m0 <- matrix(rep(1, ptilde*bncol), nrow=ptilde)
  M0 <- diag(ptilde)
  Gtmat <- matrix(c(diag(ptilde)), nrow=1)
  Wtmat <- matrix(c(diag(ptilde)), nrow=1)
  
  # The Vt is a row-covariance matrix. Sigma is the column covariance. Vt can't be spatial in this case.
  
  # Correlation between rows. Use that square distance instead. How to compute? c x c = 25 x 25 cov matrix?
  # The spatial Vt takes in a data chunk and infers the coordinate.
  Vt <- function(chunk) {
      return(Vt_mat)
  }
  
  # generates the spatial correlation matrix of size r.
  parameter_fname_new <- file.path(test.data.dir, "parameters_test.csv")
  
  # predict new data
  #print("Generating new data with existing parameters.")
  print("Predicting new data with new parameters.")
  
  param_df_new <- read.csv(parameter_fname_new)
  nsims_new <- nrow(param_df_new)
  
  e1new <- param_df_new$eta1[1:nsims_new]
  e2new <- param_df_new$eta2[1:nsims_new]
  e3new <- param_df_new$eta3[1:nsims_new]
  e4new <- param_df_new$eta4[1:nsims_new]
  
  test.out.head <- "LV_test"
  
  Y_filename_list <- file.path(test.data.dir, sprintf("%s_Y%d.%s", test.out.head, seq(AR,max(times)), out.format))
  F_filename_list <- file.path(test.data.dir, sprintf("%s_F%d.%s", test.out.head, seq(AR,max(times)), out.format))
  
  test_data_set <- big_data_set(Y_filename_list,
                                F_filename_list, 
  			              rep(nsims_new, length(times)), 
  			              rep(bncol, length(times)), 
  			              rep(AR * bncol, length(times)), 
                                nrowblock_list = nrowblock_list,
                                ncolblock_list = ncolblock_list,
                                F_ncolblock_list = AR * ncolblock_list,
                                header = TRUE,
                                split_col_F = TRUE,
                                grid.traversal.mode="rowsnake")
  
  # generate Jt
  params_old <- list(e1 = e1, e2 = e2, e3 = e3, e4 = e4)
  params_new <- list(e1 = e1new, e2 = e2new, e3 = e3new, e4 = e4new)
  
  r_new <- nsims_new
  r_old <- nsims    
  
  XXtilde_all <- matrix(ncol=r_new, nrow = r_old, data=0)
  
  for (param in names(params_old)) {
      X_c <- matrix(rep(params_old[[param]], r_new), ncol=r_new)
      Xtilde_c <- matrix(rep(params_new[[param]], r_old), ncol=r_old)
  
      XXtilde_all <- XXtilde_all + (X_c - t(Xtilde_c))^2
  }
  XXtilde_all <- sqrt(XXtilde_all)
  
  Jtmat <- exp(-beta * XXtilde_all)
  
  Jt <- function(new_chunk, old_chunk) {
      return(Jtmat)
  }
  
  # produce Vt_new
  param_df_new <- read.csv(parameter_fname_new)
  nsims_new <- nrow(param_df_new)
  
  e1new <- param_df_new$eta1
  e2new <- param_df_new$eta2
  e3new <- param_df_new$eta3
  e4new <- param_df_new$eta4
  
  # generate Vt_mat
  X2_all_new <- matrix(nrow=nsims_new, ncol=nsims_new, data=0)
  
  for (param in list(e1=e1new, e2=e2new, e3=e3new, e4=e4new)) {
      Xp <- matrix(rep(param, nsims_new), nrow=nsims_new)
      X2_all_new <- X2_all_new + (Xp - t(Xp))^2
  }
  
  X2_all_new <- sqrt(X2_all_new)
  dmax_new <- max(X2_all_new)
  #beta <- 0.4 * 3/dmax_new
  tau2 <- 1e-4 #0.04
  Vt_mat_new <- exp(-beta * X2_all_new) + diag(nsims_new) * tau2
  Vt_new <- function(chunk) {
      return(Vt_mat_new)
  }
  
  Sigma_mode <- "IW"
  print(sprintf("Running Sigma mode %s...", Sigma_mode))
  # create different out directories based on the id, IW, and IG Sigma structures.
  ffbs.out.dir <- file.path(top.dir, sprintf("vignettes/emulation/LV_train_Zparms%d_u%dv%d_AR%d_T%d_%s", 
  					       zparm_ix, u0, v0, AR, T, Sigma_mode))
  dir.create(ffbs.out.dir, showWarnings=FALSE)
  
  # run the FFBS here
  ff.out.head <- file.path(ffbs.out.dir, "FF_out")
  smooth.out.head <- file.path(ffbs.out.dir, "BS_out")
  
  t0 <- Sys.time()
  
  FF_args <- list(data_set = train_data_set, m0=m0, M0=M0,
  			            n0 = n0,
  		            Gtmat=Gtmat, Wtmat=Wtmat, Vt = Vt,
                              out.head = ff.out.head,
  	                            out.format = out.format,
  	                            fst.head = "V")
  
  #        FF_args$n0 <- n0
  if (Sigma_mode == "IW") {
      FF_args$D0 <- D0
      FF_args$d0 <- NULL
      FF_args$R <- NULL
  } else if (Sigma_mode == "IG") {
      FF_args$D0 <- NULL
      FF_args$d0 <- d0
      FF_args$R <- R
  } else if (Sigma_mode == "id") {
      FF_args$D0 <- NULL
      FF_args$d0 <- d0
      FF_args$R <- diag(bncol)
  }
  
  # Run the Forward Filter
  ff.file.paths <- do.call(FF.bigdata, FF_args)
  
  # include the BS.bigdata
  smooth.file.paths <- smooth.bigdata(ff.file.paths$mt_files, ff.file.paths$Mt_files, Gtmat, Wtmat,
  					    p = AR * bncol, c = bncol,
                                      out.head = smooth.out.head, out.format = out.format)
  
  print(sprintf("Time to run the FFBS: %.3f", Sys.time() - t0))
  
  # predicting new data here.
  if (test_new_predictions) {
      t0 <- Sys.time()
      pred.out.dir <- file.path(top.dir,
  			  sprintf("vignettes/emulation/LV_predict_Zparms%d_u%dv%d_AR%d_T%d_%s", 
  				   zparm_ix, u0, v0, AR, T, Sigma_mode))
      dir.create(pred.out.dir, showWarnings=FALSE)
      print("Predicting new data...")
  
      pred_head <- file.path(pred.out.dir,"Ytestpred")
      predict(smooth.file.paths, test_data_set, train_data_set, 
              Vt_new = Vt_new, Vt_old = Vt, Jt = Jt,
              out.head = pred_head,
              out.format = out.format,
              fst.head = "V",
              verbose=TRUE) 
      
      print(sprintf("New data predicted. Elapsed time: %.3f", Sys.time() - t0))
  }

  #####
  # Generate correlation plots.
  #####

  source(file.path(top.dir, "R/pds.R"))
  
  emu.dir <- file.path(top.dir, 
		       "vignettes/emulation",
		       sprintf("LV_train_Zparms3_u30v4_AR%d_T%d_IW", AR, T))
  #emu.dir <- sprintf("LV_train_Zparms2_u30v4_AR%d_T%d_IW", AR, T)
  
  #TODO: Read in nT and DT 
  nT <- read.csv(file.path(emu.dir, sprintf("FF_out_nt_%d.csv", T - AR + 1)), header=FALSE)
  nT <- unlist(nT[length(nT)])
  
  DT <- read.fst(file.path(emu.dir, sprintf("FF_out_Dt_uptri_%d.fst", T - AR + 1)))
  DT <- tri.to.sym(unlist(DT[2,]), S, lower=FALSE)
  
  L <- 20000
  
  set.seed(98765)
  
  Sigmas <- rinvwishart(L, nT, DT)
  
  corrs <- Sigmas[1,2,]/sqrt(Sigmas[1,1,] * Sigmas[2,2,])
  
  Sigma_11_quantiles <- quantile(Sigmas[1,1,], probs = c(0.025, 0.5, 0.975))
  Sigma_12_quantiles <- quantile(Sigmas[1,2,], probs = c(0.025, 0.5, 0.975))
  Sigma_22_quantiles <- quantile(Sigmas[2,2,], probs = c(0.025, 0.5, 0.975))
  corrs_quantiles <- quantile(corrs, probs = c(0.025, 0.5, 0.975))
  
  c_p <- ggplot(data=data.frame(correlations = corrs), aes(x=correlations)) + geom_histogram(fill="lightgreen", color="black") +
  	theme_bw() + geom_vline(xintercept=corrs_quantiles[2], color="black", linetype="solid") +
  	geom_vline(xintercept=corrs_quantiles[1], color="black", linetype="dashed") +
  	geom_vline(xintercept=corrs_quantiles[3], color="black", linetype="dashed")
  
  ggsave(c_p, file=file.path(top.dir, 
		             "figures/emulation",
			     "LV_corr_hist.png"))
  
  S11_p <- ggplot(data=data.frame(Sigma_11 = Sigmas[1,1,]), aes(x=Sigma_11)) + geom_histogram(fill="orange", color="black") +
  	theme_bw() + geom_vline(xintercept=Sigma_11_quantiles[2], color="black", linetype="solid") +
  	geom_vline(xintercept=Sigma_11_quantiles[1], color="black", linetype="dashed") +
  	geom_vline(xintercept=Sigma_11_quantiles[3], color="black", linetype="dashed")
  
  ggsave(S11_p, file=file.path(top.dir, 
		               "figures/emulation",
			       "LV_Sigma11_hist.png"))
  
  S12_p <- ggplot(data=data.frame(Sigma_12 = Sigmas[1,2,]), aes(x=Sigma_12)) + geom_histogram(fill="lightgreen", color="black") +
  	theme_bw() + geom_vline(xintercept=Sigma_12_quantiles[2], color="black", linetype="solid") +
  	geom_vline(xintercept=Sigma_12_quantiles[1], color="black", linetype="dashed") +
  	geom_vline(xintercept=Sigma_12_quantiles[3], color="black", linetype="dashed")
  
  ggsave(S12_p, file=file.path(top.dir, 
		               "figures/emulation",
			       "LV_Sigma12_hist.png"))
  
  S22_p <- ggplot(data=data.frame(Sigma_22 = Sigmas[2,2,]), aes(x=Sigma_22)) + geom_histogram(fill="lightyellow", color="black") +
  	theme_bw() + geom_vline(xintercept=Sigma_22_quantiles[2], color="black", linetype="solid") +
  	geom_vline(xintercept=Sigma_22_quantiles[1], color="black", linetype="dashed") +
  	geom_vline(xintercept=Sigma_22_quantiles[3], color="black", linetype="dashed")
  
  ggsave(S22_p, file=file.path(top.dir,
		               "figures/emulation",
			       "LV_Sigma22_hist.png"))
  
  #####
  # generate prediction plots below
  #####

  zparm_ix <- 3
  
  #data.dir <- sprintf("../../../data/LV_test_Zparms_u%dv%d_AR%d_T%d", u0, v0, AR, T)
  #emu.dir <- sprintf("../emulation/LV_train_Zparms_u%dv%d_AR%d_T%d_%s", u0, v0, AR, T, Sigma_mode)
  #pred.dir <- sprintf("LV_predict_Zparms_u%dv%d_AR%d_T%d_%s", u0, v0, AR, T, Sigma_mode)
  
  Sigma_mode <- "IW"
  data.dir <- file.path(top.dir, 
			sprintf("data/LV_test_Zparms%d_u%dv%d_AR%d_T%d", 
				zparm_ix, u0, v0, AR, T))
  emu.dir <- file.path(top.dir, 
		       "vignettes/emulation",
		       sprintf("LV_train_Zparms%d_u%dv%d_AR%d_T%d_%s", 
			       zparm_ix, u0, v0, AR, T, Sigma_mode))
  pred.dir <- file.path(top.dir,
		       "vignettes/emulation",
			sprintf("LV_predict_Zparms%d_u%dv%d_AR%d_T%d_%s", 
				zparm_ix, u0, v0, AR, T, Sigma_mode))
  
  params_test <- read.csv(file.path(data.dir, "parameters_test.csv"))
  N <- nrow(params_test)
      
  Yall.file.name <- file.path(pred.dir, sprintf("Y_all_Zparms%d.fst", zparm_ix))
  
  if (!file.exists(Yall.file.name)) {
      # join the dataframes across time into one.
      Y_all <- data.frame("t"=numeric(),
      		    "eta_ix"=numeric(), 
      		    "V1"=numeric(),
      		    "V2"=numeric())
      
      for (t in 0:T) {
          # read in the file 
          Yt <- read.fst(file.path(data.dir, sprintf("LV_test_Y%d.fst", t)))
      
          # Assign eta and t
          Yt$t <- t
          Yt$eta_ix <- seq(1, N)
      
          # append to data_all
          Y_all <- rbind(Y_all, Yt)
      }
  
      # rename to something better.
      colnames(Y_all)[colnames(Y_all) == "V1"] <- "ut"
      colnames(Y_all)[colnames(Y_all) == "V2"] <- "vt"
  
      # save it somewhere. If it already exists. Read it back in.
      write.fst(Y_all, Yall.file.name)
  } else {
      Y_all <- read.fst(Yall.file.name)
  }

  # plot the predicted Y_t for one eta alongside the ground truth for that same eta. Use the original data. Outer join???
  Y_pred_all <- data.frame("t"=numeric(),
      		    "eta_ix"=numeric(), 
      		    "utpred"=numeric(),
      		    "vtpred"=numeric(),
  		    "rowcov"=numeric(),
                      "utpredsd"=numeric(),
                      "vtpredsd"=numeric())
  
  uptri_to_diag <- function(M_uptri, p) {
      ix_arr <- rep(0, p)
      ix_next <- 1
      for (i in 1:p) {
          ix_arr[i] <- ix_next
          ix_next <- ix_next + i + 1
      }
      return(M_uptri[ix_arr])
  }
  
  # get D_T and n_T column covariance matrices.
  D_T <- read.fst(file.path(emu.dir, sprintf("FF_out_Dt_uptri_%d.fst", T - AR + 1)))
  D_T <- unlist(D_T[nrow(D_T),])
  
  n_T <- read.csv(file.path(emu.dir, sprintf("FF_out_nt_%d.csv", T - AR + 1)), header=FALSE)
  n_T <- unlist(n_T)[length(n_T)]
  
  print(sprintf("Correlation between log(u_t) and log(v_t): %.7f", D_T[2]/sqrt(D_T[1] * D_T[3])))
  
  for (t in AR:T) {
      # read in the file 
      Ytpred <- as.data.frame(matrix(unlist(read.fst(file.path(pred.dir, sprintf("Ytestpred%d.fst", t - AR + 1)))), nrow=N, ncol=2))
      names(Ytpred) <- c("utpred", "vtpred")
      Ytpred <- cbind(data.frame(t=t, eta_ix=seq(1,N)), Ytpred)
  
      #get row covariance matrix.
      #Note that it's the upper-triangle. We need just the diagonals. 
      Ytpredcovdiag <- uptri_to_diag(read.fst(file.path(pred.dir,
  						      sprintf("Ytestpred_var_%d.fst", t - AR + 1))), 
  				   N)
  
      # individual variance is (the diagonal element) x (the diagonal element of Sigma_T, depending on whether we consider ut or vt).
      Ytpred$rowcov <- unlist(Ytpredcovdiag[1,])
  
      Ytpred$utpredsd <- sqrt(Ytpred$rowcov * D_T[1] / (n_T - 3))
      Ytpred$vtpredsd <- sqrt(Ytpred$rowcov * D_T[3] / (n_T - 3))
  
      Y_pred_all <- rbind(Y_pred_all, Ytpred)
  }
  
  Y_all <- dplyr::left_join(Y_all, Y_pred_all)#, on=c("eta_ix", "t"))
  
  
  # plot for individual eta_ixs
  etas_output <- c(6,37,41)  
  for (eta_ix_plot in etas_output) {
  
      p_Yt <- ggplot(data=subset(Y_all, eta_ix == eta_ix_plot)) +#eta_ix)) +#, group=eta_ix)) + 
                     geom_point(aes(x=utpred, y=vtpred, color=rgb(1,0,0,1))) + #geom_path(aes(group=eta))
                     geom_errorbar(aes(x=utpred,
      			             y=vtpred,
      			  ymin=vtpred - qnorm(0.975) * vtpredsd,
      	                  ymax=vtpred + qnorm(0.975) * vtpredsd,
      			  width=.2,
      			  color=rgb(1,0,0,1))) +
                     geom_errorbarh(aes(#x=utpred,
      			   y=vtpred,
      			   xmin=utpred - qnorm(0.975) * utpredsd,
      	                   xmax=utpred + qnorm(0.975) * utpredsd,
      			   height=.2,
      			   color=rgb(1,0,0,1))) +
                     geom_point(aes(x=ut, y=vt, color=rgb(0,1,0,1))) +
              theme_bw() +
              theme(legend.position="none") +
  	    xlab("log(prey_pop.)") +
  	    ylab("log(predator_pop.)")
      
      ggsave(file.path(top.dir,
		       "figures/emulation",
		       sprintf("LV_plots_pred_%d.png", eta_ix_plot)), plot=p_Yt)
  
  }
  
  Y_all <- data.frame("t"=numeric(),
  		    "eta_ix"=numeric(), 
  		    "V1"=numeric(),
  		    "V2"=numeric())
  
  for (t in 0:T) {
      # read in the file 
      Yt <- read.fst(file.path(data.dir, sprintf("LV_test_Y%d.fst", t)))
  
      # Assign eta and t
      Yt$t <- t
      Yt$eta_ix <- seq(1, N)
  
      # append to data_all
      Y_all <- rbind(Y_all, Yt)
  }
  
  font.size <- 8
  
  # rename to something better.
  colnames(Y_all)[colnames(Y_all) == "V1"] <- "ut"
  colnames(Y_all)[colnames(Y_all) == "V2"] <- "vt"
  
  #Y_all <- dplyr::left_join(Y_all, Y_pred_all)#, on=c("eta_ix", "t"))
  Y_all <- merge(Y_all, Y_pred_all)#, on=c("eta_ix", "t"))
  
  Y_all$ut_covered <- (Y_all$ut > Y_all$utpred - qnorm(0.975) * Y_all$utpredsd) & 
  	            (Y_all$ut < Y_all$utpred + qnorm(0.975) * Y_all$utpredsd)
  
  Y_all$vt_covered <- (Y_all$vt > Y_all$vtpred - qnorm(0.975) * Y_all$vtpredsd) & 
  	            (Y_all$vt < Y_all$vtpred + qnorm(0.975) * Y_all$vtpredsd)
  
  # plot them separately, use the actual data for Y_all here.
  for (eta_ix_plot in etas_output) {
      # plot the "QQ"-plot. Report coverage rates in lower-right corner.
      ut_coverrate <- sum(subset(Y_all, 
  			       (eta_ix == eta_ix_plot) & 
  				       (t>= AR) & (t <= T))$ut_covered)/sum((Y_all$eta_ix == eta_ix_plot) & 
                                                                              (Y_all$t >= AR) * (Y_all$t <= T))
      vt_coverrate <- sum(subset(Y_all, 
  			       (eta_ix == eta_ix_plot) & 
  				       (t>= AR) & (t <= T))$vt_covered)/sum((Y_all$eta_ix == eta_ix_plot) & 
                                                                              (Y_all$t >= AR) * (Y_all$t <= T))
  
      x_coord <- if (eta_ix_plot == 37) 2 else 3.5
  
      p_qq_ut <- ggplot(data=subset(Y_all, (eta_ix == eta_ix_plot) & (t>= AR) & (t <= T)),
  		      aes(x=ut, y=utpred)) +
  	    geom_point(color="darkgreen") +
  	    geom_errorbar(aes(ymin = utpred - qnorm(0.975) * utpredsd,
  			      ymax = utpred + qnorm(0.975) * utpredsd),
  			  color="darkgreen", alpha=0.7) +
              geom_abline() +
  	    theme_bw() + theme(legend.position="none") +
  	    xlab("log(Population)") + ylab("est. log(Population)") +
  	    geom_label(aes(x=x_coord,y=-3,
                         label=sprintf("Coverage Rate: %.2f", ut_coverrate)
                         ), size=font.size)
  
      p_qq_vt <- ggplot(data=subset(Y_all, (eta_ix == eta_ix_plot) & (t>= AR) & (t <= T)),
  		      aes(x=vt, y=vtpred)) +
  	    geom_point(color="red") +
  	    geom_errorbar(aes(ymin = vtpred - qnorm(0.975) * vtpredsd,
  			      ymax = vtpred + qnorm(0.975) * vtpredsd),
  			  color="red", alpha=0.7) +
              geom_abline() +
  	    theme_bw() + theme(legend.position="none") +
  	    xlab("log(Population)") + ylab("est. log(Population)") +
  	    geom_label(aes(x=3,y=-1,
                         label=sprintf("Coverage Rate: %.2f", vt_coverrate)
                         ), size=font.size)
      
      ggsave(file.path(top.dir, 
		       "figures/emulation",
		       sprintf("LV_plots_ut_qq_%d.png", eta_ix_plot)), plot=p_qq_ut)
      ggsave(file.path(top.dir, 
		       "figures/emulation",
		       sprintf("LV_plots_vt_qq_%d.png", eta_ix_plot)), plot=p_qq_vt)
  }
}

