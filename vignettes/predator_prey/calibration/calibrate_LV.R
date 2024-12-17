
{
  # plot separate lynx-hare point data and the artifically-generated curves.
  library(here)
  library(fst)
  library(ggplot2)

  top.dir <- here()
  
  zreal <- read.csv(file.path(top.dir,
			      "data/lynx_hare_df.csv"))
  
  u0 <- 30
  v0 <- 4
  AR <- 1
  T <- 20
  
  N <- 50
  zparm_ix <- 3
  
  ### NOTE: The data is assumed to have been generated from ../emulate/emulate_EP_MNIW.R before this. It is advised to run it to generate both the data and the emulation before beginning this step.

  data.dir <- file.path(top.dir,
			sprintf("data/LV_train_Zparms%d_u%dv%d_AR%d_T%d", 
				zparm_ix, u0, v0, AR, T))
  
  # plot the prey curves over time as well as the hare points.
  Yall.file.name <- file.path(top.dir, 
			      "vignettes/predator_prey/calibration",
			      sprintf("Y_train_all_Zparms%d.csv", zparm_ix))
  
  if (!file.exists(Yall.file.name)) {
      #TODO: write a function (?) to join the dataframes across time into one. This is obviously not too big.
      Y_all <- data.frame("t"=numeric(),
                          "eta_ix"=numeric(),
                          "V1"=numeric(),
                          "V2"=numeric())
  
      for (t in 0:T) {
          # read in the file
          Yt <- read.fst(file.path(data.dir, sprintf("LV_train_Y%d.fst", t)))
  
          # Assign eta and t
          Yt$t <- t
          Yt$eta_ix <- seq(1, N)
  
          # append to data_all
          Y_all <- rbind(Y_all, Yt)
      }
  
      # rename to something better.
      colnames(Y_all)[colnames(Y_all) == "V1"] <- "ut"
      colnames(Y_all)[colnames(Y_all) == "V2"] <- "vt"
  
      zreal_ <- log(zreal[, c("Hare", "Lynx")])
      colnames(zreal_) <- c("ut", "vt")
      zreal_$t <- seq(0,T)
      zreal_$eta_ix <- 0
  
      Y_all <- rbind(Y_all, zreal_)
  
      # save it somewhere. If it already exists. Read it back in.
      write.csv(Y_all, Yall.file.name, row.names=FALSE)
  } else {
      Y_all <- read.csv(Yall.file.name)
  }
  
  Y_all$eta_ix <- factor(Y_all$eta_ix)
  
  #TODO: plot the curves and the real data
  up <- ggplot(data=Y_all, aes(x=t, y=ut, color=eta_ix, group=eta_ix)) + 
  	geom_line(data=subset(Y_all, eta_ix != 0)) +
  	geom_point(data=subset(Y_all, eta_ix == 0), color="black", size=5) +
  	geom_point(data=subset(Y_all, eta_ix == 0), color="green", size=2) +
  	theme_bw() + theme(legend.position = "none") + scale_color_viridis_d(option="D")
  
  ggsave(up, filename=file.path(top.dir,
				"vignettes/predator_prey/calibration",
				"LV_real_data_vs_sims_ut.png"))
  
  vp <- ggplot(data=Y_all, aes(x=t, y=vt, color=eta_ix, group=eta_ix)) + 
  	geom_line(data=subset(Y_all, eta_ix != 0)) +
  	geom_point(data=subset(Y_all, eta_ix == 0), color="black", size=5) +
  	geom_point(data=subset(Y_all, eta_ix == 0), color="orange", size=2) +
  	theme_bw() + theme(legend.position = "none") + scale_color_viridis_d(option="A")
  
  ggsave(vp, filename=file.path(top.dir,
				"vignettes/predator_prey/calibration",
				"LV_real_data_vs_sims_vt.png"))
}

{
  #####
  # Calibration run - Lotka-Volterra
  #####

  # This script tests calibration for an artificially-generated z to be calibrated.
  library(fields)
  library(matrixsampling)
  library(matrixStats)
  library(lhs)
  library(deSolve)
  library(ReacTran)
  library(readr)
  library(fst)
  
  #setwd("C:/Users/hunan/OneDrive/Documents/0_UCLA/Banerjee_GSR/projects/MNIW_Testing_temp/tests/testthat/calibration")
  
  top.dir <- here()
  
  source(file.path(top.dir, "R/MNIW.R"))
  source(file.path(top.dir, "R/big_data_set.R"))
  source(file.path(top.dir, "R/FFBS.R"))
  source(file.path(top.dir, "R/pds.R"))
  source(file.path(top.dir, "R/seed_obj.R"))
  source(file.path(top.dir, "R/calibrator.R"))
  
  #TODO: Smaller grid (e.g. 10x10 or 20x20 to make sure the spread goes far enough.)
  
  AR <- 1
  T <- 20
  times <- seq(0, T)
  f_times <- seq(AR, T)
  
  u0 <- 30
  v0 <- 4
  
  L_burn <- 0#5000
  L <- 20000
  
  zparm_ix <- 3
  
  train.data.dir <- file.path(top.dir, sprintf("data/LV_train_Zparms%d_u%dv%d_AR%d_T%d",
  					     zparm_ix, u0, v0, AR, T))
  
  Sigma_mode <- "IW"
  calib.out.dir <- file.path(top.dir,
			    "vignettes/predator_prey/calibration", 
  			 sprintf("LV_Zparms%d_zreal_u%dv%d_AR%d_%dT_%dsample_calib_%s", 
  				 zparm_ix, u0, v0, AR, T, L, Sigma_mode))
  em.out.dir <- calib.out.dir
  
  dir.create(calib.out.dir, showWarnings=FALSE)
  
  # Read in emulation files. These directories should already exist.
  eta_train <- read.csv(file.path(train.data.dir, "parameters_train.csv"))
  
  em.in.dir <- file.path(top.dir,
			 "vignettes/predator_prey/emulation",
			 sprintf("LV_train_Zparms%d_u%dv%d_AR%d_T%d_%s",
  				zparm_ix, u0, v0, AR, T, Sigma_mode))
  
  out.format <- "fst"
  
  
  set.seed(98765)
  
  mean_list <- c(0.0, 0.0, -3, 0.0, -3)
  sd_list <- c(0.5, 0.5, 0.5, 0.5, 0.5)
  
  # What am I doing now?
  logd_rho_eta <- function(rho, eta) {
      l_eta <- length(eta)
  
      logd <- 0
  
      if (!is.null(rho)) logd <- logd - (log(rho) - mean_list[1])^2/(2 * sd_list[1]^2) -
  	                       log(rho) - log(sd_list[1]) - log(2 * pi)/2
  
      for (i in 1:l_eta) {
          # The log-density.
          logd <- logd - (log(eta[i]) - mean_list[i+1])^2/(2 * sd_list[i+1]^2) -
  		log(eta[i]) - log(sd_list[i+1]) - log(2 * pi)/2
      }
  
      return(logd)
  }
  
  # rho and eta transformation functions
  rho_eta_2_Rspace <- function(rho, eta) {
      phi <- log(c(rho, eta))
      phi <- (phi - mean_list)/sd_list
  
      return(phi)
  }
  
  Rspace_2_rho_eta <- function(phi, l_rho, l_eta) {
      phi <- phi * sd_list + mean_list
      rho_eta <- exp(phi)
  
      return(rho_eta)
  }
  
  logJ_rho_eta <- function(rho, eta) {
      l_rho <- length(rho)
      l_eta <- length(eta)
  
      logJ <- -log(rho) -log(sd_list[1])
  
      for (i in 1:l_eta) logJ <- logJ - log(eta[i]) - log(sd_list[i+1])
  
      return(logJ)
  }
  
  LV <- function(t, y, parms) {
      u <- matrix(nrow = 1, ncol = 1, data = y[1])
      v <- matrix(nrow = 1, ncol = 1, data = y[2])
  
      with (as.list(parms),{
          du <- eta1 * u - eta2 * u * v
          dv <- -eta3 * v + eta4 * u * v
          list(c(du, dv))
      })
  }
  
  yini <- c(u0, v0)
  
  out_2_y <- function(ode_out) {
      return(log(ode_out))
  }
  
  N <- nrow(eta_train)
  p <- AR * 2
  
  n0z <- p + 3
  d0z <- 2*p + p
  b <- 1.0
  
  nycol <- 2
  m0z <- rep(1, nycol)
  M0z <- diag(nycol)
  
  #rho_true <- 1.5
  rho0 <- 1.5#rho_true
  DT <- read.fst(file.path(em.in.dir, sprintf("FF_out_Dt_uptri_%d.fst", T - AR + 1)))
  DT <- tri.to.sym(unlist(DT[2,]), nycol, lower=FALSE)
  
  nT <- read.csv(file.path(em.in.dir, sprintf("FF_out_nt_%d.csv", T - AR + 1)), header=FALSE)
  nT <- nT[length(nT)]
  
  ESigma <- DT/(nT - 3)
  
  tau2_true <- rep(-1, T+1)
  for (t in 1:(T+1)) tau2_true[t] <- 1/rgamma(1, shape=b^(t-1) * n0z, rate=b^(t-1) * d0z)
  
  #z_path <- file.path(calib.out.dir, "../lynx_hare_df.csv")
  z_path <- file.path(top.dir, "data/lynx_hare_df.csv")
  z <- log(as.matrix(read.csv(z_path)[,c("Hare","Lynx")]))
  
  K <- 1
  
  T_real <- 20
  
  #Y_files <- file.path(em.out.dir, sprintf("SIR_train_Ycal%d.%s", seq(1, T), out.format))
  Y_files <- file.path(train.data.dir, sprintf("LV_train_Y%d.%s", seq(1, T_real), out.format))
  F_files <- file.path(train.data.dir, sprintf("LV_train_F%d.%s", seq(1, T_real), out.format))
  
  
  # define functions for calibration
  # 
  
  # Generate AR2 Ft and return
  gen_F <- function(t, eta, ymat) {
      return(c(ymat[t-1,]))
  }
  
  # Experimenting with another random eta, instead of one that exists in the data.
  eta0 <- c(exp(rnorm(1, mean=mean_list[2], sd=sd_list[2])),
                exp(rnorm(1, mean=mean_list[3], sd=sd_list[3])),
                exp(rnorm(1, mean=mean_list[4], sd=sd_list[4])),
                exp(rnorm(1, mean=mean_list[5], sd=sd_list[5])))
  names(eta0) <- c("eta1", "eta2", "eta3", "eta4")
  
  #sigma2_file <- file.path(em.in.dir, sprintf("Bsample%d_out_sigma2_%d.csv", L, L))
  
  dists <- sqrt(rdist(eta_train, eta_train))
  
  deff <- 0.5 * max(dists)
  beta <- 3/deff
  
  Vt <- exp(-beta * dists)
  
  nT <- unlist(read.csv(file.path(em.in.dir, sprintf("FF_out_nt_%d.csv",
  						   max(f_times) - min(f_times) + 1)),
  			  header=FALSE))
  
  names(nT) <- NULL
  nT <- nT[[length(nT)]]
     
  #smooth.out.head <- file.path(em.in.dir, sprintf("Bsample0_out"))
  smooth.out.head <- file.path(em.in.dir, sprintf("BS_out"))
  st_files=paste0(smooth.out.head, "_st_", seq(1, T_real - AR + 1),".", out.format)
  St_files=paste0(smooth.out.head, "_St_uptri_", seq(1, T_real - AR + 1),".", out.format)
  
  t0 <- Sys.time()
  
  # May want to add initial y and other components.
  # Run calibrator with expressions.
  calib_out <- calibrator(z[2:(T_real+1),], Y_files, F_files, N, p,
                          n0z, d0z, b, m0z, M0z, rho0, eta0, eta_train,
                          logd_rho_eta, rho_eta_2_Rspace, Rspace_2_rho_eta, logJ_rho_eta,
                          Vt, beta,
                          NULL,
                          LV, NULL, 2, yini, out_2_y, gen_F,
                          y0 = NULL, Btmat0 = NULL,
                          Theta_files = NULL,
                          Sigma_file = NULL, sigma2_file = NULL, R = NULL,
                          seed_list_TS = c(12345), seed_list_G = c(98765),
                          nT = nT, dT = NULL, DT = DT,
                          st_files = st_files, St_files = St_files, AR = AR,
                          L = L, L_burn = L_burn, eps3 = 0.1,
                          rho_true = NULL,
                          out.head = file.path(calib.out.dir, "calib_out"),
                          save.dir = calib.out.dir,
                          out.format = out.format,
                          max_iter_size = L,
                          verbose = TRUE)
  
  
  print(sprintf("Calibrator run time for %d iterations with %d burn-in: %.3f", L, L_burn, Sys.time() - t0))

  #####
  # Plot calibration results
  #####

  library(dplyr)
  library(fst)

  S <- 2  
  p <- AR * S
  
  n0z <- p + 3
  d0z <- 2*p + p
  b <- 1.0
  
  zparm_ix <- 3
  
  calib.out.dir <- sprintf(file.path(top.dir,
				     "vignettes/predator_prey/calibration",
				     "LV_Zparms%d_zreal_u30v4_AR%d_%dT_%dsample_calib_IW"), 
  			 zparm_ix, AR, T, L)
  
  probs <- read.fst(file.path(calib.out.dir, "calib_out_metprobs.fst"))
  accepts <- read.fst(file.path(calib.out.dir, "calib_out_metaccept.fst"))
  
  print(sprintf("Expected proportion of metropolis jumps accepted: %.5f", sum(probs > 0.2)/nrow(probs)))
  print(sprintf("Actual proportion of metropolis jumps accepted: %.5f", sum(accepts > 0)/nrow(accepts)))
  
  etas <- read.fst(file.path(calib.out.dir, "calib_out_eta.fst"))
  
  medians <- rep(-1, 4)
  
  names(etas) <- paste0("eta", seq(1,4))
  
  eta_quants <- quantile(etas$eta1, probs=c(0.025, 0.5, 0.975))
  print(eta_quants)
  medians[1] <- eta_quants[2]
  
  p_etak <- ggplot(data=etas, aes(x=eta1)) + 
  	      geom_histogram(fill="orange", color="black") + theme_bw() +
  	      geom_vline(xintercept = eta_quants[1], linetype="dashed", color="black") +
  	      geom_vline(xintercept = eta_quants[2], linetype="solid", color="black") +
  	      geom_vline(xintercept = eta_quants[3], linetype="dashed", color="black") +
  	      xlim(-0.2, max(etas$eta1) + 0.2)
  
  ggsave(filename=file.path(top.dir,
			    "vignettes/predator_prey/calibration",
			    "LV_eta1_hist_zreal.png"), p_etak)
  
  eta_quants <- quantile(etas$eta2, probs=c(0.025, 0.5, 0.975))
  print(eta_quants)
  medians[2] <- eta_quants[2]
  p_etak <- ggplot(data=etas, aes(x=eta2)) + 
  	      geom_histogram(fill="lightgreen", color="black") + theme_bw() +
  	      geom_vline(xintercept = eta_quants[1], linetype="dashed", color="black") +
  	      geom_vline(xintercept = eta_quants[2], linetype="solid", color="black") +
  	      geom_vline(xintercept = eta_quants[3], linetype="dashed", color="black") +
  	      xlim(-0.2, max(etas$eta2) + 0.2)
  
  ggsave(filename=file.path(top.dir,
			    "vignettes/predator_prey/calibration",
			    "LV_eta2_hist_zreal.png"), p_etak)
  
  eta_quants <- quantile(etas$eta3, probs=c(0.025, 0.5, 0.975))
  print(eta_quants)
  medians[3] <- eta_quants[2]
  p_etak <- ggplot(data=etas, aes(x=eta3)) + 
  	      geom_histogram(fill="purple", color="black") + theme_bw() +
  	      geom_vline(xintercept = eta_quants[1], linetype="dashed", color="black") +
  	      geom_vline(xintercept = eta_quants[2], linetype="solid", color="black") +
  	      geom_vline(xintercept = eta_quants[3], linetype="dashed", color="black") +
  	      xlim(-0.2, max(etas$eta3) + 0.2)
  
  ggsave(filename=file.path(top.dir,
			    "vignettes/predator_prey/calibration",
			    "LV_eta3_hist_zreal.png"), p_etak)
  
  eta_quants <- quantile(etas$eta4, probs=c(0.025, 0.5, 0.975))
  print(eta_quants)
  medians[4] <- eta_quants[2]
  p_etak <- ggplot(data=etas, aes(x=eta4)) + 
  	      geom_histogram(fill="lightyellow", color="black") + theme_bw() +
  	      geom_vline(xintercept = eta_quants[1], linetype="dashed", color="black") +
  	      geom_vline(xintercept = eta_quants[2], linetype="solid", color="black") +
  	      geom_vline(xintercept = eta_quants[3], linetype="dashed", color="black") +
  	      xlim(-0.2, max(etas$eta4) + 0.2)
  
  ggsave(filename=file.path(top.dir,
			    "vignettes/predator_prey/calibration",
			    "LV_eta4_hist_zreal.png"), p_etak)
  
  
  write.csv(as.data.frame(medians), 
	    file = file.path(top.dir,
	                     "vignettes/predator_prey/calibration",
			     "zreal_median_etas.csv"), row.names = FALSE)
  
  # Scatter plot tau2's over time
  tau2t_df <- read.fst(file.path(calib.out.dir, "calib_out_tau2t.fst"))
  
  # Plot bias values over time for both ut and vt.
  ut_df <- read.fst(file.path(calib.out.dir, "calib_out_ut_1.fst"))
  colnames(ut_df)[colnames(ut_df) == c("V1", "V2")] <- c("ut", "vt")
  
  p_ut <- ggplot(data=ut_df, aes(x=t, y=ut)) + geom_point(color="black") +
  	theme_bw()
  
  ggsave(filename=file.path(top.dir,
			    "vignettes/predator_prey/calibration",
			    "bias_prey_zreal.png"), p_ut)
  
  p_vt <- ggplot(data=ut_df, aes(x=t, y=vt)) + geom_point(color="black") +
  	theme_bw()
  
  ggsave(filename=file.path(top.dir,
			    "vignettes/predator_prey/calibration",
			    "bias_pred_zreal.png"), p_vt)
  
}
