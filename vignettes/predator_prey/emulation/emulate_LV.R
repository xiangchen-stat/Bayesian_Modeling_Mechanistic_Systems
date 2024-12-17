# Lotka-Volterra emulation
{
  library(here)
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
  
  #locs <- expand.grid(y=1:Ny, x=1:Nx)
  #locs <- locs[,c("x", "y")]
  #locs = as.numeric(rownames(locs))
  
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
  ffbs.out.dir <- file.path(top.dir, sprintf("vignettes/predator_prey/emulation/LV_train_Zparms%d_u%dv%d_AR%d_T%d_%s", 
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
  			  sprintf("vignettes/predator_prey/emulation/LV_predict_Zparms%d_u%dv%d_AR%d_T%d_%s", 
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
		       "vignettes/predator_prey/emulation",
		       sprintf("LV_train_Zparms%d_u30v4_AR%d_T%d_IW", zparm_ix, AR, T))
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
  
  ggsave(c_p, file=file.path(emu.dir, "corr_hist.png"))
  
  S11_p <- ggplot(data=data.frame(Sigma_11 = Sigmas[1,1,]), aes(x=Sigma_11)) + geom_histogram(fill="orange", color="black") +
  	theme_bw() + geom_vline(xintercept=Sigma_11_quantiles[2], color="black", linetype="solid") +
  	geom_vline(xintercept=Sigma_11_quantiles[1], color="black", linetype="dashed") +
  	geom_vline(xintercept=Sigma_11_quantiles[3], color="black", linetype="dashed")
  
  ggsave(S11_p, file=file.path(emu.dir, "Sigma11_hist.png"))
  
  S12_p <- ggplot(data=data.frame(Sigma_12 = Sigmas[1,2,]), aes(x=Sigma_12)) + geom_histogram(fill="lightgreen", color="black") +
  	theme_bw() + geom_vline(xintercept=Sigma_12_quantiles[2], color="black", linetype="solid") +
  	geom_vline(xintercept=Sigma_12_quantiles[1], color="black", linetype="dashed") +
  	geom_vline(xintercept=Sigma_12_quantiles[3], color="black", linetype="dashed")
  
  ggsave(S12_p, file=file.path(emu.dir, "Sigma12_hist.png"))
  
  S22_p <- ggplot(data=data.frame(Sigma_22 = Sigmas[2,2,]), aes(x=Sigma_22)) + geom_histogram(fill="lightyellow", color="black") +
  	theme_bw() + geom_vline(xintercept=Sigma_22_quantiles[2], color="black", linetype="solid") +
  	geom_vline(xintercept=Sigma_22_quantiles[1], color="black", linetype="dashed") +
  	geom_vline(xintercept=Sigma_22_quantiles[3], color="black", linetype="dashed")
  
  ggsave(S22_p, file=file.path(emu.dir, "Sigma22_hist.png"))
  
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
		       "vignettes/predator_prey/emulation",
		       sprintf("LV_train_Zparms%d_u%dv%d_AR%d_T%d_%s", 
			       zparm_ix, u0, v0, AR, T, Sigma_mode))
  pred.dir <- file.path(top.dir,
		       "vignettes/predator_prey/emulation",
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
      
      ggsave(file.path(pred.dir, sprintf("LV_plots_pred_%d.png", eta_ix_plot)), plot=p_Yt)
  
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
      
      ggsave(file.path(pred.dir, sprintf("LV_plots_ut_qq_%d.png", eta_ix_plot)), plot=p_qq_ut)
      ggsave(file.path(pred.dir, sprintf("LV_plots_vt_qq_%d.png", eta_ix_plot)), plot=p_qq_vt)
  }
}

