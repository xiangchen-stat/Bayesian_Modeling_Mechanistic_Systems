read_big_csv_quick <- function(filename, rows = c(1,Inf), cols = c(1,Inf), header=TRUE) {
  
  if (is.finite(cols[2])) {
    # headers <- as.character(as.vector(read_csv(filename,
    #                                            col_names = FALSE, n_max = 1,
    #                                            col_select = cols[1]:cols[2],
    #                                            show_col_types=FALSE)))
    #TODO: How to suppress progress bar here?
    out <- read_csv(file = filename, col_names = FALSE,
                    skip = rows[1] - 1 + as.integer(header), 
                    n_max = rows[2] - rows[1] + 1,
                    col_select = cols[1]:cols[2],
                    show_col_types = FALSE)
    
  } else {
    #TODO: something is going on with this read_csv that is throwing an error, but when I run it in the R shell, it works fine. Oddly, neglecting this option does not throw an error. 
    # headers <- read_csv(filename, col_names = FALSE, n_max = 1,
    #                     show_col_types=FALSE)
    out <- read_csv(file = filename, col_names = FALSE,
                    skip = rows[1] - 1 + as.integer(header), 
                    n_max = rows[2] - rows[1] + 1,
                    show_col_types = FALSE)
  }
  
  # if (header) colnames(out) <- headers
  return(out)
}


## MNIW----

#' one time step for the Forward Filter. Computes the FF parameters given the data at the relevant time step and the relevant parameters from the last time step.
#' @param Yt The data matrix for epoch t.
#' @param Ft The matrix of covariates F_t.
#' @param Gt G_t, the beta transition matrix.
#' @param mt_1 mean of beta_{t-1} | D_{t-1}
#' @param Mt_1 left-covariance matrix of beta_{t-1} | D_{t-1}
#' @param Wt left-covariance matrix of the noise parameter of beta_{t-1} | D_{t-1}
#' @param Vt left-covariance matrix of the noise parameter of Y
#' @param nt_1 the shape parameter, or degrees of freedom, of the right-covariance matrix Sigma | D_{t-1}
#' @param Dt_1 the scale matrix of the right-covariance matrix Sigma | D_{t-1}
#' @param delta The right-variance matrix discount factor 
#' @returns The mean and covariance matrices m_t and C_t of one filtering step, the updated inverse-Wishart parameters a_t and B_t for the right-covariance matrix, plus other relevant parameters.
#' @export
# FF_1step_R_naiive <- function(Yt, Ft, Gt, Wt, Vt, mt_1, Mt_1, nt_1, Dt_1, delta = 1.0){
#     N <- dim(Yt)[1]
#     at <- Gt %*% mt_1
#     At <- Gt %*% Mt_1 %*% t(Gt) + Wt
#     qt <- Ft %*% at
#     Qt <- Ft %*% At %*% t(Ft) + Vt
#     Qtinv <- solve(Qt)
#     mt <- at + At %*% t(Ft) %*% Qtinv %*% (Yt - qt)
#     Mt <- At - At %*% t(Ft) %*% Qtinv %*% Ft %*% At
#     nt <- nt_1 + N
#     Dt <- Dt_1 + t(Yt - qt) %*% Qtinv %*% (Yt - qt)
#     
#     out <- list(nt = nt, Dt = Dt, at = at, At = At,
#                 mt = mt, Mt = Mt, delta = delta)
#     return(out)
# }




# TODO reduce return parameters
#' one time step for the Forward Filter. Computes the FF parameters given the data at the relevant time step and the relevant parameters from the last time step.
#' @param Yt The data matrix for epoch t.
#' @param Ft The matrix of covariates F_t.
#' @param Gt G_t, the beta transition matrix.
#' @param mt_1 mean of beta_{t-1} | D_{t-1}
#' @param Mt_1 left-covariance matrix of beta_{t-1} | D_{t-1}
#' @param Wt left-covariance matrix of the noise parameter of beta_{t-1} | D_{t-1}
#' @param Vt left-covariance matrix of the noise parameter of Y
#' @param nt_1 the shape parameter, or degrees of freedom, of the right-covariance matrix Sigma | D_{t-1}
#' @param Dt_1 the scale matrix of the right-covariance matrix Sigma | D_{t-1}
#' @param delta The right-variance matrix discount factor 
#' @returns The mean and covariance matrices m_t and C_t of one filtering step, the updated inverse-Wishart parameters a_t and B_t for the right-covariance matrix, plus other relevant parameters.
#' @export
# FF_1step_R <- function(Yt, Ft, Gt, Wt, Vt, mt_1, Mt_1, nt_1, Dt_1, delta = 1.0){
#   N <- dim(Yt)[1]
#   at <- Gt %*% mt_1
#   At <- tcrossprod(Gt %*% Mt_1, Gt) + Wt
#   qt <- Ft %*% at
#   FtAt <- Ft %*% At
#   Qt <- tcrossprod(FtAt, Ft) + Vt
#   Qtinv <- solve(Qt)
#   AttFtQtinv <- tcrossprod(At, Ft)%*% Qtinv
#   Yt_qt <- Yt - qt
#   mt <- at + AttFtQtinv %*% Yt_qt
#   Mt <- At - AttFtQtinv %*% FtAt
#   nt <- nt_1 + N
#   Dt <- Dt_1 + crossprod(Yt_qt, Qtinv) %*% Yt_qt
#   
#   out <- list(nt = nt, Dt = Dt, at = at, At = At,
#               mt = mt, Mt = Mt, delta = delta)
#   return(out)
# }

# TODO FF_cpp
#' Forward Filter. Computes the FF parameters given the data at the relevant time step and the relevant parameters from the last time step.
#' @param Y The data matrix.
#' @param Ft The matrix of covariates F_t.
#' @param Gt G_t, the beta transition matrix.
#' @param m0 mean of beta_{t-1} | D_{t-1}
#' @param M0 left-covariance matrix of beta_{t-1} | D_{t-1}
#' @param Wt left-covariance matrix of the noise parameter of beta_{t-1} | D_{t-1}
#' @param Vt left-covariance matrix of the noise parameter of Y
#' @param n0 the shape parameter, or degrees of freedom, of the right-covariance matrix Sigma | D_{t-1}
#' @param D0 the scale matrix of the right-covariance matrix Sigma | D_{t-1}
#' @param delta The right-variance matrix discount factor 
#' @returns The mean and covariance matrices m_t and C_t of one filtering step, the updated inverse-Wishart parameters a_t and B_t for the right-covariance matrix, plus other relevant parameters.
#' @export
FF <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, delta = 1.0){
  out <- list()
  out_name <- c()
  # initialize matirces
  Ft <- F_ls
  Gt <- G_ls
  Wt <- W_ls
  Vt <- V_ls
  
  for (i in 1:nT) {
    if(i %% round(0.1*nT) == 0){
      print(paste("FF:", i, "/", nT))
      print(Sys.time())
    }
    # Check if F_ls, G_ls, W_ls, V_ls are changing over time
    # If they are list, change accordingly, otherwise keep unchanged
    if(is.list(F_ls)){
      Ft <- as.matrix(F_ls[[i]])
    }
    if(is.list(G_ls)){
      Gt <- as.matrix(G_ls[[i]])
    }
    if(is.list(W_ls)){
      Wt <- as.matrix(W_ls[[i]])
    }
    if(is.list(V_ls)){
      Vt <- as.matrix(V_ls[[i]])
    }
    
    Yt <- Y[[i]]
    # Set T-1 values
    if (i == 1) {
      mt_1 = m0
      Mt_1 = M0
      nt_1 = n0
      Dt_1 = D0
    } else{
      mt_1 = one_step$mt
      Mt_1 = one_step$Mt
      nt_1 = one_step$nt
      Dt_1 = one_step$Dt
    }
    
    # perform 1 step forward
    one_step <- FF_1step_cpp(Yt = Yt, Ft = Ft, Gt = Gt, 
                             Wt = Wt, Vt = Vt, 
                             mt_1 = mt_1, Mt_1 = Mt_1, 
                             nt_1 = nt_1, Dt_1 = Dt_1, delta = delta)
    one_step <- one_step[c("nt", "Dt", "at", "At", "mt", "Mt")] # extract parameters
    out[[i]] <- one_step # save results
  }
  
  out_name <- paste("T", seq(1:nT), sep = "") # accumulate output list names
  names(out) <- out_name
  out$prior <- list(m0 = m0, M0 = M0, n0 = n0, D0 = D0)
  return(out)
}



# FF_bigdata_R_naiive_no_extract_wrong <- function(Y_ls, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, 
#                                 fnrow, fncol, bnrow, bncol, path_out, delta = 1.0){
#   # initialize matirces
#   Ft <- F_ls
#   Gt <- G_ls
#   Wt <- W_ls
#   Vt <- V_ls
#   
#   # generate index
#   ind <- generate.grid.rowsnake(fnrow = fnrow, fncol = fncol, bnrow = bnrow, bncol = bncol)
#   n_b <- dim(ind)[1]
#   # n_r <- fnrow / bnrow
#   # n_c <- fncol / bncol
# 
#   # loop through files
#   for (f in 1:length(Y_ls)) {
#     out <- list()
#     Y_full <- as.matrix(read.csv(Y_ls[f]))
#     # Check if F_ls, G_ls, W_ls, V_ls are changing over file
#     # If they are list, change accordingly, otherwise keep unchanged
#     if(is.list(F_ls)){
#       Ft <- as.matrix(F_ls[[f]])
#     }
#     if(is.list(G_ls)){
#       Gt <- as.matrix(G_ls[[f]])
#     }
#     if(is.list(W_ls)){
#       Wt <- as.matrix(W_ls[[f]])
#     }
#     if(is.list(V_ls)){
#       Vt <- as.matrix(V_ls[[f]])
#     }
#     
#     # loop through blocks
#     for (i in 1:n_b) {
#       Yt <- Y_full[ind[i,2]:ind[i,3], ind[i,4]:ind[i,5]]
#       # Set T+1 values
#       if (i == 1 && f == 1) {
#         mt_1 = m0
#         Mt_1 = M0
#         nt_1 = n0
#         Dt_1 = D0
#       } else{
#         mt_1 = one_step$mt
#         Mt_1 = one_step$Mt
#         nt_1 = one_step$nt
#         Dt_1 = one_step$Dt
#       }
#       
#       # perform 1 step forward
#       one_step <- FF_1step_cpp(Yt = Yt, Ft = Ft, Gt = Gt, 
#                                Wt = Wt, Vt = Vt, 
#                                mt_1 = mt_1, Mt_1 = Mt_1, 
#                                nt_1 = nt_1, Dt_1 = Dt_1, delta = delta)
#       one_step <- one_step[c("nt", "Dt", "at", "At", "mt", "Mt")] # extract parameters
#       if(i == 1){
#         out$nt <- one_step$nt
#         out$Dt <- one_step$Dt
#         out$at <- one_step$at
#         out$At <- one_step$At
#         out$mt <- one_step$mt
#         out$Mt <- one_step$Mt
#       } else{
#         out$nt <- cbind(out$nt, one_step$nt)
#         out$Dt <- cbind(out$Dt, one_step$Dt)
#         out$at <- cbind(out$at, one_step$at)
#         out$At <- cbind(out$At, one_step$At)
#         out$mt <- cbind(out$mt, one_step$mt)
#         out$Mt <- cbind(out$Mt, one_step$Mt)
#       }
#     }
#     
#     # after looping each file, save results
#     write.csv(out$nt, file = paste0(path_out, "/output_FF_nt_", f, ".csv", 
#                                          sep = ""), row.names = F)
#     write.csv(out$Dt, file = paste0(path_out, "/output_FF_Dt_", f, ".csv", 
#                                          sep = ""), row.names = F)
#     write.csv(out$mt, file = paste0(path_out, "/output_FF_mt_", f, ".csv", 
#                                          sep = ""), row.names = F)
#     write.csv(out$Mt, file = paste0(path_out, "/output_FF_MMt_", f, ".csv", 
#                                          sep = ""), row.names = F)
#     write.csv(out$at, file = paste0(path_out, "/output_FF_at_", f, ".csv", 
#                                          sep = ""), row.names = F)
#     write.csv(out$At, file = paste0(path_out, "/output_FF_AAt_", f, ".csv", 
#                                          sep = ""), row.names = F)
#     print(paste0("Finish ", f, "/", length(Y_ls)))
#     print(Sys.time())
#   }
# }

# TODO FF_cpp
#' Forward Filter. Computes the FF parameters given the data at the relevant time step and the relevant parameters from the last time step.
#' @param Y The data matrix.
#' @param Ft The matrix of covariates F_t.
#' @param Gt G_t, the beta transition matrix.
#' @param m0 mean of beta_{t-1} | D_{t-1}
#' @param M0 left-covariance matrix of beta_{t-1} | D_{t-1}
#' @param Wt left-covariance matrix of the noise parameter of beta_{t-1} | D_{t-1}
#' @param Vt left-covariance matrix of the noise parameter of Y
#' @param n0 the shape parameter, or degrees of freedom, of the right-covariance matrix Sigma | D_{t-1}
#' @param D0 the scale matrix of the right-covariance matrix Sigma | D_{t-1}
#' @param delta The right-variance matrix discount factor 
#' @returns The mean and covariance matrices m_t and C_t of one filtering step, the updated inverse-Wishart parameters a_t and B_t for the right-covariance matrix, plus other relevant parameters.
#' @export
FF_bigdata_R <- function(Y_ls, F_ls, G_ls, W_ls, V_ls, 
                                m0, M0, n0, D0, nT, 
                                fnrow, fncol, bnrow, bncol, path_out, 
                                delta = 1.0, verbose = FALSE){
  # generate index
  ind <- generate.grid.rowsnake(fnrow = fnrow, fncol = fncol, bnrow = bnrow, bncol = bncol)
  n_b <- dim(ind)[1] # number of blocks
  # n_r <- fnrow / bnrow
  # n_c <- fncol / bncol
  
  # initialize matirces if Gt and Wt never changes
  if(length(G_ls) == 1){
    Gt <- as.matrix(read_big_csv_quick(filename = G_ls))
  }
  
  if(length(W_ls) == 1){
    Wt <- as.matrix(read_big_csv_quick(filename = W_ls))
  }
  
  # loop through nT time/files
  for (f in 1:nT) {
    out <- list()
    
    # Check if G_ls, W_ls are changing over file
    # If they are array of filename, change accordingly, otherwise keep unchanged
    if(length(G_ls) != 1){
      Gt <- as.matrix(read_big_csv_quick(filename = G_ls[f]))
    }
    
    if(length(W_ls) != 1){
      Wt <- as.matrix(read_big_csv_quick(filename = W_ls[f]))
    }
    
    # indicator to avoid read Vt each time
    row_vt_old <- 0

    # loop through blocks
    for (i in 1:n_b) {
      # Extract Yt, Ft, Vt from files
      Yt <- as.matrix(read_big_csv_quick(filename = Y_ls[f], 
                                       rows = c(ind[i,2], ind[i,3]), 
                                       cols = c(ind[i,4], ind[i,5])))
      # extract Ft
      if(length(F_ls) != 1){
        # read ar1
        Ft <- as.matrix(read_big_csv_quick(filename = F_ls[f], 
                                           rows = c(ind[i,2], ind[i,3]),
                                           cols = c(ind[i,4], ind[i,5])))
        # try read ar2
        try(expr = {
            Ft_ar2 <- as.matrix(read_big_csv_quick(filename = F_ls[f], 
                                                   rows = c(ind[i,2], ind[i,3]),
                                                   cols = c((fncol + ind[i,4]), (fncol + ind[i,5]))));
            Ft <- cbind(Ft, Ft_ar2);
          }, 
          silent = T)
      } else{
        # read ar1
        Ft <- as.matrix(read_big_csv_quick(filename = F_ls, 
                                           rows = c(ind[i,2], ind[i,3]),
                                           cols = c(ind[i,4], ind[i,5])))
        # try read ar2
        try(expr = {
          Ft_ar2 <- as.matrix(read_big_csv_quick(filename = F_ls, 
                                                 rows = c(ind[i,2], ind[i,3]),
                                                 cols = c((fncol + ind[i,4]), (fncol + ind[i,5]))));
          Ft <- cbind(Ft, Ft_ar2);
        }, 
        silent = T)
      }
      
      # extract Vt
      row_vt_new <- ind[i,2]
      if(row_vt_new != row_vt_old){
        # extract Vt
        if(length(V_ls) != 1){
          Vt <- as.matrix(read_big_csv_quick(filename = V_ls[f], 
                                             rows = c(ind[i,2], ind[i,3]), 
                                             cols = c(ind[i,2], ind[i,3])))
        } else{
          Vt <- as.matrix(read_big_csv_quick(filename = V_ls, 
                                             rows = c(ind[i,2], ind[i,3]), 
                                             cols = c(ind[i,2], ind[i,3])))
        }
      }
      row_vt_old <- row_vt_new # update indicator for Ft and Vt 1st row
      
      
      # Set T-1 values
      if (i == 1 && f == 1) {
        mt_1 = m0
        Mt_1 = M0
        nt_1 = n0
        Dt_1 = D0
      } else{
        mt_1 = one_step$mt
        Mt_1 = one_step$Mt
        nt_1 = one_step$nt
        Dt_1 = one_step$Dt
      }
      
      # perform 1 step forward
      one_step <- FF_1step_cpp(Yt = Yt, Ft = Ft, Gt = Gt, 
                               Wt = Wt, Vt = Vt, 
                               mt_1 = mt_1, Mt_1 = Mt_1, 
                               nt_1 = nt_1, Dt_1 = Dt_1, delta = delta)
      # one_step <- one_step[c("nt", "Dt", "at", "At", "mt", "Mt")] # extract parameters
      if(i == 1){
        out$nt <- one_step$nt
        out$Dt <- one_step$Dt
        out$at <- one_step$at
        out$At <- one_step$At
        out$mt <- one_step$mt
        out$Mt <- one_step$Mt
      } else{
        out$nt <- cbind(out$nt, one_step$nt)
        out$Dt <- cbind(out$Dt, one_step$Dt)
        out$at <- cbind(out$at, one_step$at)
        out$At <- cbind(out$At, one_step$At)
        out$mt <- cbind(out$mt, one_step$mt)
        out$Mt <- cbind(out$Mt, one_step$Mt)
      }
    }
    
    # after looping one time/file, save results
    write.csv(out$nt, file = paste0(path_out, "/output_FF_nt_", f, ".csv", 
                                    sep = ""), row.names = F)
    write.csv(out$Dt, file = paste0(path_out, "/output_FF_Dt_", f, ".csv", 
                                    sep = ""), row.names = F)
    write.csv(out$mt, file = paste0(path_out, "/output_FF_mt_", f, ".csv", 
                                    sep = ""), row.names = F)
    write.csv(out$Mt, file = paste0(path_out, "/output_FF_MMt_", f, ".csv", 
                                    sep = ""), row.names = F)
    write.csv(out$at, file = paste0(path_out, "/output_FF_at_", f, ".csv", 
                                    sep = ""), row.names = F)
    write.csv(out$At, file = paste0(path_out, "/output_FF_AAt_", f, ".csv", 
                                    sep = ""), row.names = F)
    if(verbose == TRUE){
      print(paste0("Forward filter", f, "/", nT))
      print(Sys.time())
      # if(f %% 10 == 0){
      #   print(paste("Finish FF at time", f, "/", nT))
      #   print(Sys.time())
      # }
    }
  }
  
  # Save last Ft and Vt for BS 1st sample
  write.csv(Ft, file = paste0(path_out, "/output_FF_FT.csv", 
                              sep = ""), row.names = F)
  write.csv(Vt, file = paste0(path_out, "/output_FF_VT.csv", 
                              sep = ""), row.names = F)
  
  file_ls <- list()
  file_ls$nt_ls <- paste0(path_out, "/output_FF_nt_", seq(1:nT), ".csv", sep = "")
  file_ls$Dt_ls <- paste0(path_out, "/output_FF_Dt_", seq(1:nT), ".csv", sep = "")
  file_ls$mt_ls <- paste0(path_out, "/output_FF_mt_", seq(1:nT), ".csv", sep = "")
  file_ls$Mt_ls <- paste0(path_out, "/output_FF_MMt_", seq(1:nT), ".csv", sep = "")
  file_ls$at_ls <- paste0(path_out, "/output_FF_at_", seq(1:nT), ".csv", sep = "")
  file_ls$At_ls <- paste0(path_out, "/output_FF_AAt_", seq(1:nT), ".csv", sep = "")
  file_ls$FT <- paste0(path_out, "/output_FF_FT.csv", sep = "")
  file_ls$VT <- paste0(path_out, "/output_FF_VT.csv", sep = "")
  
  return(file_ls)
}

##sigma2R----
#' one time step for the Forward Filter. Computes the FF parameters given the data at the relevant time step and the relevant parameters from the last time step.
#' @param Yt The data matrix for epoch t.
#' @param Ft The matrix of covariates F_t.
#' @param Gt G_t, the beta transition matrix.
#' @param mt_1 mean of beta_{t-1} | D_{t-1}
#' @param Mt_1 left-covariance matrix of beta_{t-1} | D_{t-1}
#' @param Wt left-covariance matrix of the noise parameter of beta_{t-1} | D_{t-1}
#' @param Vt left-covariance matrix of the noise parameter of Y
#' @param nt_1 the shape parameter, or degrees of freedom, of the right-covariance matrix Sigma | D_{t-1}
#' @param Dt_1 the scale matrix of the right-covariance matrix Sigma | D_{t-1}
#' @param delta The right-variance matrix discount factor 
#' @returns The mean and covariance matrices m_t and C_t of one filtering step, the updated inverse-Wishart parameters a_t and B_t for the right-covariance matrix, plus other relevant parameters.
#' @export
# FF_1step_R_naiive_sigma2R <- function(Yt, Ft, Gt, Wt, Vt, mt_1, Mt_1, nt_1, Dt_1, Rinv, delta = 1.0){
#   N <- dim(Yt)[1]
#   S <- dim(Yt)[2] 
#   at <- Gt %*% mt_1
#   At <- Gt %*% Mt_1 %*% t(Gt) + Wt
#   qt <- Ft %*% at
#   Qt <- Ft %*% At %*% t(Ft) + Vt
#   Qtinv <- solve(Qt)
#   mt <- at + At %*% t(Ft) %*% Qtinv %*% (Yt - qt)
#   Mt <- At - At %*% t(Ft) %*% Qtinv %*% Ft %*% At
#   nt <- nt_1 + N * S / 2
#   Dt <- Dt_1 + 1/2 * matrixcalc::matrix.trace(t(Yt - qt) %*% Qtinv %*% (Yt - qt) %*% Rinv)
#   
#   out <- list(nt = nt, Dt = Dt, at = at, At = At,
#               mt = mt, Mt = Mt, delta = delta)
#   return(out)
# }


# TODO reduce return parameters
#' one time step for the Forward Filter. Computes the FF parameters given the data at the relevant time step and the relevant parameters from the last time step.
#' @param Yt The data matrix for epoch t.
#' @param Ft The matrix of covariates F_t.
#' @param Gt G_t, the beta transition matrix.
#' @param mt_1 mean of beta_{t-1} | D_{t-1}
#' @param Mt_1 left-covariance matrix of beta_{t-1} | D_{t-1}
#' @param Wt left-covariance matrix of the noise parameter of beta_{t-1} | D_{t-1}
#' @param Vt left-covariance matrix of the noise parameter of Y
#' @param nt_1 the shape parameter, or degrees of freedom, of the right-covariance matrix Sigma | D_{t-1}
#' @param Dt_1 the scale matrix of the right-covariance matrix Sigma | D_{t-1}
#' @param delta The right-variance matrix discount factor 
#' @returns The mean and covariance matrices m_t and C_t of one filtering step, the updated inverse-Wishart parameters a_t and B_t for the right-covariance matrix, plus other relevant parameters.
#' @export
FF_1step_R_sigma2R <- function(Yt, Ft, Gt, Wt, Vt, mt_1, Mt_1, nt_1, Dt_1, Rinv, delta = 1.0){
  N <- dim(Yt)[1]
  S <- dim(Yt)[2] 
  at <- Gt %*% mt_1
  At <- tcrossprod(Gt %*% Mt_1, Gt) + Wt
  qt <- Ft %*% at
  FtAt <- Ft %*% At
  Qt <- tcrossprod(FtAt, Ft) + Vt
  Qtinv <- solve(Qt)
  AttFtQtinv <- tcrossprod(At, Ft)%*% Qtinv
  Yt_qt <- Yt - qt
  mt <- at + AttFtQtinv %*% Yt_qt
  Mt <- At - AttFtQtinv %*% FtAt
  nt <- nt_1 + N * S / 2
  Dt <- Dt_1 + 1/2 * matrixcalc::matrix.trace(crossprod(Yt_qt, Qtinv) %*% Yt_qt %*% Rinv)
  
  out <- list(nt = nt, Dt = Dt, at = at, At = At,
              mt = mt, Mt = Mt, delta = delta)
  return(out)
}


# TODO FF_cpp
#' Forward Filter. Computes the FF parameters given the data at the relevant time step and the relevant parameters from the last time step.
#' @param Y The data matrix.
#' @param Ft The matrix of covariates F_t.
#' @param Gt G_t, the beta transition matrix.
#' @param m0 mean of beta_{t-1} | D_{t-1}
#' @param M0 left-covariance matrix of beta_{t-1} | D_{t-1}
#' @param Wt left-covariance matrix of the noise parameter of beta_{t-1} | D_{t-1}
#' @param Vt left-covariance matrix of the noise parameter of Y
#' @param n0 the shape parameter, or degrees of freedom, of the right-covariance matrix Sigma | D_{t-1}
#' @param D0 the scale matrix of the right-covariance matrix Sigma | D_{t-1}
#' @param delta The right-variance matrix discount factor 
#' @returns The mean and covariance matrices m_t and C_t of one filtering step, the updated inverse-Wishart parameters a_t and B_t for the right-covariance matrix, plus other relevant parameters.
#' @export
FF_sigma2R <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, n0, D0, nT, R, delta = 1.0){
  if (!is.null(dim(D0)) || length(D0) != 1) {
    stop("D0 is not a number")
  }
  
  out <- list()
  out_name <- c()
  # initialize matirces
  Ft <- F_ls
  Gt <- G_ls
  Wt <- W_ls
  Vt <- V_ls
  
  # Rinv = solve(R)
  RR <- chol(R)
  Rinv <- chol2inv(RR)
  
  for (i in 1:nT) {
    if(i %% round(0.1*nT) == 0){
      print(paste("FF:", i, "/", nT))
      print(Sys.time())
    }
    # Check if F_ls, G_ls, W_ls, V_ls are changing over time
    # If they are list, change accordingly, otherwise keep unchanged
    if(is.list(F_ls)){
      Ft <- as.matrix(F_ls[[i]])
    }
    if(is.list(G_ls)){
      Gt <- as.matrix(G_ls[[i]])
    }
    if(is.list(W_ls)){
      Wt <- as.matrix(W_ls[[i]])
    }
    if(is.list(V_ls)){
      Vt <- as.matrix(V_ls[[i]])
    }
    
    Yt <- Y[[i]]
    # Set T+1 values
    if (i == 1) {
      mt_1 = m0
      Mt_1 = M0
      nt_1 = n0
      Dt_1 = D0
    } else{
      mt_1 = one_step$mt
      Mt_1 = one_step$Mt
      nt_1 = one_step$nt
      Dt_1 = one_step$Dt
    }
    
    # perform 1 step forward
    one_step <- FF_1step_R_sigma2R(Yt = Yt, Ft = Ft, Gt = Gt, 
                                   Wt = Wt, Vt = Vt, 
                                   mt_1 = mt_1, Mt_1 = Mt_1, 
                                   nt_1 = nt_1, Dt_1 = Dt_1, 
                                   Rinv = Rinv, delta = delta)
    # one_step <- one_step[c("nt", "Dt", "at", "At", "mt", "Mt")] # extract parameters
    out[[i]] <- one_step # save results
  }
  
  out_name <- paste("T", seq(1:nT), sep = "") # accumulate output list names
  names(out) <- out_name
  out$prior <- list(m0 = m0, M0 = M0, n0 = n0, D0 = D0)
  return(out)
}

# I----
#' @export
FF_1step_R_I <- function(Yt, Ft, Gt, Wt, Vt, mt_1, Mt_1, delta = 1.0){
  N <- dim(Yt)[1]
  S <- dim(Yt)[2] 
  at <- Gt %*% mt_1
  At <- tcrossprod(Gt %*% Mt_1, Gt) + Wt
  qt <- Ft %*% at
  FtAt <- Ft %*% At
  Qt <- tcrossprod(FtAt, Ft) + Vt
  Qtinv <- solve(Qt)
  AttFtQtinv <- tcrossprod(At, Ft)%*% Qtinv
  Yt_qt <- Yt - qt
  mt <- at + AttFtQtinv %*% Yt_qt
  Mt <- At - AttFtQtinv %*% FtAt
  
  out <- list(at = at, At = At,
              mt = mt, Mt = Mt, delta = delta)
  return(out)
}


#' @export
FF_I <- function(Y, F_ls, G_ls, W_ls, V_ls, m0, M0, nT, delta = 1.0){
  out <- list()
  out_name <- c()
  # initialize matirces
  Ft <- F_ls
  Gt <- G_ls
  Wt <- W_ls
  Vt <- V_ls
  
  for (i in 1:nT) {
    if(i %% round(0.1*nT) == 0){
      print(paste("FF:", i, "/", nT))
      print(Sys.time())
    }
    # Check if F_ls, G_ls, W_ls, V_ls are changing over time
    # If they are list, change accordingly, otherwise keep unchanged
    if(is.list(F_ls)){
      Ft <- as.matrix(F_ls[[i]])
    }
    if(is.list(G_ls)){
      Gt <- as.matrix(G_ls[[i]])
    }
    if(is.list(W_ls)){
      Wt <- as.matrix(W_ls[[i]])
    }
    if(is.list(V_ls)){
      Vt <- as.matrix(V_ls[[i]])
    }
    
    Yt <- Y[[i]]
    # Set T+1 values
    if (i == 1) {
      mt_1 = m0
      Mt_1 = M0
    } else{
      mt_1 = one_step$mt
      Mt_1 = one_step$Mt
    }
    
    # perform 1 step forward
    one_step <- FF_1step_R_I(Yt = Yt, Ft = Ft, Gt = Gt, 
                                   Wt = Wt, Vt = Vt, 
                                   mt_1 = mt_1, Mt_1 = Mt_1, 
                                   delta = delta)
    # one_step <- one_step[c("nt", "Dt", "at", "At", "mt", "Mt")] # extract parameters
    out[[i]] <- one_step # save results
  }
  
  out_name <- paste("T", seq(1:nT), sep = "") # accumulate output list names
  names(out) <- out_name
  out$prior <- list(m0 = m0, M0 = M0)
  return(out)
}
