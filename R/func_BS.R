## MNIW----

#' A single step of the backwards smoothing algorithm. The parameters are taken
#'   from the forward filter case. This computes the parameters of beta_t | D_T.
#'   beta_{t+1} | D_T is assumed to be normal with mean s and covariance matrix S.
#' 
#' @param s mean of beta_{t+1} | D_T.
#' @param S left-covariance matrix of beta_{t+1} | D_T
#' @param m mean of beta_{t} | D_t, computed in the FF step
#' @param G G_t, the beta transition matrix.
#' @param C left-covariance matrix of beta_{t} | D_t, computed in the FF step
#' @param uRt upper-triangular Cholesky factor of R_t computed in the FF step.
#' @param delta right-variance matrix discount factor from the FF
#' @returns The mean and left-covariance matrix of beta_{t} | D_T
#' @export
# BS_1step_R_naiive <- function(mt, Mt, st1, St1, at1, At1, Gt1, delta = 1){
#   out <- list()
#   At1inv <- solve(At1)
#   st <- mt + Mt %*% t(Gt1) %*% At1inv %*% (st1 - at1)
#   St <- Mt - Mt %*% t(Gt1) %*% At1inv %*% (At1 - St1) %*% At1inv %*% Gt1 %*% Mt
#   
#   out <- list(st = st, St = St)
#   return(out)
# }


#' A single step of the backwards smoothing algorithm. The parameters are taken
#'   from the forward filter case. This computes the parameters of beta_t | D_T.
#'   beta_{t+1} | D_T is assumed to be normal with mean s and covariance matrix S.
#' 
#' @param s mean of beta_{t+1} | D_T.
#' @param S left-covariance matrix of beta_{t+1} | D_T
#' @param m mean of beta_{t} | D_t, computed in the FF step
#' @param G G_t, the beta transition matrix.
#' @param C left-covariance matrix of beta_{t} | D_t, computed in the FF step
#' @param uRt upper-triangular Cholesky factor of R_t computed in the FF step.
#' @param delta right-variance matrix discount factor from the FF
#' @returns The mean and left-covariance matrix of beta_{t} | D_T
#' @export
# BS_1step_R <- function(mt, Mt, st1, St1, at1, At1, Gt1, delta = 1){
#   out <- list()
#   At1inv <- solve(At1)
#   MttGt1At1inv <- tcrossprod(Mt, Gt1) %*% At1inv
#   st <- mt + MttGt1At1inv %*% (st1 - at1)
#   St <- Mt - MttGt1At1inv %*% tcrossprod((At1 - St1), MttGt1At1inv)
#   
#   out <- list(st = st, St = St)
#   return(out)
# }


# TODO MNIW_sampler_cpp and rmn_chol_cppBackward sampling
#' @return nsam samples from FFBS
#' @export
# BS_samp <- function(nsam, res_ff, F_ls, V_ls, G_ls, nT, delta = 1){
#   out <- list()
#   
#   # Firstly, sample from t = nT
#   # initialize matirces
#   if(is.list(F_ls)){
#     Ft <- as.matrix(F_ls[[nT]])
#   } else{
#     Ft <- F_ls
#   }
#   if(is.list(V_ls)){
#     Vt <- as.matrix(V_ls[[nT]])
#   } else{
#     Vt <- V_ls
#   }
#   
#   para_nt <- res_ff[[nT]]
#   st1 <- para_nt$mt
#   St1 <- check_pds(para_nt$Mt)
#   para_nt$Dt <- check_pds(para_nt$Dt)
#   post_nT <- MNIW_sampler(nsam = nsam, X = Ft, H = Vt, v = para_nt$nt, 
#                           S = para_nt$Dt, C = st1, Vb = St1)
#   out[[nT]] <- post_nT$B
#   dim_array <- dim(post_nT$B)
#   
#   Gt1 <- G_ls
#   # then backward t
#   for (i in (nT-1):0) {
#     if(i %% 10 == 0){
#       print(paste("Start BS at time", i, "/", nT))
#       print(Sys.time())
#     }
#     # Check if G_ls are changing over time
#     # If they are list, change accordingly, otherwise keep unchanged
#     if(is.list(G_ls)){
#       Gt1 <- as.matrix(G_ls[[i+1]])
#     }
#     
#     if(i == 0){
#       out$T0 <- array(dim = dim_array)
#       para_nt <- BS_1step_cpp(mt = res_ff$prior$m0, Mt = res_ff$prior$M0, 
#                               st1 = st1, St1 = St1,
#                               at1 = res_ff[[i+1]]$at, At1 = res_ff[[i+1]]$At, 
#                               Gt1 = Gt1, delta = delta)
#       # get nsam samples for different sigma 
#       # update st1, St1
#       st1 <- para_nt$st
#       St1 <- check_pds(para_nt$St)
#       for (j in 1:nsam) {
#         Sigma_j <- post_nT$sigma[,,j]
#         # out$T0[,,j] <- rmn_chol(m = para_nt$st, RM = chol(para_nt$St), RSigma = chol(Sigma_j))
#         # out$T0[,,j] <- rmatrixnormal(1, M = st1, U = St1, V = Sigma_j, checkSymmetry = F, keep = F)
#         out$T0[,,j] <- rmn_cpp(m = st1, U = St1, V = Sigma_j)
#       }
#     } else{
#       out[[i]] <- array(dim = dim_array)
#       para_nt <- BS_1step_cpp(mt = res_ff[[i]]$mt, Mt = res_ff[[i]]$Mt, 
#                               st1 = st1, St1 = St1,
#                               at1 = res_ff[[i+1]]$at, At1 = res_ff[[i+1]]$At, 
#                               Gt1 = Gt1, delta = delta)
#       # update st1, St1
#       st1 <- para_nt$st
#       St1 <- check_pds(para_nt$St)
#       # get nsam samples for different sigma 
#       for (j in 1:nsam) {
#           Sigma_j <- post_nT$sigma[,,j]
#           # out[[i]][,,j] <- rmn_chol(m = para_nt$st, RM = chol(para_nt$St), RSigma = chol(Sigma_j))
#           # out[[i]][,,j] <- rmatrixnormal(1, M = st1, U = St1, V = Sigma_j, checkSymmetry = F, keep = F)
#           out[[i]][,,j] <- rmn_cpp(m = st1, U = St1, V = Sigma_j)
#           
#       }
#       if(i == 1){
#         names(out) <- paste0("T", 1:nT)
#       }
#     }
#     
#   }
#   out$Sigma <- post_nT$sigma
#   return(out)
# }


# TODO MNIW_sampler_cpp and rmn_chol_cpp
#' Backward sampling
#' @return st, St
#' @export
BS <- function(res_ff, G_ls, nT, delta = 1){
  out <- list()
  out$st <- array(dim = c(dim(res_ff[[nT]]$mt), nT))
  out$St <- array(dim = c(dim(res_ff[[nT]]$Mt), nT))
  
  # initialize matirces
  Gt1 <- G_ls

  # Firstly, calculate t = nT
  st1 <- res_ff[[nT]]$mt
  St1 <- res_ff[[nT]]$Mt
  out$st[,,nT] <- st1
  out$St[,,nT] <- St1
  
  # then backward t
  for (i in (nT-1):1) {
    if(i %% round(0.1*nT) == 0){
      print(paste("BS:", i, "/", nT))
      print(Sys.time())
    }
    # Check if G_ls are changing over time
    # If they are list, change accordingly, otherwise keep unchanged
    if(is.list(G_ls)){
      Gt1 <- as.matrix(G_ls[[i+1]])
    }
    
    para_nt <- BS_1step_cpp(mt = res_ff[[i]]$mt, Mt = res_ff[[i]]$Mt, 
                            st1 = st1, St1 = St1,
                            at1 = res_ff[[i+1]]$at, At1 = res_ff[[i+1]]$At, 
                            Gt1 = Gt1, delta = delta)
    # update st1, St1
    st1 <- para_nt$st
    St1 <- para_nt$St
    out$st[,,i] <- st1
    out$St[,,i] <- St1
  }
  
  return(out)
}


# TODO MNIW_sampler_cpp and rmn_chol_cpp
#' Backward sampling
# BS_bigdata_R_naiive_no_extract <- function(nt_ls, Dt_ls, mt_ls, Mt_ls, at_ls, At_ls, 
#                                 F_ls, V_ls, G_ls, m0, M0, fnrow, fncol,
#                                 bnrow, bncol, path_out, nT, delta = 1){
#   nsam = 1 # force to get only one sample
#   
#   # Firstly, sample from last file, last block
#   # initialize matirces
#   if(is.list(F_ls)){
#     Ft <- as.matrix(F_ls[[nT]])
#   } else{
#     Ft <- F_ls
#   }
#   if(is.list(V_ls)){
#     Vt <- as.matrix(V_ls[[nT]])
#   } else{
#     Vt <- V_ls
#   }
#   
#   # read data
#   nt_ls_nT <- as.matrix(read.csv(nt_ls[nT]))
#   Dt_ls_nT <- as.matrix(read.csv(Dt_ls[nT]))
#   mt_ls_nT <- as.matrix(read.csv(mt_ls[nT]))
#   Mt_ls_nT <- as.matrix(read.csv(Mt_ls[nT]))
#   
#   # get parameters
#   para_nt <- list()
#   n_b <- length(nt_ls_nT)
#   para_nt$nt <- nt_ls_nT[,dim(nt_ls_nT)[2]] # get last element
#   para_nt$Dt <- Dt_ls_nT[,(dim(Dt_ls_nT)[2] - dim(Dt_ls_nT)[1] + 1):dim(Dt_ls_nT)[2]] # symmetric
#   para_nt$mt <- mt_ls_nT[,(dim(mt_ls_nT)[2] - bncol + 1):dim(mt_ls_nT)[2]]
#   para_nt$Mt <- Mt_ls_nT[,(dim(Mt_ls_nT)[2] - dim(Mt_ls_nT)[1] + 1):dim(Mt_ls_nT)[2]] # symmetric
#   
#   post_nT <- MNIW_sampler(nsam = nsam, X = Ft, H = Vt, v = para_nt$nt, 
#                           S = check_pds(para_nt$Dt), C = para_nt$mt, Vb = check_pds(para_nt$Mt))
#   
#   Sigma <- post_nT$sigma[,,1]
#   write.csv(Sigma, file = paste0(path_out, "/output_BS_Sigma.csv", 
#                                   sep = ""), row.names = F)
#   out <- post_nT$B[,,1]
#   
#   # loop backward for files
#   for (f in nT:1) {
#     if(f != nT){
#       out <- list()
#     }
#     # Check if G_ls are changing over time
#     # If they are list, change accordingly, otherwise keep unchanged
#     if(is.list(G_ls)){
#       Gt1 <- as.matrix(G_ls[[f+1]])
#     } else{
#       Gt1 <- G_ls
#     }
#     
#     # read data
#     mt_ls_f <- as.matrix(read.csv(mt_ls[f]))
#     Mt_ls_f <- as.matrix(read.csv(Mt_ls[f]))
#     at_ls_f <- as.matrix(read.csv(at_ls[f]))
#     At_ls_f <- as.matrix(read.csv(At_ls[f]))
#     if(f != nT){
#       mt_ls_f_1 <- as.matrix(read.csv(mt_ls[f + 1]))
#       Mt_ls_f_1 <- as.matrix(read.csv(Mt_ls[f + 1]))
#       at_ls_f_1 <- as.matrix(read.csv(at_ls[f + 1]))
#       At_ls_f_1 <- as.matrix(read.csv(At_ls[f + 1]))
#     }
#     
#     # loop over blocks
#     for (i in n_b:1) {
#       if(f!= nT || i != n_b){ # skip the last block since already dealed
#         # very first block should be dealt with separately
#         if(i == n_b){
#           # get parameters
#           temp_para <- list()
#           temp_para$mt <- mt_ls_f[,((i - 1) * bncol + 1):(i * bncol)]
#           temp_para$Mt <- Mt_ls_f[,((i - 1) * dim(Mt_ls_f)[1] + 1):(i * dim(Mt_ls_f)[1])] # symmetric
#           temp_para$mt_1 <- mt_ls_f_1[,1:bncol]
#           temp_para$Mt_1 <- Mt_ls_f_1[,1:dim(Mt_ls_f_1)[1]] # symmetric
#           temp_para$at_1 <- at_ls_f_1[,1:bncol]
#           temp_para$At_1 <- At_ls_f_1[,1:dim(At_ls_f_1)[1]] # symmetric
#           
#           para_nt <- BS_1step_cpp(mt = temp_para$mt, Mt = temp_para$Mt, 
#                                   st1 = temp_para$mt_1, St1 = temp_para$Mt_1,
#                                   at1 = temp_para$at_1, At1 = temp_para$At_1, 
#                                   Gt1 = Gt1, delta = delta)
#           
#           # get nsam samples for different sigma 
#           # rmn_sam <- rmn_chol(m = para_nt$st, RM = chol(para_nt$St), RSigma = chol(Sigma))
#           rmn_sam <- rmatrixnormal(1, M = para_nt$st, U = check_pds(para_nt$St), V = Sigma, checkSymmetry = F, keep = F)
#           out <- rmn_sam
#           # out$T0 <- array(dim = dim_array)
#           # para_nt <- BS_1step_cpp(mt = res_ff$prior$m0, Mt = res_ff$prior$M0, 
#           #                         st1 = res_ff[[i+1]]$mt, St1 = res_ff[[i+1]]$Mt,
#           #                         at1 = res_ff[[i+1]]$at, At1 = res_ff[[i+1]]$At, 
#           #                         Gt1 = Gt1, delta = delta)
#           # # get nsam samples for different sigma 
#           # for (j in 1:nsam) {
#           #   Sigma_j <- post_nT$sigma[,,j]
#           #   out$T0[,,j] <- rmn_chol(m = para_nt$st, RM = chol(para_nt$St), RSigma = chol(Sigma_j))
#           # }
#         } else{
#           # get parameters
#           temp_para <- list()
#           temp_para$mt <- mt_ls_f[,((i - 1) * bncol + 1):(i * bncol)]
#           temp_para$Mt <- Mt_ls_f[,((i - 1) * dim(Mt_ls_f)[1] + 1):(i * dim(Mt_ls_f)[1])] # symmetric
#           temp_para$mt_1 <- mt_ls_f[,((i) * bncol + 1):((i + 1) * bncol)]
#           temp_para$Mt_1 <- Mt_ls_f[,((i) * dim(Mt_ls_f)[1] + 1):((i + 1) * dim(Mt_ls_f)[1])] # symmetric
#           temp_para$at_1 <- at_ls_f[,((i) * bncol + 1):((i + 1) * bncol)]
#           temp_para$At_1 <- At_ls_f[,((i) * dim(At_ls_f)[1] + 1):((i + 1) * dim(At_ls_f)[1])] # symmetric
#           
#           para_nt <- BS_1step_cpp(mt = temp_para$mt, Mt = temp_para$Mt, 
#                                   st1 = temp_para$mt_1, St1 = temp_para$Mt_1,
#                                   at1 = temp_para$at_1, At1 = temp_para$At_1, 
#                                   Gt1 = Gt1, delta = delta)
#           
#           # get nsam samples for different sigma 
#           # rmn_sam <- rmn_chol(m = para_nt$st, RM = chol(para_nt$St), RSigma = chol(Sigma))
#           rmn_sam <- rmatrixnormal(1, M = para_nt$st, U = check_pds(para_nt$St), V = Sigma, checkSymmetry = F, keep = F)
#           out <- cbind(out, rmn_sam)
#         }
#       }
#     }
#     # save theta in csv
#     write.csv(out, file = paste0(path_out, "/output_BS_theta_", f, ".csv", 
#                                     sep = ""), row.names = F)
#   }
# }


# BS_bigdata_R_naiive_wrong <- function(nt_ls, Dt_ls, mt_ls, Mt_ls, at_ls, At_ls, FT, VT, G_ls,  
#                                 m0, M0, fnrow, fncol, bnrow, bncol, 
#                                 path_out, nT, delta = 1, verbose = FALSE){
#   # generate index
#   ind <- generate.grid.rowsnake(fnrow = fnrow, fncol = fncol, bnrow = bnrow, bncol = bncol)
#   n_b <- dim(ind)[1] # number of blocks
#   
#   # Firstly, sample from last file, last block
#   Ft <- as.matrix(read_big_csv_quick(FT)) # results from FF last step
#   Vt <- as.matrix(read_big_csv_quick(VT))
#   # Check if G_ls are changing over time
#   # If they are list, change accordingly, otherwise keep unchanged
#   if(length(G_ls) != 1){
#     Gt1 <- as.matrix(read_big_csv_quick(filename = G_ls[nT]))
#   } else{
#     Gt1 <- as.matrix(read_big_csv_quick(filename = G_ls))
#   }
#   
#   # read data
#   nt_ls_nT <- as.matrix(read.csv(nt_ls[nT]))
#   Dt_ls_nT <- as.matrix(read.csv(Dt_ls[nT]))
#   mt_ls_nT <- as.matrix(read.csv(mt_ls[nT]))
#   Mt_ls_nT <- as.matrix(read.csv(Mt_ls[nT]))
#   
#   # get parameters for the last block of nT
#   para_nt <- list()
#   para_nt$nt <- nt_ls_nT[, n_b] # get last element
#   para_nt$Dt <- Dt_ls_nT[,(dim(Dt_ls_nT)[2] - dim(Dt_ls_nT)[1] + 1):dim(Dt_ls_nT)[2]] # symmetric
#   # check with checks <- Dt_ls_nT[,(dim(Dt_ls_nT)[2] - bncol + 1):dim(Dt_ls_nT)[2]]
#   para_nt$mt <- mt_ls_nT[,(dim(mt_ls_nT)[2] - bncol + 1):dim(mt_ls_nT)[2]]
#   para_nt$Mt <- Mt_ls_nT[,(dim(Mt_ls_nT)[2] - dim(Mt_ls_nT)[1] + 1):dim(Mt_ls_nT)[2]] # symmetric
#   # check with checks <- Mt_ls_nT[,(dim(Mt_ls_nT)[2] - dim(Ft)[2] + 1):dim(Mt_ls_nT)[2]]
#   
#   # try, if failing, make matrix symmetric and p.d.
#   e <- simpleError("test error")
#   tryCatch(
#     expr = {post_nT <- MNIW_sampler(nsam = nsam, X = Ft, H = Vt, v = para_nt$nt, 
#                                     S = para_nt$Dt, C = para_nt$mt, Vb = para_nt$Mt)},
#     error = function(e) e, 
#     finally = {post_nT <- MNIW_sampler(nsam = nsam, X = Ft, H = Vt, v = para_nt$nt, 
#                                        S = check_pds(para_nt$Dt), C = para_nt$mt, Vb = check_pds(para_nt$Mt))})
#   
#   
#   # TODO save nsam of sigma
#   Sigma <- post_nT$sigma
#   # write.csv(Sigma, file = paste0(path_out, "/output_BS_Sigma.csv", 
#   #                                sep = ""), row.names = F)
#   save(Sigma, file = here(path_out, "output_BS_Sigma.RData"))
#   out <- post_nT$B
#   
#   # loop backward for files
#   for (f in nT:1) {
#     if(f != nT){
#       out <- list()
#       if(length(G_ls) != 1){
#         Gt1 <- as.matrix(read_big_csv_quick(filename = G_ls[f+1]))
#       }
#     }
#     
#     # read data for later extraction
#     mt_ls_f <- as.matrix(read.csv(mt_ls[f]))
#     Mt_ls_f <- as.matrix(read.csv(Mt_ls[f]))
#     at_ls_f <- as.matrix(read.csv(at_ls[f]))
#     At_ls_f <- as.matrix(read.csv(At_ls[f]))
#     if(f != nT){
#       mt_ls_f_1 <- as.matrix(read.csv(mt_ls[f + 1]))
#       Mt_ls_f_1 <- as.matrix(read.csv(Mt_ls[f + 1]))
#       at_ls_f_1 <- as.matrix(read.csv(at_ls[f + 1]))
#       At_ls_f_1 <- as.matrix(read.csv(At_ls[f + 1]))
#     }
#     
#     # loop over blocks
#     for (i in n_b:1) {
#       if(f!= nT || i != n_b){ # skip the last block since already calculate
#         # very first block should be dealt with separately
#         if(i == n_b){
#           # get parameters
#           temp_para <- list()
#           temp_para$mt <- mt_ls_f[,((i - 1) * bncol + 1):(i * bncol)]
#           temp_para$Mt <- Mt_ls_f[,((i - 1) * dim(Mt_ls_f)[1] + 1):(i * dim(Mt_ls_f)[1])] # symmetric
#           temp_para$mt_1 <- mt_ls_f_1[,1:bncol]
#           temp_para$Mt_1 <- Mt_ls_f_1[,1:dim(Mt_ls_f_1)[1]] # symmetric
#           temp_para$at_1 <- at_ls_f_1[,1:bncol]
#           temp_para$At_1 <- At_ls_f_1[,1:dim(At_ls_f_1)[1]] # symmetric
#           
#           para_nt <- BS_1step_cpp(mt = temp_para$mt, Mt = temp_para$Mt, 
#                                   st1 = temp_para$mt_1, St1 = temp_para$Mt_1,
#                                   at1 = temp_para$at_1, At1 = temp_para$At_1, 
#                                   Gt1 = Gt1, delta = delta)
#           
#           # get nsam samples for different sigma 
#           # rmn_sam <- rmn_chol(m = para_nt$st, RM = chol(para_nt$St), RSigma = chol(Sigma))
#           # rmn_sam <- rmatrixnormal(1, M = para_nt$st, U = check_pds(para_nt$St), V = Sigma, checkSymmetry = F, keep = F)
#           rmn_sam <- rmn_cpp(m = para_nt$st, U = check_pds(para_nt$St), V = Sigma)
#           
#           out <- rmn_sam
#           # out$T0 <- array(dim = dim_array)
#           # para_nt <- BS_1step_cpp(mt = res_ff$prior$m0, Mt = res_ff$prior$M0, 
#           #                         st1 = res_ff[[i+1]]$mt, St1 = res_ff[[i+1]]$Mt,
#           #                         at1 = res_ff[[i+1]]$at, At1 = res_ff[[i+1]]$At, 
#           #                         Gt1 = Gt1, delta = delta)
#           # # get nsam samples for different sigma 
#           # for (j in 1:nsam) {
#           #   Sigma_j <- post_nT$sigma[,,j]
#           #   out$T0[,,j] <- rmn_chol(m = para_nt$st, RM = chol(para_nt$St), RSigma = chol(Sigma_j))
#           # }
#         } else{
#           # get parameters: parameters are stacked by column in a wider matrix
#           temp_para <- list()
#           temp_para$mt <- mt_ls_f[,((i - 1) * bncol + 1):(i * bncol)]
#           temp_para$Mt <- Mt_ls_f[,((i - 1) * dim(Mt_ls_f)[1] + 1):(i * dim(Mt_ls_f)[1])] # symmetric
#           temp_para$mt_1 <- mt_ls_f[,((i) * bncol + 1):((i + 1) * bncol)]
#           temp_para$Mt_1 <- Mt_ls_f[,((i) * dim(Mt_ls_f)[1] + 1):((i + 1) * dim(Mt_ls_f)[1])] # symmetric
#           temp_para$at_1 <- at_ls_f[,((i) * bncol + 1):((i + 1) * bncol)]
#           temp_para$At_1 <- At_ls_f[,((i) * dim(At_ls_f)[1] + 1):((i + 1) * dim(At_ls_f)[1])] # symmetric
#           
#           para_nt <- BS_1step_cpp(mt = temp_para$mt, Mt = temp_para$Mt, 
#                                   st1 = temp_para$mt_1, St1 = temp_para$Mt_1,
#                                   at1 = temp_para$at_1, At1 = temp_para$At_1, 
#                                   Gt1 = Gt1, delta = delta)
#           
#           # get nsam samples for different sigma 
#           # rmn_sam <- rmn_chol(m = para_nt$st, RM = chol(para_nt$St), RSigma = chol(Sigma))
#           # rmn_sam <- rmatrixnormal(1, M = para_nt$st, U = check_pds(para_nt$St), V = Sigma, checkSymmetry = F, keep = F)
#           rmn_sam <- rmn_cpp(m = para_nt$st, U = check_pds(para_nt$St), V = Sigma)
#           out <- cbind(out, rmn_sam)
#         }
#       }
#     }
#     # save theta in csv
#     # write.csv(out, file = paste0(path_out, "/output_BS_theta_", f, ".csv", 
#     #                              sep = ""), row.names = F)
#     save(out, file = here(path_out, paste0("output_BS_theta_", f, ".csv", sep = "")))
#     
#     if(verbose == T){
#       print(paste0("Backward sampling", f, "/", nT))
#     }
#   }
#   
#   file_ls <- list()
#   file_ls$st_ls <- paste0(path_out, "/output_BS_st_", seq(1:nT), ".csv", sep = "")
#   file_ls$St_ls <- paste0(path_out, "/output_BS_SSt_", seq(1:nT), ".csv", sep = "")
#   
#   return(file_ls)
# }

# TODO HOW to save and read Dt; MNIW_sampler_cpp and rmn_chol_cpp
#' Backward sampling
#' @export
# BS_bigdata_R <- function(nt_ls, Dt_ls, mt_ls, Mt_ls, at_ls, At_ls, G_ls,  
#                                 m0, M0, fnrow, fncol, bnrow, bncol, 
#                                 path_out, nT, delta = 1, verbose = FALSE){
#   # generate index
#   ind <- generate.grid.rowsnake(fnrow = fnrow, fncol = fncol, bnrow = bnrow, bncol = bncol)
#   n_b <- dim(ind)[1] # number of blocks
#   
#   # Check if G_ls are changing over time
#   # If they are list, change accordingly, otherwise keep unchanged
#   if(length(G_ls) != 1){
#     Gt1 <- as.matrix(read_big_csv_quick(filename = G_ls[nT]))
#   } else{
#     Gt1 <- as.matrix(read_big_csv_quick(filename = G_ls))
#   }
#  
#   # read data
#   mt_ls_nt <- as.matrix(read.csv(mt_ls[nT]))
#   Mt_ls_nt <- as.matrix(read.csv(Mt_ls[nT]))
#   
#   # get parameters for the last block of nT
#   para_nt <- list()
#   para_nt$mt <- mt_ls_nt[,(dim(mt_ls_nt)[2] - bncol + 1):dim(mt_ls_nt)[2]]
#   para_nt$Mt <- Mt_ls_nt[,(dim(Mt_ls_nt)[2] - dim(Mt_ls_nt)[1] + 1):dim(Mt_ls_nt)[2]] # symmetric
#   st1 <- para_nt$mt
#   St1 <- para_nt$Mt
#   out_st <- st1
#   out_St <- St1
#   
#   # loop backward for files
#   for (f in nT:1) {
#     if(f != nT){
#       out_st <- matrix()
#       out_St <- matrix()
#       if(length(G_ls) != 1){
#         Gt1 <- as.matrix(read_big_csv_quick(filename = G_ls[f+1]))
#       }
#     }
#     
#     # read data for later extraction
#     mt_ls_f <- as.matrix(read.csv(mt_ls[f]))
#     Mt_ls_f <- as.matrix(read.csv(Mt_ls[f]))
#     at_ls_f <- as.matrix(read.csv(at_ls[f]))
#     At_ls_f <- as.matrix(read.csv(At_ls[f]))
#     if(f != nT){
#       mt_ls_f_1 <- as.matrix(read.csv(mt_ls[f + 1]))
#       Mt_ls_f_1 <- as.matrix(read.csv(Mt_ls[f + 1]))
#       at_ls_f_1 <- as.matrix(read.csv(at_ls[f + 1]))
#       At_ls_f_1 <- as.matrix(read.csv(At_ls[f + 1]))
#     }
#     
#     # loop over blocks
#     for (i in n_b:1) {
#       if(f!= nT || i != n_b){ # skip the last block since already calculate
#         # very first block should be dealt with separately
#         if(i == n_b){
#           # get parameters
#           temp_para <- list()
#           temp_para$mt <- mt_ls_f[,((i - 1) * bncol + 1):(i * bncol)]
#           temp_para$Mt <- Mt_ls_f[,((i - 1) * dim(Mt_ls_f)[1] + 1):(i * dim(Mt_ls_f)[1])] # symmetric
#           temp_para$mt_1 <- mt_ls_f_1[,1:bncol]
#           temp_para$Mt_1 <- Mt_ls_f_1[,1:dim(Mt_ls_f_1)[1]] # symmetric
#           temp_para$at_1 <- at_ls_f_1[,1:bncol]
#           temp_para$At_1 <- At_ls_f_1[,1:dim(At_ls_f_1)[1]] # symmetric
#           
#           para_nt <- BS_1step_cpp(mt = temp_para$mt, Mt = temp_para$Mt, 
#                                   st1 = temp_para$mt_1, St1 = temp_para$Mt_1,
#                                   at1 = temp_para$at_1, At1 = temp_para$At_1, 
#                                   Gt1 = Gt1, delta = delta)
#         } else{
#           # get parameters: parameters are stacked by column in a wider matrix
#           temp_para <- list()
#           temp_para$mt <- mt_ls_f[,((i - 1) * bncol + 1):(i * bncol)]
#           temp_para$Mt <- Mt_ls_f[,((i - 1) * dim(Mt_ls_f)[1] + 1):(i * dim(Mt_ls_f)[1])] # symmetric
#           temp_para$mt_1 <- mt_ls_f[,((i) * bncol + 1):((i + 1) * bncol)]
#           temp_para$Mt_1 <- Mt_ls_f[,((i) * dim(Mt_ls_f)[1] + 1):((i + 1) * dim(Mt_ls_f)[1])] # symmetric
#           temp_para$at_1 <- at_ls_f[,((i) * bncol + 1):((i + 1) * bncol)]
#           temp_para$At_1 <- At_ls_f[,((i) * dim(At_ls_f)[1] + 1):((i + 1) * dim(At_ls_f)[1])] # symmetric
#           
#           para_nt <- BS_1step_cpp(mt = temp_para$mt, Mt = temp_para$Mt, 
#                                   st1 = temp_para$mt_1, St1 = temp_para$Mt_1,
#                                   at1 = temp_para$at_1, At1 = temp_para$At_1, 
#                                   Gt1 = Gt1, delta = delta)
#         }
#       }
#     }
#     # save theta in csv
#     # write.csv(out, file = paste0(path_out, "/output_BS_theta_", f, ".csv", 
#     #                              sep = ""), row.names = F)
#     save(out, file = here(path_out, paste0("output_BS_theta_", f, ".csv", sep = "")))
#     
#     if(verbose == T){
#       print(paste0("Backward sampling", f, "/", nT))
#     }
#   }
#   
#   file_ls <- list()
#   file_ls$st_ls <- paste0(path_out, "/output_BS_st_", seq(1:nT), ".csv", sep = "")
#   file_ls$St_ls <- paste0(path_out, "/output_BS_SSt_", seq(1:nT), ".csv", sep = "")
#   
#   return(file_ls)
# }


# sigma2R----
# TODO MNIW_sampler_cpp and rmn_chol_cpp
#' Backward sampling
#' @export
# BS_sigma2R <- function(nsam, res_ff, F_ls, V_ls, G_ls, nT, R, delta = 1){
#   out <- list()
# 
#   # Firstly, sample from t = nT
#   # initialize matirces
#   if(is.list(F_ls)){
#     Ft <- as.matrix(F_ls[[nT]])
#   } else{
#     Ft <- F_ls
#   }
#   if(is.list(V_ls)){
#     Vt <- as.matrix(V_ls[[nT]])
#   } else{
#     Vt <- V_ls
#   }
# 
#   para_nt <- res_ff[[nT]]
#   post_nT <- MNIG_sampler(nsam = nsam, X = Ft, v = para_nt$nt,
#                           S = para_nt$Dt, C = para_nt$mt, Vb = check_pds(para_nt$Mt), R = R)
#   out[[nT]] <- post_nT$B
#   dim_array <- dim(post_nT$B)
# 
#   # then backward t
#   for (i in (nT-1):0) {
#     # print(i)
#     # print(paste("BS step", i))
#     # Check if G_ls are changing over time
#     # If they are list, change accordingly, otherwise keep unchanged
#     if(is.list(G_ls)){
#       Gt1 <- as.matrix(G_ls[[i+1]])
#     } else{
#       Gt1 <- G_ls
#     }
# 
#     if(i == 0){
#       out$T0 <- array(dim = dim_array)
#       para_nt <- BS_1step_cpp(mt = res_ff$prior$m0, Mt = res_ff$prior$M0,
#                               st1 = res_ff[[i+1]]$mt, St1 = res_ff[[i+1]]$Mt,
#                               at1 = res_ff[[i+1]]$at, At1 = res_ff[[i+1]]$At,
#                               Gt1 = Gt1, delta = delta)
#       # get nsam samples for different sigma
#       para_nt$St <- check_pds(para_nt$St)
#       for (j in 1:nsam) {
#         Sigma_j <- post_nT$sigma2[j] * R
#         # out$T0[,,j] <- rmn_chol(m = para_nt$st, RM = chol(para_nt$St), RSigma = chol(Sigma_j))
#         # out$T0[,,j] <- rmatrixnormal(1, M = para_nt$st, U = para_nt$St, V = Sigma_j, checkSymmetry = F, keep = F)
#         out$T0[,,j] <- rmn_cpp(m = para_nt$st, U = para_nt$St, V = Sigma_j)
#       }
#     } else{
#       out[[i]] <- array(dim = dim_array)
#       para_nt <- BS_1step_cpp(mt = res_ff[[i]]$mt, Mt = res_ff[[i]]$Mt,
#                               st1 = res_ff[[i+1]]$mt, St1 = res_ff[[i+1]]$Mt,
#                               at1 = res_ff[[i+1]]$at, At1 = res_ff[[i+1]]$At,
#                               Gt1 = Gt1, delta = delta)
# 
#       # get nsam samples for different sigma
#       para_nt$St <- check_pds(para_nt$St)
#       for (j in 1:nsam) {
#         Sigma_j <- post_nT$sigma2[j] * R
#         # out[[i]][,,j] <- rmn_chol(m = para_nt$st, RM = chol(para_nt$St), RSigma = chol(Sigma_j))
#         # out[[i]][,,j] <- rmatrixnormal(1, M = para_nt$st, U = para_nt$St, V = Sigma_j, checkSymmetry = F, keep = F)
#         out[[i]][,,j] <- rmn_cpp(m = para_nt$st, U = para_nt$St, V = Sigma_j)
#       }
#       if(i == 1){
#         names(out) <- paste0("T", 1:nT)
#       }
#     }
# 
#     if(i %% 10 == 0){
#       print(paste("Finish BS at time", i))
#       print(Sys.time())
#     }
#   }
#   out$Sigma <- post_nT$sigma
#   return(out)
# }