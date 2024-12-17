#' Generates the model comparison statistic of Gelfand and Ghosh (1998) for the DLM. 
#' @param object The object for the independent replicate model comparison to analyze. The default treats object as a formula to be applied to the dataset over which the original model was fitted. The function also supports the big_data_set object.
#' @param data The data set that was fitted by the FFBS.
#' @param 

GPD <- function(object, ...) {
    UseMethod("GPD")
}

GPD.big_data_set <- function(object, mean_files, left_cov_files,
					  nt_files = NULL, Dt_files = NULL,
					  dt_files = NULL, R = NULL,
					  Sigma_mode = "IW",
					  Vt = function(chunk) return(diag(nrow(chunk))),
					  NITER=NULL,
			                  save_dir = NULL,
					  out.format = "csv",
					  fst.head = "V") {
    #preprocess across different Sigma_mode's here
    if((Sigma_mode == "IW") & ( is.null(nt_files) | is.null(Dt_files))) stop("nt_files and Dt_files must be provided for Sigma mode IW (inverse-Wishart).")
    if ((Sigma_mode == "IG") & (is.null(nt_files) | is.null(dt_files) | is.null(R))) stop("nt_files, dt_files, and R must be provided for Sigma mode IG (inverse-gamma)!")
    if (!(out.format %in% c("csv", "fst"))) stop("The MCIR only supports the 'csv' or 'fst' data formats!")
    
    # Call the specific model comparison function here.
    if (Sigma_mode == "id") {
        return(GPD.big_data_set.id(object, mean_files, left_cov_files, Vt = Vt,
						NITER=NITER,
						save_dir = save_dir,
						out.format = out.format, fst.head = fst.head))

    } else if (Sigma_mode == "IW") {
        return(GPD.big_data_set.IW(object, mean_files, left_cov_files,
			         nt_files, Dt_files, Vt = Vt,
				 NITER=NITER,
			         save_dir = save_dir,
				 out.format = out.format, fst.head = fst.head))
    
    } else if (Sigma_mode == "IG") {
        return(GPD.big_data_set.IG(object, mean_files, left_cov_files,
				    nt_files, dt_files, R, Vt = Vt,
				    NITER=NITER,
				    save_dir = save_dir,
				    out.format = out.format, fst.head = fst.head))
    } else {
        stop(sprintf("Sigma mode %s not recognized!", Sigma_mode))
    }
}


GPD.big_data_set.id <- function(big_data_set, mean_files, left_cov_files, 
					     Vt = function(chunk) diag(nrow(chunk)), NITER=NULL,
			                     save_dir = NULL,
					     out.format = "csv", fst.head = "V") {
    T <- length(mean_files)

    G <- 0
    P <- 0

    # record per-episode GPD
    sum_K <- big_data_set$get("sum_K")
    GPD_perep <- data.frame(tk_ix = seq(1, sum_K),
			    G = rep(-1, sum_K),
                            P = rep(-1, sum_K),
                            D = rep(-1, sum_K))

    tk_ix <- 1
    for (t in 1:T) {
        Y_file <- big_data_set$get(paste0("Y",t))
        file_name_Y <- Y_file$file_name

        F_file <- big_data_set$get(paste0("F",t))
        file_name_F <- F_file$file_name

	Ynrow <- Y_file$nrow
        Yncol <- Y_file$ncol
	Ybnrow <- Y_file$nrowblock
	Ybncol <- Y_file$ncolblock

	Fnrow <- F_file$nrow
        Fncol <- F_file$ncol
	Fbnrow <- F_file$nrowblock
	Fbncol <- F_file$ncolblock

        st_file <- if (out.format == "csv") {
	       	       read.csv(mean_files[t])
	           } else if (out.format == "fst") {
	       	       read.fst(mean_files[t])
	           }
	St_file <- if (out.format == "csv") {
		       read.csv(left_cov_files[t])
	           } else if (out.format == "fst") {
		       read.fst(left_cov_files[t])
	           }

	K <- Y_file$K

	for (k in 1:K) {
            # read in Yk
	    ud_Y <- unlist(Y_file$grid[Y_file$grid$block == k, c("U", "D")])
	    lr_Y <- unlist(Y_file$grid[Y_file$grid$block == k, c("L", "R")])
	    Y_tk <- if (out.format == "csv") {
		        as.matrix(read_big_csv(file_name_Y, rows = ud_Y, cols = lr_Y,
				  	       has_header=Y_file$header,
					       ignore_header=TRUE))
	            } else if (out.format == "fst") {
		        as.matrix(read_big_fst(file_name_Y, rows = ud_Y,
					       cols = paste0(fst.head, seq(lr_Y[1],lr_Y[2]))))
	            }
            
	    # read in Fk
	    ud_F <- unlist(F_file$grid[F_file$grid$block == k, c("U", "D")])
	    lr_F <- unlist(F_file$grid[F_file$grid$block == k, c("L", "R")])
	    F_tk <- if (out.format == "csv") {
		        as.matrix(read_big_csv(file_name_F, rows = ud_F, cols = lr_F,
					       has_header=F_file$header,
					       ignore_header=TRUE))
	            } else if (out.format == "fst") {
		        as.matrix(read_big_fst(file_name_F, rows = ud_F, 
					       cols = paste0(fst.head, seq(lr_F[1], lr_F[2]))))
	            }
	    
	    # read in sk and Sk
	    s_tk <- matrix(unlist(st_file[k,]), nrow=Fbncol, ncol=Ybncol)
	    S_tk <- tri.to.sym(unlist(St_file[k,]), Fbncol, lower=FALSE)
            
            # compute WAIC terms
	    Y_diff <- if (is.null(NITER)) {
		          Y_tk - F_tk %*% s_tk
	              } else {
                          Y_tk - G_sample_Ymean(F_tk %*% s_tk, Vt(F_tk), "id", NITER=NITER)
	              }
	    
	    # calculate the terms G and P.
            G_tk <- sum(Y_diff * Y_diff)
	    P_tk <- Ybncol * sum(diag(S_tk))

            GPD_perep[tk_ix, "G"] <- G_tk
            GPD_perep[tk_ix, "P"] <- P_tk

	    G <- G + G_tk
	    P <- P + P_tk

	    tk_ix <- tk_ix + 1
       	}
    }

    if (!is.null(save_dir)) {
        GPD_perep$D <- GPD_perep$G + GPD_perep$P
        write.csv(GPD_perep, file=file.path(save_dir, "GPD_perep.csv"), row.names=FALSE)
    }

    D <- G + P

    return(list(D = D, G = G, P = P))
}

GPD.big_data_set.IW <- function(big_data_set, mean_files, left_cov_files,
			         nt_files, Dt_files, Vt = function(chunk) diag(nrow(chunk)), NITER=NULL,
			         save_dir = NULL,
			         out.format = "csv", fst.head = "V") {

    update_ntDt <- (!is.numeric(nt_files)) & (length(nt_files) > 1)
    if (!update_ntDt) {
        # load them naively if they do not need to be updated.
        nt <- nt_files
        Dt <- Dt_files
    }

    T <- length(mean_files)

    G <- 0
    P <- 0

    # record per-episode GPD
    sum_K <- big_data_set$get("sum_K")
    GPD_perep <- data.frame(tk_ix = seq(1, sum_K),
			    G = rep(-1, sum_K),
                            P = rep(-1, sum_K),
                            D = rep(-1, sum_K))

    tk_ix <- 1
    for (t in 1:T) {
        Y_file <- big_data_set$get(paste0("Y",t))
        file_name_Y <- Y_file$file_name

        F_file <- big_data_set$get(paste0("F",t))
        file_name_F <- F_file$file_name

	Ynrow <- Y_file$nrow
        Yncol <- Y_file$ncol
	Ybnrow <- Y_file$nrowblock
	Ybncol <- Y_file$ncolblock

	Fnrow <- F_file$nrow
        Fncol <- F_file$ncol
	Fbnrow <- F_file$nrowblock
	Fbncol <- F_file$ncolblock

        st_file <- if (out.format == "csv") {
	       	       read.csv(mean_files[t])
	           } else if (out.format == "fst") {
	       	       read.fst(mean_files[t])
	           }
	St_file <- if (out.format == "csv") {
		       read.csv(left_cov_files[t])
	           } else if (out.format == "fst") {
		       read.fst(left_cov_files[t])
	           }

	# load the nt and Dt files here
	if (update_ntDt) {
            nt_list <- unlist(read.csv(nt_files[t], header=FALSE))
	    names(nt_list) <- NULL

	    Dt_list <- if (out.format == "csv") {
		           read.csv(Dt_files[t])
	               } else if (out.format == "fst") {
                           read.fst(Dt_files[t])
	               }
	}

	K <- Y_file$K

	for (k in 1:K) {
            # read in Yk
	    ud_Y <- unlist(Y_file$grid[Y_file$grid$block == k, c("U", "D")])
	    lr_Y <- unlist(Y_file$grid[Y_file$grid$block == k, c("L", "R")])
	    Y_tk <- if (out.format == "csv") {
		        as.matrix(read_big_csv(file_name_Y, rows = ud_Y, cols = lr_Y,
				  	       has_header=Y_file$header,
					       ignore_header=TRUE))
	            } else if (out.format == "fst") {
		        as.matrix(read_big_fst(file_name_Y, rows = ud_Y,
					       cols = paste0(fst.head, seq(lr_Y[1],lr_Y[2]))))
	            }
            
	    # read in Fk
	    ud_F <- unlist(F_file$grid[F_file$grid$block == k, c("U", "D")])
	    lr_F <- unlist(F_file$grid[F_file$grid$block == k, c("L", "R")])
	    F_tk <- if (out.format == "csv") {
		        as.matrix(read_big_csv(file_name_F, rows = ud_F, cols = lr_F,
					       has_header=F_file$header,
					       ignore_header=TRUE))
	            } else if (out.format == "fst") {
		        as.matrix(read_big_fst(file_name_F, rows = ud_F, 
					       cols = paste0(fst.head, seq(lr_F[1], lr_F[2]))))
	            }
	    
	    # read in sk and Sk
	    s_tk <- matrix(unlist(st_file[k,]), nrow=Fbncol, ncol=Ybncol)
	    S_tk <- tri.to.sym(unlist(St_file[k,]), Fbncol, lower=FALSE)

            if (update_ntDt) {
                nt <- nt_list[k]
	        Dt <- tri.to.sym(unlist(Dt_list[k,]), Ybncol, lower=FALSE)
	    }	   

	    Y_diff <- if (is.null(NITER)) {
		          Y_tk - F_tk %*% s_tk
	              } else {
                          Y_tk - G_sample_Ymean(F_tk %*% s_tk, Vt(F_tk), "IW", nt = nt, Dt = Dt, 
			                    NITER=NITER)
	              }

	    # calculate the terms G and P.
            G_tk <- sum(Y_diff * Y_diff)
	    P_tk <- sum(diag(S_tk)) * sum(diag(Dt)) /(nt - 2)

            GPD_perep[tk_ix, "G"] <- G_tk
            GPD_perep[tk_ix, "P"] <- P_tk

	    G <- G + G_tk

	    #TODO: If NITER is not null, randomly sample Yreps and take the trace of the variance of its vectorization.
	    P <- P + P_tk

	    tk_ix <- tk_ix + 1
       	}
    }

    if (!is.null(save_dir)) {
        GPD_perep$D <- GPD_perep$G + GPD_perep$P
        write.csv(GPD_perep, file=file.path(save_dir, "GPD_perep.csv"), row.names=FALSE)
    }

    D <- G + P

    return(list(D = D, G = G, P = P))
}

GPD.big_data_set.IG <- function(big_data_set, mean_files, left_cov_files,
				 nt_files, dt_files, R,
			         Vt = function(chunk) diag(nrow(chunk)), NITER=NULL,
			         save_dir = NULL,
			         out.format = "csv", fst.head = "V") {

    update_ntdt <- (!is.numeric(nt_files)) & (length(nt_files) > 1)
    if (!update_ntdt) {
        # load them naively if they do not need to be updated.
        nt <- nt_files
        dt <- dt_files
    }

    T <- length(mean_files)

    G <- 0
    P <- 0

    # record per-episode GPD
    sum_K <- big_data_set$get("sum_K")
    GPD_perep <- data.frame(tk_ix = seq(1, sum_K),
			    G = rep(-1, sum_K),
                            P = rep(-1, sum_K),
                            D = rep(-1, sum_K))

    tk_ix <- 1

    for (t in 1:T) {
        Y_file <- big_data_set$get(paste0("Y",t))
        file_name_Y <- Y_file$file_name

        F_file <- big_data_set$get(paste0("F",t))
        file_name_F <- F_file$file_name

	Ynrow <- Y_file$nrow
        Yncol <- Y_file$ncol
	Ybnrow <- Y_file$nrowblock
	Ybncol <- Y_file$ncolblock

	Fnrow <- F_file$nrow
        Fncol <- F_file$ncol
	Fbnrow <- F_file$nrowblock
	Fbncol <- F_file$ncolblock

        st_file <- if (out.format == "csv") {
	       	       read.csv(mean_files[t])
	           } else if (out.format == "fst") {
	       	       read.fst(mean_files[t])
	           }
	St_file <- if (out.format == "csv") {
		       read.csv(left_cov_files[t])
	           } else if (out.format == "fst") {
		       read.fst(left_cov_files[t])
	           }

        if (update_ntdt) {
            nt_list <- unlist(read.csv(nt_files[t], header=FALSE))
	    names(nt_list) <- NULL
            dt_list <- unlist(read.csv(dt_files[t], header=FALSE))
	    names(dt_list) <- NULL
	}

	K <- Y_file$K

	for (k in 1:K) {
            # read in Yk
	    ud_Y <- unlist(Y_file$grid[Y_file$grid$block == k, c("U", "D")])
	    lr_Y <- unlist(Y_file$grid[Y_file$grid$block == k, c("L", "R")])
	    Y_tk <- if (out.format == "csv") {
		        as.matrix(read_big_csv(file_name_Y, rows = ud_Y, cols = lr_Y,
				  	       has_header=Y_file$header,
					       ignore_header=TRUE))
	            } else if (out.format == "fst") {
		        as.matrix(read_big_fst(file_name_Y, rows = ud_Y,
					       cols = paste0(fst.head, seq(lr_Y[1],lr_Y[2]))))
	            }
            
	    # read in Fk
	    ud_F <- unlist(F_file$grid[F_file$grid$block == k, c("U", "D")])
	    lr_F <- unlist(F_file$grid[F_file$grid$block == k, c("L", "R")])
	    F_tk <- if (out.format == "csv") {
		        as.matrix(read_big_csv(file_name_F, rows = ud_F, cols = lr_F,
					       has_header=F_file$header,
					       ignore_header=TRUE))
	            } else if (out.format == "fst") {
		        as.matrix(read_big_fst(file_name_F, rows = ud_F, 
					       cols = paste0(fst.head, seq(lr_F[1], lr_F[2]))))
	            }
	    
	    # read in sk and Sk
	    s_tk <- matrix(unlist(st_file[k,]), nrow=Fbncol, ncol=Ybncol)
	    S_tk <- tri.to.sym(unlist(St_file[k,]), Fbncol, lower=FALSE)
            
	    if (update_ntdt) {
                nt <- nt_list[k]
	        dt <- dt_list[k]
	    }

	    Y_diff <- if (is.null(NITER)) {
		          Y_tk - F_tk %*% s_tk
	              } else {
                          Y_tk - G_sample_Ymean(F_tk %*% s_tk, Vt(F_tk), "IG", nt = nt, dt = dt, R = R,
			                    NITER=NITER)
	              }

	    # calculate the terms G and P.
            G_tk <- sum(Y_diff * Y_diff)
	    P_tk <- sum(diag(S_tk)) * dt/(nt - 1)

            GPD_perep[tk_ix, "G"] <- G_tk
            GPD_perep[tk_ix, "P"] <- P_tk

	    G <- G + G_tk
	    #TODO: If NITER is not null, randomly sample Yreps and take the trace of the variance of its vectorization.
	    P <- P + P_tk

	    tk_ix <- tk_ix + 1
       	}
    }

    if (!is.null(save_dir)) {
        GPD_perep$P <- GPD_perep$P * sum(diag(R))
        GPD_perep$D <- GPD_perep$G + GPD_perep$P
        write.csv(GPD_perep, file=file.path(save_dir, "GPD_perep.csv"), row.names=FALSE)
    }

    P <- P * sum(diag(R))

    D <- G + P
    return(list(D = D, G = G, P = P))
}

# Sample the underlying mu_rep for Y_t NITER times. Takes in the ymean = Ft %*% st, Vt, and proceed.
G_sample_Ymean <- function(ymean, Vt, Sigma_mode, nt = NULL, Dt = NULL, dt = NULL, R = NULL,
			  NITER = 1000) {
    Theta_samples <- if (Sigma_mode == "id") {
	rmn(NITER, ymean, Vt, diag(ncol(ymean)))
    } else if (Sigma_mode == "IW") {
        rmatrixt(NITER, ymean, Vt, nt, Dt)
    } else if (Sigma_mode == "IG") {
        rmatrixtscS(NITER, ymean, Vt, nt, dt, R)
    }
    return(apply(Theta_samples, c(1,2), mean))        
}


