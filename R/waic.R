library(loo)
library(fst)

#' Computes the Watanabe-Akaike Information Criteria (Widely-Applicable Information Criteria) for the MNIW ensemble of objects.
#' @param object The object for the waic to analyze. The default treats object as a formula to be applied to the dataset over which the original model was fitted. The function also supports the big_data_set object.
#' @param data The data set that was fitted by the FFBS.
#' @param 


waic <- function(object, ...) {
    UseMethod("waic")
}

#NOTE: 3 cases each: id (Identity), IW (Inverse-Wishart), and IG (inverse-Gamma with a scale matrix).

##TODO: The typical MNIW FFBS case. waic.default here is the y_formula.
#waic.default <- function(object, x_formula, data, means, left_covariances,
#			      nt_list = NULL, Dt_list = NULL, Vt = function(chunk) diag(nrow(chunk)),
#			      tvar = "epoch", Sigma_mode = "IW", dt_list = NULL, R = NULL,
#			      NITER = 1000, verbose = FALSE) {
#    if ((Sigma_mode == "IW") & (is.null(df_list) | is.null(nt_list) | is.null(Dt_list))) stop("df_list, nt_list, and Dt_list must be provided for Sigma mode IW (inverse-Wishart).")
#    if ((Sigma_mode == "IG") & (is.null(df_list) | is.null(nt_list) | is.null(dt_list) | is.null(R))) stop("df_list, nt_list, dt_list, and R must be provided for Sigma mode IG (inverse-gamma)!")
#   
#    data_file_list <- unlist(data_set$keys())
#    data_file_list <- str_subset(data_file_list, "Y")
#    T <- length(data_file_list)
#
#    p <- sqrt(ncol(left_covariances))
#    S <- ncol(means)/p
#    m.cols <- setdiff(colnames(means), c(tvar))
#    M.cols <- setdiff(colnames(left_covariances), c(tvar))
#
#    form_x_t <- update(formula_x, paste("~. +",tvar))
#    mfX_t <- stats::model.frame(form_x_t, data)
#    mtX_t <- attr(x = mfX_t, which = "terms")
#
#    X_aug <- as.data.frame(model.matrix(object = mtX_t, data = mfX_t)) %>%
#            rename(epoch = all_of(tvar))
#    xnames <- colnames(X_aug)
#    xnames <- xnames[xnames != "epoch"]
#
#    #NOTE: object here refers to the y_formula.
#    if (!("0" %in% all.vars(object))) object <- update(object, "~. + 0")
#
#    #TODO: transform data the same way.
#    form_y_t <- update(formula_y, paste("~. +",tvar))
#    mfY_t <- stats::model.frame(form_y_t, data)
#    mtY_t <- attr(x = mfY_t, which = "terms")
#    Y_aug <- as.data.frame(model.matrix(object = mtY_t, data = mfY_t)) %>%
#            rename(epoch = all_of(tvar))
#    ynames <- colnames(Y_aug)
#    ynames <- ynames[ynames != "epoch"]
#
#    #TODO: update correctly
#    update_ntDt <- length(nt_list) > 1
#    if (!update_ntDt) {
#        nt <- nt_list[1]
#        dt <- dt_list[1]
#        
#    }
#   
#    WAIC_each <- rep(Inf, T)
#    WAIC <- 0 
#    for (t in 1:T) {
#        if (update_df) {
#            nt <- nt_list[t]
#	    dt <- dt_list[t]
#	}
#        if (verbose) print(sprintf("epoch: %d", t))
#	st <- matrix(unlist(means[means[[tvar]] == t, m.cols]), nrow=S, ncol=p)
#	St <- matrix(unlist(left_covariances[left_covariances[[tvar]] == t, M.cols]), nrow=p)
#
#	# compute Ft
#	N <- sum(data[[tvar]] == t)
#
#	Ft <- as.matrix(X_aug %>% filter(epoch == t) %>% select(-c(epoch)))
#	Yt <- Y_aug %>% filter(epoch == t) %>% select(-c(epoch))
#
#        #TODO: Split it based on the Sigma_mode.
#	if (Sigma_mode == "id") {
#
#	} else if (Sigma_mode == "IW") {
#
#	} else if (Sigma_mode == "IG") {
#
#	}
#    }
#
#    return(WAIC)
#}


#TODO: add direct_sample option in a future update. Not high priority.

#' Computes the Watanabe-Akaike Information Criterion (WAIC) for the FFBS transfer learning setup. 
#' @param object The big_data_set object that wraps around the filenames containing Y_t. See big_data_set.R for further information on this object. 
#' @param mean_files The smoothing coefficients generated at the end of the smoothing step of the FFBS. 
#' @param left_cov_files The left covariance moments generated at the end of the smoothing step of the FFBS.
#' @param nt_files 
#' @param Dt_files 
#' @param Vt 
#' @param Sigma_mode Set "IW" for the inverse-Wishart covariance structure, "IG" for the inverse Gamma covariance structure, or "id" for the identity covariance structure. Defaults to "IW".
#' @param dt_files 
#' @param R 
#' @param NITER
#' @param verbose 
#' @param save_dir 
#' @param out.format 
#' @param fst.head 
#' @returns A dataframe containing the WAIC, lppd, analytic lppd ("lppd_analytic"), and p_WAIC values, along with their sample standard errors across all episodes.

waic.big_data_set <- function(object, mean_files, left_cov_files,
			      nt_files = NULL, Dt_files = NULL, Vt = function(chunk) diag(nrow(chunk)),
			      Sigma_mode = "IW",dt_files = NULL, R = NULL,
			      NITER = 1000, verbose = FALSE,
			      save_dir = NULL,
			      out.format = "csv",
			      fst.head = "V") {
    if((Sigma_mode == "IW") & ( is.null(nt_files) | is.null(Dt_files))) stop("nt_files and Dt_files must be provided for Sigma mode IW (inverse-Wishart).")
    if ((Sigma_mode == "IG") & (is.null(nt_files) | is.null(dt_files) | is.null(R))) stop("nt_files, dt_files, and R must be provided for Sigma mode IG (inverse-gamma)!")
    
    # Call the specific WAIC function here.
    if (Sigma_mode == "id") {
        return(waic.big_data_set.id(object, mean_files, left_cov_files, Vt = Vt,
			      NITER = NITER, verbose = verbose, save_dir = save_dir,
			      out.format = out.format, fst.head = fst.head))

    } else if (Sigma_mode == "IW") {
        return(waic.big_data_set.IW(object, mean_files, left_cov_files,
			         nt_files, Dt_files, Vt = Vt,
			         NITER = NITER, verbose = verbose, save_dir = save_dir,
				 out.format = out.format, fst.head = fst.head))
    
    } else if (Sigma_mode == "IG") {
        return(waic.big_data_set.IG(object, mean_files, left_cov_files,
				    nt_files, dt_files, R,
			            Vt = Vt,
			            NITER = NITER, verbose = verbose, save_dir = save_dir,
				    out.format = out.format, fst.head = fst.head))
    } else {
        stop(sprintf("Sigma mode %s not recognized!", Sigma_mode))
    }
}


waic.big_data_set.id <- function(big_data_set, mean_files, left_cov_files, Vt = function(chunk) diag(nrow(chunk)),
			      NITER = 1000, verbose = FALSE, save_dir = NULL, out.format = "csv", fst.head = "V") {
    T <- length(mean_files)

    loglikmat <- matrix(rep(-Inf, NITER * big_data_set$get("sum_K")), nrow=NITER)

    WAIC_each <- rep(Inf, big_data_set$get("sum_K"))
    lppd_each <- rep(-Inf, big_data_set$get("sum_K"))
    p_waic2_each <- rep(Inf, big_data_set$get("sum_K"))

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
            if (verbose) cat(sprintf("Running time %d and block %d...", t, k))
            # read in Yk
	    ud_Y <- unlist(Y_file$grid[Y_file$grid$block == k, c("U", "D")])
	    lr_Y <- unlist(Y_file$grid[Y_file$grid$block == k, c("L", "R")])

	    Y_tk <- if (out.format == "csv") {
		        as.matrix(read_big_csv(file_name_Y, rows = ud_Y, cols = lr_Y,
		  	 		       has_header=Y_file$header,
					       ignore_header=TRUE))
	            } else if (out.format == "fst") {
                        as.matrix(read_big_fst(file_name_Y, rows = ud_Y, 
					       cols = paste0(fst.head, seq(lr_Y[1], lr_Y[2]))))
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
	    s_tk <- matrix(unlist(st_file[k+1,]), nrow=Fbncol, ncol=Ybncol)
	    S_tk <- tri.to.sym(unlist(St_file[k+1,]), Fbncol, lower=FALSE)
            
            # compute WAIC terms
	    #TODO: Save Theta_tki in an Rdata object if save_dir is specified. (i.e. as waic_rmatrixnorm_samples_t%d_k%d.R)
            Theta_tki <- rmn(NITER, s_tk, S_tk, diag(ncol(s_tk)))
            for (i in 1:NITER) {
                loglikmat[i,tk_ix] <- dmn(Y_tk, F_tk %*% Theta_tki[,,i],
	                                  Vt(F_tk), diag(ncol(Y_tk)), log=TRUE)
	    }

            lppd <- dmn(Y_tk, mu = F_tk %*% s_tk, 
			V = F_tk %*% S_tk %*% t(F_tk) + Vt(F_tk), 
			S = diag(Ybncol), log=TRUE)

#            p_waic2 <- p_waic2_id_1t(NITER, s_tk, S_tk, Y_tk, F_tk, Vt(F_tk))

	    lppd_each[tk_ix] <- lppd
#	    p_waic2_each[tk_ix] <- p_waic2
#	    WAIC_each[tk_ix] <- lppd - p_waic2
	    #if (!is.finite(WAIC_each_t[k])) break

	    tk_ix <- tk_ix + 1
	    if (verbose) cat(sprintf("Time %d episode %d completed.\n", t, k))
       	}
    }

#    WAIC <- -2 *sum(WAIC_each)
#    SE <- sqrt(var(WAIC_each) * T) * 2


    lppd_dfz <- sum(lppd_each) 
    lppd_emp_each <- apply(loglikmat, 2, logSumExp) - log(NITER)
    lppd_dfz_emp <- sum(lppd_emp_each)
    p_waic2_each <- apply(loglikmat, 2, var)
    p_waic2_dfz <- sum(p_waic2_each)

    if (!is.null(save_dir)) {
        waic_df_perep <- data.frame(WAIC = -2 * (lppd_emp_each - p_waic2_each),
                                    SE = 2 * sqrt(var(lppd_emp_each - p_waic2_each)),
                                    p_WAIC = p_waic2_each,
                                    p_WAIC_SE = sqrt(var(p_waic2_each)),
                                    lppd = lppd_emp_each,
                                    lppd_SE = sqrt(var(lppd_emp_each)),
                                    lppd_analytic = lppd_each,
                                    lppd_analytic_SE = sqrt(var(lppd_each)))
        #T save log-likelihood matrix in 
        saveRDS(loglikmat, file=file.path(save_dir, "loglikmat.rds"))
        write.csv(waic_df_perep, file=file.path(save_dir, "WAIC_perep.csv"), row.names=FALSE)
    }

    loo_obj <- loo::waic(loglikmat)

    return(list(WAIC = loo_obj$estimates["waic", "Estimate"], 
		SE = loo_obj$estimates["waic", "SE"],
		p_WAIC = loo_obj$estimates["p_waic", "Estimate"],
		p_WAIC_SE = loo_obj$estimates["p_waic", "SE"],
		lppd = loo_obj$estimates["elpd_waic", "Estimate"] + loo_obj$estimates["p_waic", "Estimate"],
		lppd_SE = sqrt(var(lppd_emp_each) * T),
		lppd_analytic = lppd_dfz,
		lppd_analytic_SE = sqrt(var(lppd_each) * T)
		))
}

waic.big_data_set.IW <- function(big_data_set, mean_files, left_cov_files,
			         nt_files, Dt_files, Vt = function(chunk) diag(nrow(chunk)),
			         NITER = 1000, verbose = FALSE, save_dir = NULL,
				 out.format = "csv", fst.head = "V") {

    update_ntDt <- (!is.numeric(nt_files)) & (length(nt_files) > 1)
    if (!update_ntDt) {
        # load them naively if they do not need to be updated.
        nt <- nt_files
        Dt <- Dt_files
    }

    T <- length(mean_files)

    loglikmat <- matrix(rep(-Inf, NITER * big_data_set$get("sum_K")), nrow=NITER)

    lppd_each <- rep(-Inf, big_data_set$get("sum_K"))
    p_waic2_each <- rep(Inf, big_data_set$get("sum_K"))
    WAIC_each <- rep(Inf, big_data_set$get("sum_K"))

    tk_ix <- 1
#    WAIC <- 0
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
#        WAIC_each_t <- c(Inf, K)

	for (k in 1:K) {
            if (verbose) cat(sprintf("Running time %d and block %d...", t, k))
            # read in Yk
	    ud_Y <- unlist(Y_file$grid[Y_file$grid$block == k, c("U", "D")])
	    lr_Y <- unlist(Y_file$grid[Y_file$grid$block == k, c("L", "R")])
	    Y_tk <- if (out.format == "csv") {
		        as.matrix(read_big_csv(file_name_Y, rows = ud_Y, cols = lr_Y,
		  	 		       has_header=Y_file$header,
					       ignore_header=TRUE))
	            } else if (out.format == "fst") {
                        as.matrix(read_big_fst(file_name_Y, rows = ud_Y, 
					       cols = paste0(fst.head, seq(lr_Y[1], lr_Y[2]))))
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
	    s_tk <- matrix(unlist(st_file[k+1,]), nrow=Fbncol, ncol=Ybncol)
	    S_tk <- tri.to.sym(unlist(St_file[k+1,]), Fbncol, lower=FALSE)

            if (update_ntDt) {
                nt <- nt_list[k+1]
	        Dt <- tri.to.sym(unlist(Dt_list[k+1,]), Ybncol, lower=FALSE)
	    }	    

            Theta_tki <- rmatrixt(NITER, s_tk, S_tk, nt, Dt)
            for (i in 1:NITER) {
                loglikmat[i, tk_ix] <- dmatrixt(Y_tk, F_tk %*% Theta_tki[,,i], Vt(F_tk), 
        			                nt, Dt, log=TRUE)
            } 

            lppd <- dmatrixt(Y_tk, F_tk %*% s_tk, F_tk %*% S_tk %*% t(F_tk) + Vt(F_tk), 
			     nt, Dt, log=TRUE)

#            p_waic2 <- p_waic2_IW_1t(NITER, s_tk, S_tk, nt, Dt, Y_tk, F_tk, Vt(F_tk))

#	    WAIC_each_t[k] <- lppd - p_waic2
	    #if (!is.finite(WAIC_each_t[k])) break
	    lppd_each[tk_ix] <- lppd
#	    p_waic2_each[tk_ix] <- p_waic2
#	    WAIC_each[tk_ix] <- lppd - p_waic2
	    
	    tk_ix <- tk_ix + 1
	    if (verbose) cat(sprintf("Time %d episode %d completed.\n", t, k))
       	}
#	WAIC_each <- c(WAIC_each, WAIC_each_t)
    }

#    WAIC <- -2 *sum(WAIC_each)
#    SE <- sqrt(var(WAIC_each) * T) * 2

    loo_obj <- loo::waic(loglikmat)

    lppd_dfz <- sum(lppd_each)
    lppd_emp_each <- apply(loglikmat, 2, logSumExp) - log(NITER)
    lppd_dfz_emp <- sum(lppd_emp_each)
    p_waic2_each <- apply(loglikmat, 2, var)
    p_waic2_dfz <- sum(p_waic2_each)

    if (!is.null(save_dir)) {
        waic_df_perep <- data.frame(WAIC = -2 * (lppd_emp_each - p_waic2_each),
                                    SE = 2 * sqrt(var(lppd_emp_each - p_waic2_each)),
                                    p_WAIC = p_waic2_each,
                                    p_WAIC_SE = sqrt(var(p_waic2_each)),
                                    lppd = lppd_emp_each,
                                    lppd_SE = sqrt(var(lppd_emp_each)),
                                    lppd_analytic = lppd_each,
                                    lppd_analytic_SE = sqrt(var(lppd_each)))
        #T save log-likelihood matrix in 
        saveRDS(loglikmat, file=file.path(save_dir, "loglikmat.rds"))
        write.csv(waic_df_perep, file=file.path(save_dir, "WAIC_perep.csv"), row.names=FALSE)
    }

    loo_obj <- loo::waic(loglikmat)

    return(list(WAIC = loo_obj$estimates["waic", "Estimate"], 
		SE = loo_obj$estimates["waic", "SE"],
		p_WAIC = loo_obj$estimates["p_waic", "Estimate"],
		p_WAIC_SE = loo_obj$estimates["p_waic", "SE"],
		lppd = loo_obj$estimates["elpd_waic", "Estimate"] + loo_obj$estimates["p_waic", "Estimate"],
		lppd_SE = sqrt(var(lppd_emp_each) * T),
		lppd_analytic = lppd_dfz,
		lppd_analytic_SE = sqrt(var(lppd_each) * T)
		))
    #return(list(WAIC = WAIC, SE = SE))
}

waic.big_data_set.IG <- function(big_data_set, mean_files, left_cov_files,
				 nt_files, dt_files, R,
			         Vt = function(chunk) diag(nrow(chunk)),
			         NITER = 1000, verbose = FALSE, save_dir = NULL,
			         out.format = "csv", fst.head = "V") {

    update_ntdt <- (!is.numeric(nt_files)) & (length(nt_files) > 1)
    if (!update_ntdt) {
        # load them naively if they do not need to be updated.
        nt <- nt_files
        dt <- dt_files
    }

    T <- length(mean_files)

    loglikmat <- matrix(rep(-Inf, NITER * big_data_set$get("sum_K")), nrow=NITER)
    tk_ix <- 1

#    WAIC <- 0
    lppd_each <- rep(-Inf, big_data_set$get("sum_K"))
    p_waic2_each <- rep(Inf, big_data_set$get("sum_K"))
    WAIC_each <- rep(Inf, big_data_set$get("sum_K"))
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
        #WAIC_each_t <- c(Inf, K)

	for (k in 1:K) {
            if (verbose) cat(sprintf("Running time %d and block %d...", t, k))
            # read in Yk
	    ud_Y <- unlist(Y_file$grid[Y_file$grid$block == k, c("U", "D")])
	    lr_Y <- unlist(Y_file$grid[Y_file$grid$block == k, c("L", "R")])
	    Y_tk <- if (out.format == "csv") {
		        as.matrix(read_big_csv(file_name_Y, rows = ud_Y, cols = lr_Y,
		  	 		       has_header=Y_file$header,
					       ignore_header=TRUE))
	            } else if (out.format == "fst") {
                        as.matrix(read_big_fst(file_name_Y, rows = ud_Y, 
					       cols = paste0(fst.head, seq(lr_Y[1], lr_Y[2]))))
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
	    s_tk <- matrix(unlist(st_file[k+1,]), nrow=Fbncol, ncol=Ybncol)
	    S_tk <- tri.to.sym(unlist(St_file[k+1,]), Fbncol, lower=FALSE)
            
	    if (update_ntdt) {
                nt <- nt_list[k+1]
	        dt <- dt_list[k+1]
	    }

	    #TODO: Save Theta_tki in an Rdata object if save_dir is specified. (i.e. as waic_rinvgamma_samples_t%d_k%d.R)
            sigma2s <- 1/rgamma(NITER, nt, rate=dt)

	    #TODO: find some way to parallelize this over iterations. Ideally offload it to C? How to parallelize?

	    #TODO: How about offloading it to C?

            for (i in 1:NITER) {
	        #TODO: Save Theta_tki in an Rdata object if save_dir is specified. (i.e. as waic_rmatrixnorm_samples_t%d_k%d.R)
                Theta_tki <- rmn(1, s_tk, S_tk, sigma2s[i] * R)[,,1]
                loglikmat[i,tk_ix] <- dmn(Y_tk, F_tk %*% Theta_tki, Vt(F_tk), sigma2s[i] * R, log=TRUE)
            }

	    lppd <- dmatrixtscS(Y_tk, F_tk %*% s_tk, F_tk %*% S_tk %*% t(F_tk) + Vt(F_tk), 
	                        nt, dt, R)
            p_waic2 <- p_waic2_IG_1t(NITER, s_tk, S_tk, nt, dt, R, Y_tk, F_tk, Vt(F_tk))

	    lppd_each[tk_ix] <- lppd
	    p_waic2_each[tk_ix] <- p_waic2
	    WAIC_each[tk_ix] <- lppd - p_waic2
	    #WAIC_each_t[k] <- lppd - p_waic2
	    #if (!is.finite(WAIC_each_t[k])) break
	    tk_ix <- tk_ix + 1
	    if (verbose) cat(sprintf("Time %d episode %d completed.\n", t, k))
       	}
	#WAIC_each <- c(WAIC_each, WAIC_each_t)
    }

    #WAIC <- -2 *sum(WAIC_each)
    #SE <- sqrt(var(WAIC_each) * T) * 2

    if (!is.null(save_dir)) {
        #T save log-likelihood matrix in 
        saveRDS(loglikmat, file=file.path(save_dir, "loglikmat.rds"))
    }

    loo_obj <- loo::waic(loglikmat)

    lppd_dfz <- sum(lppd_each) 
    lppd_emp_each <- apply(loglikmat, 2, logSumExp) - log(NITER)
    lppd_dfz_emp <- sum(lppd_emp_each)
    p_waic2_each <- apply(loglikmat, 2, var)
    p_waic2_dfz <- sum(p_waic2_each)

    if (!is.null(save_dir)) {
        waic_df_perep <- data.frame(WAIC = -2 * (lppd_emp_each - p_waic2_each),
                                    SE = 2 * sqrt(var(lppd_emp_each - p_waic2_each)),
                                    p_WAIC = p_waic2_each,
                                    p_WAIC_SE = sqrt(var(p_waic2_each)),
                                    lppd = lppd_emp_each,
                                    lppd_SE = sqrt(var(lppd_emp_each)),
                                    lppd_analytic = lppd_each,
                                    lppd_analytic_SE = sqrt(var(lppd_each)))
        #T save log-likelihood matrix in 
        saveRDS(loglikmat, file=file.path(save_dir, "loglikmat.rds"))
        write.csv(waic_df_perep, file=file.path(save_dir, "WAIC_perep.csv"), row.names=FALSE)
    }

    loo_obj <- loo::waic(loglikmat)

    # TODO: return the table of estimates a'la loo?
#    out <- data.frame(col.names=c("Estimate", "SE"), row.names = c("waic", "p_waic", "elpd"))
#    out["waic"] <- loo_obj$estimates["waic", "Estimate"]

    return(list(WAIC = loo_obj$estimates["waic", "Estimate"], 
		SE = loo_obj$estimates["waic", "SE"],
		p_WAIC = loo_obj$estimates["p_waic", "Estimate"],
		p_WAIC_SE = loo_obj$estimates["p_waic", "SE"],
		lppd = loo_obj$estimates["elpd_waic", "Estimate"] + loo_obj$estimates["p_waic", "Estimate"],
		lppd_SE = sqrt(var(lppd_emp_each) * T),
		lppd_analytic = lppd_dfz,
		lppd_analytic_SE = sqrt(var(lppd_each) * T)
		))
    #return(list(WAIC = WAIC, SE = SE))
}

#TODO: update the below with fst as with waic.

# Compute the total lppd. Split it between the three options as before.
lppd.big_data_set <- function(object, mean_files, left_cov_files,
			      nt_files = NULL, Dt_files = NULL, Vt = function(chunk) diag(nrow(chunk)),
			      Sigma_mode = "IW",dt_files = NULL, R = NULL,
			      verbose = FALSE,
			      out.format = "csv",
			      fst.head = "V") {
    if((Sigma_mode == "IW") & ( is.null(nt_files) | is.null(Dt_files))) stop("nt_files and Dt_files must be provided for Sigma mode IW (inverse-Wishart).")
    if ((Sigma_mode == "IG") & (is.null(nt_files) | is.null(dt_files) | is.null(R))) stop("nt_files, dt_files, and R must be provided for Sigma mode IG (inverse-gamma)!")
    
    # Call the specific lppd function here.
    if (Sigma_mode == "id") {
        return(lppd.big_data_set.id(object, mean_files, left_cov_files, Vt = Vt,
			            verbose = verbose,
				    out.format = out.format, fst.head = fst.head))

    } else if (Sigma_mode == "IW") {
        return(lppd.big_data_set.IW(object, mean_files, left_cov_files,
			         nt_files, Dt_files, Vt = Vt,
				 verbose = verbose,
				 out.format = out.format, fst.head = fst.head))
    
    } else if (Sigma_mode == "IG") {
        return(lppd.big_data_set.IG(object, mean_files, left_cov_files,
				    nt_files, dt_files, R,
			            Vt = Vt,
				    verbose = verbose,
				    out.format = out.format, fst.head = fst.head))
    } else {
        stop(sprintf("Sigma mode %s not recognized!", Sigma_mode))
    }

}

lppd.big_data_set.id <- function(big_data_set, mean_files, left_cov_files, Vt = function(chunk) diag(nrow(chunk)),
			         verbose = FALSE,
			         out.format = "csv",
			         fst.head = "V") {
    T <- length(mean_files)

    lppd_each <- rep(-Inf, big_data_set$get("sum_K"))
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
					       cols = paste0(fst.head, seq(lr_Y[1], lr_Y[2]))))
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
	    s_tk <- matrix(unlist(st_file[k+1,]), nrow=Fbncol, ncol=Ybncol)
	    S_tk <- tri.to.sym(unlist(St_file[k+1,]), Fbncol, lower=FALSE)
 
            lppd_each[tk_ix] <- lppd_id_1t(Y_tk, F_tk, s_tk, S_tk, Vt(F_tk))

	    tk_ix <- tk_ix + 1
       	}
    }

    return(sum(lppd_each))
}


lppd.big_data_set.IW <- function(big_data_set, mean_files, left_cov_files,
			         nt_files, Dt_files, Vt = function(chunk) diag(nrow(chunk)),
				 verbose = FALSE,
			         out.format = "csv",
			         fst.head = "V") {
    update_ntDt <- (!is.numeric(nt_files)) & (length(nt_files) > 1)
    if (!update_ntDt) {
        # load them naively if they do not need to be updated.
        nt <- nt_files
        Dt <- Dt_files
    }

    T <- length(mean_files)

    lppd_each <- rep(-Inf, big_data_set$get("sum_K"))
    tk_ix <- 1
#    WAIC <- 0
#    WAIC_each <- c()
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
					       cols = paste0(fst.head, seq(lr_Y[1], lr_Y[2]))))
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
	    s_tk <- matrix(unlist(st_file[k+1,]), nrow=Fbncol, ncol=Ybncol)
	    S_tk <- tri.to.sym(unlist(St_file[k+1,]), Fbncol, lower=FALSE)

            if (update_ntDt) {
                nt <- nt_list[k+1]
	        Dt <- tri.to.sym(unlist(Dt_list[k+1,]), Ybncol, lower=FALSE)
	    }	    

            lppd_each[tk_ix] <- lppd_IW_1t(Y_tk, F_tk, s_tk, S_tk, Vt(F_tk), nt, Dt)

	    tk_ix <- tk_ix + 1
       	}
    }

    return(sum(lppd_each))
}

lppd.big_data_set.IG <- function(big_data_set, mean_files, left_cov_files,
				 nt_files, dt_files, R,
			         Vt = function(chunk) diag(nrow(chunk)),
			         verbose = FALSE, 
				 discount_factor=FALSE,
			         out.format = "csv",
			         fst.head = "V") {

    update_ntdt <- (!is.numeric(nt_files)) & (length(nt_files) > 1)
    if (!update_ntdt) {
        # load them naively if they do not need to be updated.
        nt <- nt_files
        dt <- dt_files
    }

    T <- length(mean_files)

    lppd_each <- rep(-Inf, big_data_set$get("sum_K"))
    tk_ix <- 1

#    WAIC <- 0
#    WAIC_each <- c()
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
        #WAIC_each_t <- c(Inf, K)

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
					       cols = paste0(fst.head, seq(lr_Y[1], lr_Y[2]))))
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
	    s_tk <- matrix(unlist(st_file[k+1,]), nrow=Fbncol, ncol=Ybncol)
	    S_tk <- tri.to.sym(unlist(St_file[k+1,]), Fbncol, lower=FALSE)
            
	    if (update_ntdt) {
                nt <- nt_list[k+1]
	        dt <- dt_list[k+1]
	    }

            lppd_each[tk_ix] <- lppd_IG_1t(Y_tk, F_tk, s_tk, S_tk, Vt(F_tk), nt, dt, R)

	    tk_ix <- tk_ix + 1
       	}
	#WAIC_each <- c(WAIC_each, WAIC_each_t)
    }
    return(sum(lppd_each))
}

#TODO: should I write a function the per-episode WAIC??? Yes, for now. But place it in its own script.


# Generate the analytic lppd
lppd_id_1t <- function (Y_t, F_t, s_t, S_t, V_t) {
    lppd <- dmn(Y_t, mu = F_t %*% s_t, 
		V = F_t %*% S_t %*% t(F_t) + V_t, 
		S = diag(ncol(Y_t)), log=TRUE)
    return(lppd)
}

lppd_IW_1t <- function(Y_t, F_t, s_t, S_t, V_t, n_t, D_t) {
    lppd <- dmatrixt(Y_t, F_t %*% s_t, F_t %*% S_t %*% t(F_t) + V_t, 
		     n_t, D_t, log=TRUE)
    return(lppd)
}

lppd_IG_1t <- function(Y_t, F_t, s_t, S_t, V_t, n_t, d_t, R) {
    lppd <- dmatrixtscS(Y_t, F_t %*% s_t, F_t %*% S_t %*% t(F_t) + V_t, 
                        n_t, d_t, R)
    return(lppd)
}

# p_waic2's to be computed for each time point with NITER times.
p_waic2_id_1t <- function(NITER, mt, Mt, Yt, Ft, Vt) {
    logpS <- rep(Inf, NITER)
    
    Theta_tki <- rmn(NITER, mt, Mt, diag(ncol(mt)))

    for (i in 1:NITER) {
        logpS[i] <- dmn(Yt, Ft %*% Theta_tki[,,i],
	                    Vt, diag(ncol(Yt)), log=TRUE)
    }

    return(var(logpS))
}

p_waic2_IW_1t <- function(NITER, mt, Mt, nt, Dt, Yt, Ft, Vt) {
    logpS <- rep(Inf, NITER)
    Sigma_tki <- rinvwishart(NITER, nu=nt, Omega=Dt, epsilon=0, checkSymmetry=F)

    for (i in 1:NITER) {
        Theta_tki <- rmatrixt(1, mt, Mt, nt, Dt)[,,1]
        logpS[i] <- dmn(Yt, Ft %*% Theta_tki, Vt, 
			      Sigma_tki[,,i], log=TRUE)
    } 

    return(var(logpS))
}

p_waic2_IG_1t <- function(NITER, mt, Mt, nt, dt, R, Yt, Ft, Vt) {
    logpS <- rep(Inf, NITER)

    sigma2s <- 1/rgamma(NITER, nt, rate=dt)

    for (i in 1:NITER) {
        Theta_tki <- rmn(1, mt, Mt, sigma2s[i] * R)[,,1]
        logpS[i] <- dmn(Yt, Ft %*% Theta_tki, Vt, sigma2s[i] * R, log=TRUE)
    }
    return(var(logpS))
}


