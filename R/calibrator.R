library(collections)
library(rlang)
library(mcmc)

##TODO: Write function to process big data F_tk Theta_tk. Specify how these data structures are arranged.
#
##' Process [F_tk %*% Theta_tk]_{k=1:K}, so that we end up with a single matrix for the episode-season model.
##' @param F_tk_df The data frame containing values of F_tk.
##' @param Theta_tk_df The data frame containing values of Theta_tk_df.
##' @param N The number of rows of F_t; the number of distinct sets of parameters used in the training set.
##' @param p The number of parameters encoded in F_t.
##' @param S The number of spatial coordinates encoded in Y_t, and also the number of rows in Theta_t.
##' @param bnrow The number of rows in one block of Y_t.
##' @param bncol The number of columns in one block Y_tk, also the number of rows in one block Theta_tk.
##' @returns Yhat, the matrix composed of horizontally-stacked elements [F_tk %*% Theta_tk]_{k=1,...,K}.
#
##TODO: Pass in the file containing F_tk instead. Sample a Theta_tk_df with a set seed here?
##TODO: Sample over spatial columns, specifically for Theta_tk_df. Subset with 
##TODO: We will need the grid here???
#big_data_stack_mult <- function(F_tk_df, st_file, St_file, N, p, S, bncol,
#				sp_coords_1ix, ) {
#    
#    Stilde <- length(sp_coords_1ix)
#
#    K <- N * S/bncol
#
#    Yhat <- matrix(nrow=N, ncol=S, data=-1)
#
#    #TODO: 
#
#    #TODO: What to do here???
#    for (k in 1:K) {
#        #TODO: Read in the relevant k-slice of F_t to get F_tk
#        #TODO: Sample 
#        
#    }
#
#    return(Yhat)
#}



#' Sample the values of tau2 for the Gibbs sampler.
#' @param n0z The shape parameter of the prior for tau^-2_0
#' @param d0z The initial rate parameter of the prior for tau^-2_0 (and consequently scale parameter for tau2_0).
#' @param b The "discount factor" for the priors for tau^2_1:T.
#' @param umat A (T+1) x S matrix of the values of u computed using the FFBS_calibration. Each row consists of u_t, with t ranging from 0 to T; each column contains the value of u_t(s_i) at the ith spatial index.
#' @param zdiffmat A T x S matrix of the values of z_t(s_i) - y_t(eta, s_i). Each row consists of one value of z_t - y_t(eta); each column contains the value of z_t(s_i) - y_t(eta, s_i) at the ith spatial index.
#' @param U The S x S correlation matrix for each u_t, sans scalar tau2_t.
#' @returns A single sample of each value of tau2_t.

sample_tau2s <- function(n0z, d0z, b, umat, zdiffmat, U) {
    T <- nrow(umat) - 1
    S <- ncol(umat)

    tau2s <- rep(-1, T+1)

    # Get tau2_0
    tau2s[1] <- 1/rgamma(1, shape=n0z, rate=d0z)

    uU <- chol(U)

    for (t in 1:T) {
	udiff <- umat[t+1,] - umat[t,]

        udiff <- backsolve(uU, udiff, transpose=TRUE)
#	.Call("trsm", 'L', 'U', 'T', 'N', S, 1, 1.0, uU, S, udiff, S)

	nt <- b^t * n0z + S

	zdiffbias <- zdiffmat[t,] - umat[t+1,]
	dt <- b^t * d0z + (sum(udiff^2) + sum(zdiffbias^2))/2

	#TODO: debug. Blows up after T = 100???
#	print(c(t, nt, dt))
#	print(c(sum(udiff^2), sum(zdiffbias^2)))

	tau2s[t+1] <- 1/rgamma(1, shape=nt, rate=dt)
    }
    
    return(tau2s)
}

#' Get the lth sample from each file containing Theta_t
#' @param Theta_files The list of files containing Theta_t. Each file contains L samples.
#' @param l The sample number to access. l <= L.
#' @param p The number of rows in each Theta_t
#' @param Stilde The number of columns in each Theta_t.
#' @param out.format The format of the files to read and write. Supports 'csv' and 'fst'. 'fst' is a fast file I/O format developed by Hadley Wickham. Defaults to 'csv'.


#TODO: Not actually being used in my use cases.
#TODO: Adjust for direct sampling??? Or do that elsewhere???
get_Theta_L <- function(Theta_files, l, p, Stilde, out.format = "csv") {
    T <- length(Theta_files)

    Theta_L <- matrix(nrow=T, ncol=p*Stilde, data=0)

    if (out.format == "csv") {
        for (t in 1:T) {
            Theta_L[t,] <- read.csv(Theta_files[t])[l,]
	}

    } else if (out.format == "fst") {
        for (t in 1:T) {
            Theta_L[t,] <- read.fst(Theta_files[t])[l,]
	}
    }
    
    return(Theta_L)
}

##TODO: Generate Bt and bt given one set of variables 
#gen_Btbt <- function(ft_eta, Theta_t, Jt, Vt, Ytrain, Ftrain, Sigma_block, tau2_t, zt, ut) {
#    scale_factor <- 0
#    Sigma_eta <- Sigma_block
#
#}


##TODO: code a minimal predict here. Needed to generate Ftilde for AR2.
#predict_calib <- function(Y_files, F_files,
#			  N, S, Stilde, p, K,
#			  Theta_L = NULL, st_files = NULL, St_files = NULL) {
#    #TODO: 
#    T <- length(Y_files)
#    sample_Theta_here <- (!is.null(st_files) & !is.null(St_files)) & is.null(Theta_L) 
#
#    #TODO: continue here.
#
#}

#' Generate y_t(eta) where eta is not one of the given parameters used for training the algorithm.

# toggles whether or not to generate the analytic y_t(eta). 
# TODO: Somehow the pure analytic case makes the jump probability worse??? especially for the prey population for LV: u_t.
analytic_Y <- FALSE #TRUE

# Fully analytic ODE generator. Use it here to help debug your MCMC.
gen_ys_analytic <- function(eta, eta_train, T, AR, sp_coords_1ix, sp_coords_1ix_lks, PDEini, pde_fun, n_PDEstates, grid_dims,
		   PDEsol_2_Y, gen_F,
		   seed_obj_TS = NULL, seed_obj_G = NULL) {

    #TODO: debug. Profiling this run
    t0 <- Sys.time()

    # Generate y1(eta_new) and y2(eta_new) using the PDE solver.
    #TODO: Is there a faster way to do this than with deSolve? We really only need one time point.
    #TODO: Should I just write my own function? Or specify a way that lets me only generate one time point???
    #NOTE: Rows are time, columns are spatial coordinates.

    is_spatial <- !is.null(grid_dims)
#    S <- grid_dims[1] * grid_dims[2]

    Y0s <- if (is_spatial) {
	       ode.2D(y = PDEini, parms = eta,
                     func = pde_fun,
                     nspec = n_PDEstates, dimens = grid_dims, times = seq(0,T),
                     lrw = 200000000, names=paste0("X", seq(1, n_PDEstates)))[,-1]
           } else {
               ode(y = PDEini, parms= eta,
	             func = pde_fun, times = seq(0,T))[,-1]
           }

    #TODO: debug
    t1 <- Sys.time()
    print(sprintf("ODE generation run: %.3f", t1 - t0))

    Stilde <- length(sp_coords_1ix)

    # If PDE cannot generate stable solutions.
    if (any(is.nan(Y0s))) {
        # break and return NaN
        print("Resulting parameters yielded invalid ODE initial solutions.")
        return(NaN)
    }

    # Transform into y's. e.g. Y = log(I + 1).
    Ytmat <- PDEsol_2_Y(Y0s)[2:(T+1), sp_coords_1ix]

    # If PDE cannot generate stable solutions.
    if (any(is.nan(Ytmat))) {
        # break and return NaN
        print("Resulting parameters yielded invalid ODE initial solutions.")
        return(NaN)
    }


    Btmat <- matrix(nrow=T-AR, ncol=Stilde*Stilde, data=0)
    Btbtmat <- Ytmat[(AR+1):T,]

    return(list(Y = Ytmat, Bt = Btmat, Btbt = Btbtmat,
		seed_obj_TS = seed_obj_TS, seed_obj_G = seed_obj_G))

}

# TODO: Whether or not to "cheat" when generating Ftilde. May help in developing and debugging.
#TODO: Do according to the paper.
cheat_Ft <- TRUE

#' Generate the Yt(eta) and Ft(eta) priors

gen_ys_prior <- function(eta, eta_train, sp_coords_1ix, sp_coords_1ix_lks, PDEini, pde_fun, n_PDEstates, grid_dims,
		   PDEsol_2_Y, gen_F, tau2s, z, u, Y_files, F_files, Sigma_block,
		   uVt, beta,
		   N, S, Stilde, p, K, AR = 2,
		   Theta_L = NULL,
		   st_files = NULL, St_files = NULL,
		   seed_obj_TS = NULL, seed_obj_G = NULL,
		   out.format = "csv",
		   fst.head = "V") {

    #TODO: mandate that Y0 be included in Y_files???
    T <- length(Y_files)
    #T <- nrow(z)
    
    sample_Theta_here <- (!is.null(st_files) & !is.null(St_files)) & is.null(Theta_L) 

    #TODO: Temporarily stopping direct sampling of Theta_tk.
    if (!sample_Theta_here) {
        stop("Temporarily rejecting existing Theta samples, especially for a large number of iterations. Supply the files containing the moments s_t and S_t generated with the emulator, as well as a range of preset seeds, to generate the set of Theta_tk's in a reproducible manner.")
    }

    #TODO: debug. Profiling this run
    t0 <- Sys.time()

    # Generate y1(eta_new) and y2(eta_new) using the PDE solver.
    #TODO: Is there a faster way to do this than with deSolve? We really only need one time point.
    #TODO: Should I just write my own function? Or specify a way that lets me only generate one time point???
    #NOTE: Rows are time, columns are spatial coordinates.

    is_spatial <- !is.null(grid_dims)
#    S <- grid_dims[1] * grid_dims[2]

    # NOTE: A minimum of 2 times is required for this PDE solver.
    Y0s <- tryCatch({
	       if (is_spatial) {
                   ode.2D(y = PDEini, parms = eta,
                          func = pde_fun,
                          nspec = n_PDEstates, dimens = grid_dims, times = c(0,1),
                          lrw = 200000000, names=paste0("X", seq(1, n_PDEstates)))[,-1]
               } else {
                   ode(y = PDEini, parms= eta,
    	             func = pde_fun, times = c(0,1))[,-1]
               }
           }, error = function(e) {
                  return(NaN)
              },
	      warning = function(e) {
                  return(NaN)
              }
           )

    #TODO: debug
    t1 <- Sys.time()
    print(sprintf("ODE generation run: %.3f", t1 - t0))

    # PDE cannot generate stable solutions.
    if (any(is.nan(Y0s))) {
        # break and return NaN
        print("Resulting parameters yielded invalid ODE initial solutions.")
        return(NaN)
    }

    # Transform into y's. e.g. Y = log(I + 1).
    Y0s <- PDEsol_2_Y(Y0s)

    # PDE cannot generate stable solutions.
    if (any(is.nan(Y0s))) {
        # break and return NaN
        print("Resulting parameters yielded invalid ODE initial solutions.")
        return(NaN)
    }
    #rbind(PDEsol_2_Y(matrix(PDEini, nrow=1)), PDEsol_2_Y(Y0s))

    bncol <- S/K
    fbncol <- p/K

    #TODO: We don't need to save all the Ftk. Just save the last Ftk.
    # NOTE: Ftmat if AR2 requires the full Y0s.
#    Ftkmat <- matrix(nrow=(T-AR)*K, ncol=fbncol, data=0)
    Ftkmat <- matrix(nrow=T*K, ncol=fbncol, data=0)
    # TODO: populate this instead??? Less storage overall. But you will not be able to "cheat".
    Ftk <- matrix(nrow=K, ncol=fbncol)

    for (k in 1:K) {
        Ftkmat[(AR-1)*K+k,] <- gen_F(AR+1, unlist(eta), Y0s[,((k-1)*bncol+1):(k*bncol)])
#        Ftkmat[AR*K+k,] <- gen_F(AR+2, unlist(eta), Y0s[,((k-1)*bncol+1):(k*bncol)])
        Ftk[k,] <- gen_F(AR+1, unlist(eta), Y0s[,((k-1)*bncol+1):(k*bncol)])
    }


    # to generate Ytmat_full, per episode, we need to store the full partitions per episode, but only consecutive ones
    ytm2 <- Y0s[1,] 
    ytm1 <- Y0s[2,] 

    # Subset with sp_coords
#    Y0s <- Y0s[, sp_coords_1ix]

    # Take advantage of the extra one generated.
    Ytmat <- matrix(nrow=T, ncol=Stilde, data=0)
    Ytmat[1,] <- Y0s[2, sp_coords_1ix]

    Jt <- exp(-beta * sqrt(rdist(eta_train, matrix(eta, nrow=1))))
    scale_factor <- duplicate(Jt, shallow=F)

#    .Call("trsm", 'L', 'U', 'T', 'N', N, 1, 1.0, uVt, N, scale_factor, N)

    scale_factor <- backsolve(uVt, Jt, transpose=TRUE)
    scale_factor <- 1 - sum(scale_factor^2) #dot(scale_factor, scale_factor)    

    Btmat <- matrix(nrow=T, ncol=Stilde*Stilde, data=0)
    Btbtmat <- matrix(nrow=T, ncol=Stilde, data=-1)

    if (scale_factor == 0) {
        # Skip the Bt and bt computation and set to the corresponding y, or at least the y corresponding to the closest distance numerically
        # get argmax of Jt
        eta_ix <- which.max(Jt)

#        for (t in (AR+1):T) {
#            Yt <- if (out.format == "csv") {
#		                 read.csv(Y_files[t])[eta_ix, sp_coords_1ix]
#	                     } else if (out.format == "fst") {
#                                 read.fst(Y_files[t])[eta_ix, sp_coords_1ix]
#	                     }

#        }

        t0 <- 1#AR
        k_ix0 <- 0 # track episode-wise indices
	for (t_k in (AR*K+1):((T+1)*K)) {
	    t <- floor((t_k - 1)/K)
            k <- ((t_k - 1) %% K) + 1

            if (t > t0) {
                # to accommodate big data, we will have to loop by K
                Ytmat[t0,] <- Btbtmat[t0,]

                k_ix0 <- 0
		t0 <- t
            }

            Y_tk_train <- if (out.format == "csv") {
                              as.matrix(read_big_csv(Y_files[t], cols = c((k-1)*bncol + 1, k*bncol),
                                                     has_header=TRUE,
                                                     ignore_header=TRUE))[eta_ix,]
                          } else if (out.format == "fst") {
                              as.matrix(read_big_fst(Y_files[t], 
                                                     cols = paste0(fst.head,
    				   	             seq((k-1)*bncol + 1, k*bncol))))[eta_ix,]
                          }

            l_k <- sp_coords_1ix_lks[k] #length(sp_coords_1ix_k)

            if (l_k == 0) next

            sp_coords_1ix_k <- sp_coords_1ix[(k_ix0 + 1):(k_ix0 + l_k)] - bncol * (k-1)
            Btbtmat[t,(k_ix0 + 1):(k_ix0 + l_k)] <- Y_tk_train[sp_coords_1ix_k]

	    k_ix0 <- k_ix0 + l_k
	}

    } else {
        #TODO: debug and development only?
        if (cheat_Ft) {
            #TODO: try/catch this the same as with the data generation.
            Yt_eta <- if (is_spatial) {
        	          ode.2D(y = PDEini, parms = eta, func = pde_fun,
                             nspec = n_PDEstates, dimens = grid_dims, times = seq(0,T),
                             lrw = 200000000, names=paste0("X", seq(1, n_PDEstates)))[,-1]
                      } else {
                          ode(y = PDEini, parms = eta,
        	              func = pde_fun, times = seq(0,T))[,-1]
                      }

            Yt_eta <- PDEsol_2_Y(Yt_eta)

            for (t_k in (AR*K+1):(T*K)) {
	        t <- floor((t_k - 1)/K)
                k <- ((t_k - 1) %% K) + 1
	        #t <- floor(t_k/K)
                #k <- (t_k %% K) + 1

                Ftkmat[t*K+k,] <- gen_F(t+AR+1, eta, Yt_eta[,((k-1)*bncol+1):(k*bncol)])
#		Ftk <- gen_F(t+1, eta, Yt_eta[,((k-1)*bncol+1):(k*bncol)])
	    }

	}

        # Construct Sigmatilde_eta outside of loop. The locations to sample should be the same over time.
        Sigma_eta <- scale_factor * Sigma_block
        Sigmatilde_eta <- matrix(nrow=Stilde, ncol=Stilde, data=0)
    
        k_ix0 <- 0
        for (k in 1:K) {
            l_k <- sp_coords_1ix_lks[k] #length(sp_coords_1ix_k)
    
            if (l_k == 0) next
    
            sp_coords_1ix_k <- sp_coords_1ix[(k_ix0 + 1):(k_ix0 + l_k)] - bncol * (k-1)
    		    #t(Ftkmat[(t-2)*K+(k-1),] %*%
        		                                           #Theta_tk_l[, sp_coords_1ix_k])
    #		    mu_tilde_t[(k_ix0 + 1):(k_ix0 + l_k)] <- mu_tilde_t[(k_ix0 + 1):(k_ix0 + l_k)] +
    #			    Ytkerr[,sp_coords_1ix_k]
      
            Sigmatilde_eta[(k_ix0 + 1):(k_ix0 + l_k), (k_ix0 + 1):(k_ix0 + l_k)] <-
                    Sigma_eta[sp_coords_1ix_k, sp_coords_1ix_k]
        
            k_ix0 <- k_ix0 + l_k
        }
     
        uSigmatilde_eta <- chol(Sigmatilde_eta)
        #TODO: Run the predictor here to generate the remaining y's.

        t0 <- AR-1
        k_ix0 <- 0

        # one index for both t and k? Will that speed things up?
        for (t_k in (AR*K + 1):((T+1)*K)) {
	    t <- floor((t_k - 1)/K) #+ AR
            k <- ((t_k - 1) %% K) + 1
	    #t <- floor(t_k/K)
            #k <- (t_k %% K) + 1

        #TODO: Compare with generating the full Y through ODE. At the rate I am going, would it seriously be slower???

#        for (t in 3:T) {
            if (t > t0) {
                ytm2 <- duplicate(ytm1, shallow=FALSE)
                if (sample_Theta_here) {
                    mu_tilde_t <- rep(0, Stilde) 
      
		    # t ranges from 3 to 50. st ranges from 1 to 49. Align st so that s_T aligns with T.
                    st_file <- if (out.format == "csv") {
                                   read.csv(st_files[t])
                               } else if (out.format == "fst") {
                                   read.fst(st_files[t])
                               }
        
                    St_file <- if (out.format == "csv") {
                                   read.csv(St_files[t])
                               } else if (out.format == "fst") {
                                   read.fst(St_files[t])
                               }
		}
                k_ix0 <- 0
		t0 <- t
            }
    
            if (sample_Theta_here) {
                #TODO: Store the full N x S Yterr and build the full y, not just ytilde. 
#                Yterr_full <- matrix(nrow=1, ncol=S, data=0)	
    
                #TODO: copy this from the Btbt and Bt initialization case before the loop
#                Yterr <- matrix(nrow=N, ncol=Stilde, data=0)
    
                # construct Yterr and mu_tilde_t from the episode-season model
 #               for (k in 1:K) {
                    Y_tk_train <- if (out.format == "csv") {
                                      as.matrix(read_big_csv(Y_files[t], cols = c((k-1)*bncol + 1, k*bncol),
                                                             has_header=TRUE,
                                                             ignore_header=TRUE))
                                  } else if (out.format == "fst") {
                                      as.matrix(read_big_fst(Y_files[t], 
                                                             cols = paste0(fst.head,
    								       seq((k-1)*bncol + 1, k*bncol))))
                                  }
    
                    F_tk_train <- if (out.format == "csv") {
                                      as.matrix(read_big_csv(F_files[t], cols = c((k-1)*fbncol+1, k*fbncol),
                                                             has_header=TRUE,
                                                             ignore_header=TRUE))
                                  } else if (out.format == "fst") {
                                      as.matrix(read_big_fst(F_files[t], 
                                                             cols = paste0(fst.head,
          								 seq((k-1)*fbncol+1, k*fbncol))))
                                  }
    
                    s_tk <- matrix(unlist(st_file[k+1,]), nrow=fbncol, ncol=bncol)
                    S_tk <- tri.to.sym(unlist(St_file[k+1,]), fbncol, lower=FALSE)
    
                    seed_obj_TS <- update_seed_newcopy(seed_obj_TS, seed_obj_TS$seed_step)

                    Theta_tk_l <- rmn(1, s_tk, S_tk, Sigma_block)[,,1]
#		    as.matrix(rmn(1, s_tk, S_tk, Sigma_block)[,,1],
#    		                            nrow=fbncol, ncol=bncol)
  
		    #TODO: Offload it to pre-computation??? 
                    Ytkerr <- Y_tk_train - F_tk_train %*% Theta_tk_l

		    Ytkerr <- backsolve(uVt, backsolve(uVt, Ytkerr, transpose=TRUE))
#                    .Call("trsm", 'L', 'U', 'T', 'N', N, bncol, 1.0, uVt, N, Ytkerr, N)
#                    .Call("trsm", 'L', 'U', 'N', 'N', N, bncol, 1.0, uVt, N, Ytkerr, N)

                    mu_tilde_tk_full <- t(Ftkmat[(t-1)*K+k,] %*% Theta_tk_l + 
			                  crossprod(Jt, Ytkerr))
#                    mu_tilde_tk_full <- t(Ftk[k,] %*% Theta_tk_l + 
#			                  crossprod(Jt, Ytkerr))

#Yterr_full[,((k-1)*bncol+1):(k*bncol)]
    
#                    seed_obj_G <- update_seed_newcopy(seed_obj_G, seed_obj_G$seed_step)

		    # TODO: Xiang uses the actual PDE solutions to generate Ftkmat and its values. This is technically cheating but should we just do it???
		    # TODO: AR2-specific implementation. Generalize later.
                    # sample the full set of y_tk's for Ftkmat.
#                    ytm1[((k-1)*bncol+1):(k*bncol)] <- rmvnorm(1, mu_tilde_tk_full, Sigma_eta)
                    ytm1[((k-1)*bncol+1):(k*bncol)] <- mu_tilde_tk_full
 
		    # NOTE: No need to store past this time. 
                    # sample the new y_tk and store part of y_tk in ytm1
                    if (!cheat_Ft) Ftkmat[t*K+k,] <- gen_F(AR+1, eta, 
    					                    rbind(ytm2[((k-1)*bncol+1):(k*bncol)],
    						            ytm1[((k-1)*bncol+1):(k*bncol)]))
                    if (!cheat_Ft) Ftk[k,] <- gen_F(AR+1, eta, 
    					            rbind(ytm2[((k-1)*bncol+1):(k*bncol)],
    						          ytm1[((k-1)*bncol+1):(k*bncol)]))
    
		    #TODO: We can compute these outside the Gibbs loop and store them in a list or dictionary,
		    #      and then access them for each k.  
                    # Store the coordinate subsets
#                    sp_coords_1ix_k <- sp_coords_1ix_k_all$get(k)
#			    sp_coords_1ix[(sp_coords_1ix < k * bncol) & 
#						     (sp_coords_1ix >= (k-1) * bncol + 1)] - 
#						     bncol * (k-1)
#                    l_k <- sum((sp_coords_1ix < k * bncol) & (sp_coords_1ix >= (k-1) * bncol + 1))
		    l_k <- sp_coords_1ix_lks[k] #length(sp_coords_1ix_k)

                    if (l_k == 0) next

                    sp_coords_1ix_k <- sp_coords_1ix[(k_ix0 + 1):(k_ix0 + l_k)] - bncol * (k-1)
                    mu_tilde_t[(k_ix0 + 1):(k_ix0 + l_k)] <- mu_tilde_tk_full[sp_coords_1ix_k]

#		    #t(Ftkmat[(t-2)*K+(k-1),] %*%
#    		                                           #Theta_tk_l[, sp_coords_1ix_k])
##		    mu_tilde_t[(k_ix0 + 1):(k_ix0 + l_k)] <- mu_tilde_t[(k_ix0 + 1):(k_ix0 + l_k)] +
##			    Ytkerr[,sp_coords_1ix_k]
#  
#                    Sigmatilde_eta[(k_ix0 + 1):(k_ix0 + l_k), (k_ix0 + 1):(k_ix0 + l_k)] <-
#                            Sigma_eta[sp_coords_1ix_k, sp_coords_1ix_k]
#    
                    k_ix0 <- k_ix0 + l_k
#                }
    
    #            Yterr_full <- crossprod(backsolve(uVt, backsolve(uVt, Yterr_full, transpose=TRUE)), Jt)
    
#                mu_tilde_t <- mu_tilde_t + Yterr_full[, sp_coords_1ix]
            } else {
                #TODO: Only applicable where the whole Theta_L is given, not episode-season. Re-adapt here to episode-season?  
                Theta_t_l <- matrix(Theta_L[t,], nrow=p, ncol=S)[,sp_coords_1ix]
                Yterr <- crossprod(backsolve(uVt, backsolve(uVt, Yt_train - Ft_train %*% Theta_t_l, transpose=TRUE)), Jt)[,sp_coords_1ix]
                #TODO: Theta_t_l %*% Ftmat[t-2,] will need to be processed separately, especially for the big data case.
                mu_tilde_t <- t(matrix(t(Ftkmat[((t-1)*K+(k-1)):((t-1)*K),]), nrow=1) %*% Theta_t_l) + Yterr
            }
   
            #TODO: Do the below after all k's have had their run

	    if (k_ix0 == Stilde) {
                Btbtmat[t,] <- mu_tilde_t
                Btmat[t,] <- c(Sigmatilde_eta)
        
                seed_obj_G <- update_seed_newcopy(seed_obj_G, seed_obj_G$seed_step)
                
                # Yt must be sampled from N(Bt * bt, Bt)
                Ytmat[t,] <- rmvnorm(1, mu_tilde_t, Sigmatilde_eta)
	    }
        }
    }

    return(list(Y = Ytmat, #F = Ftkmat,
		Bt = Btmat, Btbt = Btbtmat, 
		seed_obj_TS = seed_obj_TS, seed_obj_G = seed_obj_G))


}



#' Generate the values of Yt(eta) and Ft(eta)
#' @param eta The value of eta for which to generate the new values of y_t. A named list with each element of the list equal to the variable of the PDE in question.
#' @param eta_train The N x d matrix of eta's used to train the emulator. 
#' @param sp_coords_1ix The single-index spatial coordinates with which to index the PDE solutions.
#' @param sp_coords_1ix_lks An array of the number of spatial coordinates per block.
#' @param PDEini A grid_dim[1] x (n_PDEstates * grid_dim[2]) matrix containing the initial spatial configurations of each state of the PDE.
#' @param pde_fun The PDE function to be passed into the PDE solver to generate synthetic data.
#' @param n_PDEstates The number of states in the PDE to keep track of.
#' @param grid_dims A 2-element vector containing the dimensions of the full grid of the PDE function. Note that the product of the entries of grid_dims must be >= Stilde.
#' @param PDEsol_2_Y A function that takes the time index and the raw solutions of the PDE and outputs them to the desired value of Y.
#' @param gen_F A function that takes in the time index, parameter values eta, and Y values to produce a value of F_t as a set of arguments in that order. This function will be used to produce F_t-tilde, corresponding to Y(eta-tilde).
#' @param tau2s A length-T array containing the values of tau2_t sampled at the tau2_t Gibbs step from t = 1 to T.
#' @param z The T x Stilde field data.
#' @param u The (T + 1) x Stilde bias terms computed with FFBS_calibration at the u_t Gibbs step.
#' @param Y_files The list of Yt files used to train the emulator.
#' @param F_files The list of Ft files used to train the emulator.
#' @param Sigma_block The sampled matrix from the IW or IG * scale matrix generated outside of gen_ys.
#' @param uVt 
#' @param beta 
#' @param N The number of different sets of eta used to train the emulator.
#' @param S 
#' @param Stilde The number of spatial coordinates to calibrate to.
#' @param p The number of rows in each Theta_t.
#' @param K The number of "episodes" per "season" t.
#' @param AR The autoregression number used in the data and model. Required to adjust the calibrator to a certain offset. Default: 2.
#' @param Theta_L A T x (p * Stilde) matrix containing the lth sampled values of Theta_t. Each row contains the flattened lth sample at time index t. 
#' @param st_files The array of files containing the T posterior means generated from smoothing. Relevant where L samples each of ThetaSt needs to be generated in the calibrator. Default: NULL.
#' @param St_files The array of files containing the T posterior row-covariance matrices generated from smoothing. Relevant where L samples each of Theta_t needs to be generated in the calibrator. Default: NULL.
#' @param seed_obj_TS The seed object to facilitate random number generation for Theta_kt and Sigma (or sigma^2).
#' @param seed_obj_G The seed object to facilitate random number generation for y_tk.
#' @param out.format The output format of the files for the calibrator. Supports both 'csv' and 'fst'. Defaults to 'csv'.
#' @returns A list containing yt and ft flattened into matrices, with each row contaninig one yt or ft, and the random number generator, to maintain control over randomness, where relevant.

gen_ys <- function(eta, eta_train, sp_coords_1ix, sp_coords_1ix_lks, PDEini, pde_fun, n_PDEstates, grid_dims,
		   PDEsol_2_Y, gen_F, tau2s, z, u, Y_files, F_files, Sigma_block,
		   uVt, beta,
		   N, S, Stilde, p, K, AR = 2,
		   Theta_L = NULL,
		   st_files = NULL, St_files = NULL,
		   seed_obj_TS = NULL, seed_obj_G = NULL,
		   out.format = "csv",
		   fst.head = "V") {

#    T <- length(Y_files)
    T <- nrow(z)

    sample_Theta_here <- (!is.null(st_files) & !is.null(St_files)) & is.null(Theta_L) 

    #TODO: Temporarily stopping direct sampling of Theta_tk.
    if (!sample_Theta_here) {
        stop("Temporarily rejecting existing Theta samples, especially for a large number of iterations. Supply the files containing the moments s_t and S_t generated with the emulator, as well as a range of preset seeds, to generate the set of Theta_tk's in a reproducible manner.")
    }

    #TODO: debug. Profiling this run
    t0 <- Sys.time()

    # Generate y1(eta_new) and y2(eta_new) using the PDE solver.
    #TODO: Is there a faster way to do this than with deSolve? We really only need one time point.
    #TODO: Should I just write my own function? Or specify a way that lets me only generate one time point???
    #NOTE: Rows are time, columns are spatial coordinates.

    is_spatial <- !is.null(grid_dims)
#    S <- grid_dims[1] * grid_dims[2]

    # NOTE: A minimum of 2 times is required for this PDE solver.
    Y0s <- tryCatch({
              if (is_spatial) {
                  ode.2D(y = PDEini, parms = eta,
                         func = pde_fun,
                         nspec = n_PDEstates, dimens = grid_dims, times = c(0,1),
                         lrw = 200000000, names=paste0("X", seq(1, n_PDEstates)))[,-1]
              } else {
                  ode(y = PDEini, parms= eta,
	                 func = pde_fun, times = c(0,1))[,-1]
              }
           }, error = function(e) {
                  return(NaN)
               },
	       warning = function(e) {
                  return(NaN)
               }
           )

    #TODO: debug
    t1 <- Sys.time()
    print(sprintf("ODE generation run: %.3f", t1 - t0))

    # PDE cannot generate stable solutions.
    if (any(is.nan(Y0s))) {
        # break and return NaN
        print("Resulting parameters yielded invalid ODE initial solutions.")
        return(NaN)
    }

    # Transform into y's. e.g. Y = log(I + 1).
    Y0s <- PDEsol_2_Y(Y0s)

    # If PDE cannot generate stable solutions.
    if (any(is.nan(Y0s))) {
        # break and return NaN
        print("Resulting parameters yielded invalid ODE initial solutions.")
        return(NaN)
    }
	    #rbind(PDEsol_2_Y(matrix(PDEini, nrow=1)), PDEsol_2_Y(Y0s))

    bncol <- S/K
    fbncol <- p/K

    #TODO: We don't need to save all the Ftk. Just save the last Ftk.
    # NOTE: Ftmat if AR2 requires the full Y0s.
#    Ftkmat <- matrix(nrow=(T-AR)*K, ncol=fbncol, data=0)
    Ftkmat <- matrix(nrow=T*K, ncol=fbncol, data=0)
    #TODO: populate this instead??? Less storage overall.
    Ftk <- matrix(nrow=K, ncol=fbncol)
    for (k in 1:K) {
        Ftkmat[k,] <- gen_F(AR+1, unlist(eta), Y0s[,((k-1)*bncol+1):(k*bncol)])
        Ftk[k,] <- gen_F(AR+1, unlist(eta), Y0s[,((k-1)*bncol+1):(k*bncol)])
    }


    # to generate Ytmat_full, per episode, we need to store the full partitions per episode, but only consecutive ones
    ytm2 <- Y0s[1,] 
    ytm1 <- Y0s[2,] 

    # Subset with sp_coords
#    Y0s <- Y0s[, sp_coords_1ix]

    # Take advantage of the extra one generated.
    Ytmat <- matrix(nrow=T, ncol=Stilde, data=0)
#    Ytmat[1:AR,] <- Y0s[1:(AR), sp_coords_1ix]
    Ytmat[1,] <- Y0s[2, sp_coords_1ix]

    Jt <- exp(-beta * sqrt(rdist(eta_train, matrix(eta, nrow=1))))
    scale_factor <- duplicate(Jt, shallow=F)

#    .Call("trsm", 'L', 'U', 'T', 'N', N, 1, 1.0, uVt, N, scale_factor, N)

    scale_factor <- backsolve(uVt, Jt, transpose=TRUE)
    scale_factor <- 1 - sum(scale_factor^2) #dot(scale_factor, scale_factor)
#    scale_factor <- 1 - crossprod(scale_factor)[1,1]

    #NOTE: Based on the setup of our PDE generator, but is certainly AR2-specific.
    Btmat <- matrix(nrow=T, ncol=Stilde*Stilde, data=0)
    Btbtmat <- matrix(nrow=T, ncol=Stilde, data=-1)

    #TODO: debug
    t0 <- Sys.time()

    #TODO: implement this also for intro, so we can progress from there.
    if (scale_factor == 0) {
        # Skip the Bt and bt computation and set to the corresponding y, or at least the y corresponding to the closest distance numerically
        # get argmax of Jt
        eta_ix <- which.max(Jt)

#        for (t in (AR+1):T) {
#            Yt <- if (out.format == "csv") {
#		                 read.csv(Y_files[t])[eta_ix, sp_coords_1ix]
#	                     } else if (out.format == "fst") {
#                                 read.fst(Y_files[t])[eta_ix, sp_coords_1ix]
#	                     }

#        }

        t0 <- AR 
        k_ix0 <- 0
	for (t_k in (AR*K+1):((T+1)*K)) {
	    t <- floor((t_k - 1)/K)
            k <- ((t_k - 1) %% K) + 1

            if (t > t0) {
                # to accommodate big data, we will have to loop by K
                Ytmat[t0,] <- Btbtmat[t0,]

                k_ix0 <- 0
		t0 <- t
            }

            Y_tk_train <- if (out.format == "csv") {
                              as.matrix(read_big_csv(Y_files[t], cols = c((k-1)*bncol + 1, k*bncol),
                                                     has_header=TRUE,
                                                     ignore_header=TRUE))[eta_ix,]
                          } else if (out.format == "fst") {
                              as.matrix(read_big_fst(Y_files[t], 
                                                     cols = paste0(fst.head,
    				   	             seq((k-1)*bncol + 1, k*bncol))))[eta_ix,]
                          }

            l_k <- sp_coords_1ix_lks[k] #length(sp_coords_1ix_k)

            if (l_k == 0) next

            sp_coords_1ix_k <- sp_coords_1ix[(k_ix0 + 1):(k_ix0 + l_k)] - bncol * (k-1)
            Btbtmat[t,(k_ix0 + 1):(k_ix0 + l_k)] <- Y_tk_train[sp_coords_1ix_k]

	    k_ix0 <- k_ix0 + l_k
	}

    } else {
        #TODO: debug and development only?
        if (cheat_Ft) {
            Yt_eta <- tryCatch({
        	          if (is_spatial) {
                              ode.2D(y = PDEini, parms = eta,
                                     func = pde_fun,
                                     nspec = n_PDEstates, dimens = grid_dims, times = seq(0,T),
                                     lrw = 200000000, names=paste0("X", seq(1, n_PDEstates)))[,-1]
                          } else {
                              ode(y = PDEini, parms= eta,
            	                 func = pde_fun, times = seq(0,T))[,-1]
                          }
                      }, error = function(e) {
                          return(NaN)
                      },
        	      warning = function(e) {
                          return(NaN)
                      }
                   )

            # PDE cannot generate stable solutions.
            if (any(is.nan(Yt_eta))) {
                # break and return NaN
                print("Resulting parameters yielded invalid ODE initial solutions where we \"cheat\" F_t(eta). Either don't \"cheat\" or pick new parameters.")
                return(NaN)
            }

            Yt_eta <- PDEsol_2_Y(Yt_eta)

            if (any(is.nan(Yt_eta))) {
                # break and return NaN
                print("Resulting parameters yielded invalid ODE initial solutions where we \"cheat\" F_t(eta). Either don't \"cheat\" or pick new parameters.")
                return(NaN)
            }
	    #rbind(PDEsol_2_Y(matrix(PDEini, nrow=1)), PDEsol_2_Y(Y0s))

	    # TODO: I think it's an indexing problem for F_t
	    # Note that the first Ftkmat index has already been taken care of.
            for (t_k in (AR*K+1):((T+1)*K)) {
	        t <- floor((t_k - 1)/K)
                k <- ((t_k - 1) %% K) + 1

                Ftkmat[(t-1)*K+k,] <- gen_F(t+AR, eta, Yt_eta[,((k-1)*bncol+1):(k*bncol)])
	    }
	}

        # Construct Sigmatilde_eta outside of loop. The locations to sample should be the same over time.
        Sigma_eta <- scale_factor * Sigma_block
        Sigmatilde_eta <- matrix(nrow=Stilde, ncol=Stilde, data=0)
   
        k_ix0 <- 0
        for (k in 1:K) {
            l_k <- sp_coords_1ix_lks[k] #length(sp_coords_1ix_k)
    
            if (l_k == 0) next
    
            sp_coords_1ix_k <- sp_coords_1ix[(k_ix0 + 1):(k_ix0 + l_k)] - bncol * (k-1)
    		    #t(Ftkmat[(t-2)*K+(k-1),] %*%
        		                                           #Theta_tk_l[, sp_coords_1ix_k])
    #		    mu_tilde_t[(k_ix0 + 1):(k_ix0 + l_k)] <- mu_tilde_t[(k_ix0 + 1):(k_ix0 + l_k)] +
    #			    Ytkerr[,sp_coords_1ix_k]
      
            Sigmatilde_eta[(k_ix0 + 1):(k_ix0 + l_k), (k_ix0 + 1):(k_ix0 + l_k)] <-
                    Sigma_eta[sp_coords_1ix_k, sp_coords_1ix_k]
        
            k_ix0 <- k_ix0 + l_k
        }
     
        uSigmatilde_eta <- chol(Sigmatilde_eta)
        #TODO: Run the predictor here to generate the remaining y's.

        # one index for both t and k? Will that speed things up?
        t0 <- AR - 1
        k_ix0 <- 0

	#TODO: make sure that the time point samples 1 as well.
        for (t_k in (AR*K+1):((T+1)*K)) {
	    t <- floor((t_k - 1)/K)
            k <- ((t_k - 1) %% K) + 1

        #TODO: Compare with generating the full Y through ODE. At the rate I am going, would it seriously be slower???

#        for (t in 3:T) {
            if (t > t0) {
                ytm2 <- duplicate(ytm1, shallow=FALSE)
                if (sample_Theta_here) {
                    mu_tilde_t <- rep(0, Stilde) 
      
		    # t ranges from 3 to 50. st ranges from 1 to 49. Align st so that s_T aligns with T.
                    st_file <- if (out.format == "csv") {
                                   read.csv(st_files[t])
                               } else if (out.format == "fst") {
                                   read.fst(st_files[t])
                               }
        
                    St_file <- if (out.format == "csv") {
                                   read.csv(St_files[t])
                               } else if (out.format == "fst") {
                                   read.fst(St_files[t])
                               }
		}
                k_ix0 <- 0
		t0 <- t
            }
    
            if (sample_Theta_here) {
                #TODO: Store the full N x S Yterr and build the full y, not just ytilde. 
#                Yterr_full <- matrix(nrow=1, ncol=S, data=0)	
    
                #TODO: copy this from the Btbt and Bt initialization case before the loop
#                Yterr <- matrix(nrow=N, ncol=Stilde, data=0)
    
                # construct Yterr and mu_tilde_t from the episode-season model
 #               for (k in 1:K) {
                    Y_tk_train <- if (out.format == "csv") {
                                      as.matrix(read_big_csv(Y_files[t], cols = c((k-1)*bncol + 1, k*bncol),
                                                             has_header=TRUE,
                                                             ignore_header=TRUE))
                                  } else if (out.format == "fst") {
                                      as.matrix(read_big_fst(Y_files[t], 
                                                             cols = paste0(fst.head,
    								       seq((k-1)*bncol + 1, k*bncol))))
                                  }
   
                    F_tk_train <- if (out.format == "csv") {
                                      as.matrix(read_big_csv(F_files[t], cols = c((k-1)*fbncol+1, k*fbncol),
                                                             has_header=TRUE,
                                                             ignore_header=TRUE))
                                  } else if (out.format == "fst") {
                                      as.matrix(read_big_fst(F_files[t], 
                                                             cols = paste0(fst.head,
          								 seq((k-1)*fbncol+1, k*fbncol))))
                                  }
    
                    s_tk <- matrix(unlist(st_file[k+1,]), nrow=fbncol, ncol=bncol)
                    S_tk <- tri.to.sym(unlist(St_file[k+1,]), fbncol, lower=FALSE)
    
                    seed_obj_TS <- update_seed_newcopy(seed_obj_TS, seed_obj_TS$seed_step)

		    #TODO: Theta_tk_l is sampled here and the loop should end here. Otherwise, read from file and get the next iteration of the loop from file. (requires Gibbs index)
                    Theta_tk_l <- rmn(1, s_tk, S_tk, Sigma_block)[,,1]
#		    as.matrix(rmn(1, s_tk, S_tk, Sigma_block)[,,1],
#    		                            nrow=fbncol, ncol=bncol)
  
		    #TODO: Offload it to pre-computation??? 
                    Ytkerr <- Y_tk_train - F_tk_train %*% Theta_tk_l

		    Ytkerr <- backsolve(uVt, backsolve(uVt, Ytkerr, transpose=TRUE))
#                    .Call("trsm", 'L', 'U', 'T', 'N', N, bncol, 1.0, uVt, N, Ytkerr, N)
#                    .Call("trsm", 'L', 'U', 'N', 'N', N, bncol, 1.0, uVt, N, Ytkerr, N)

#                    Yterr_full[,((k-1)*bncol+1):(k*bncol)] <- crossprod(backsolve(uVt,
#    									backsolve(uVt,
#    							   		    Y_tk_train - F_tk_train %*% Theta_tk_l,
#    									    transpose=TRUE)), Jt)
   
#                    mu_tilde_tk_full <- t(Ftkmat[(t-AR+1)*K+k-1,] %*% Theta_tk_l + 
#			                  crossprod(Jt, Ytkerr))

                    mu_tilde_tk_full <- t(Ftkmat[(t-1)*K+k,] %*% Theta_tk_l + 
			                  crossprod(Jt, Ytkerr))

#Yterr_full[,((k-1)*bncol+1):(k*bncol)]
    
#                    seed_obj_G <- update_seed_newcopy(seed_obj_G, seed_obj_G$seed_step)

		    # TODO: Xiang uses the actual PDE solutions to generate Ftkmat and its values. This is technically cheating but should we just do it???
		    # TODO: AR2-specific implementation. Generalize later.
                    # sample the full set of y_tk's for Ftkmat.
#                    ytm1[((k-1)*bncol+1):(k*bncol)] <- rmvnorm(1, mu_tilde_tk_full, Sigma_eta)
                    ytm1[((k-1)*bncol+1):(k*bncol)] <- mu_tilde_tk_full
   
                    # sample the new y_tk and store part of y_tk in ytm1
                    if ((!cheat_Ft) & (t < T - AR)) Ftkmat[t*K+k,] <- gen_F(AR+1, eta, 
    					            rbind(ytm2[((k-1)*bncol+1):(k*bncol)],
    						          ytm1[((k-1)*bncol+1):(k*bncol)]))
                    if (!cheat_Ft) Ftk[k,] <- gen_F(AR+1, eta, 
    					            rbind(ytm2[((k-1)*bncol+1):(k*bncol)],
    						          ytm1[((k-1)*bncol+1):(k*bncol)]))
    
		    #TODO: We can compute these outside the Gibbs loop and store them in a list or dictionary,
		    #      and then access them for each k.  
                    # Store the coordinate subsets
#                    sp_coords_1ix_k <- sp_coords_1ix_k_all$get(k)
#			    sp_coords_1ix[(sp_coords_1ix < k * bncol) & 
#						     (sp_coords_1ix >= (k-1) * bncol + 1)] - 
#						     bncol * (k-1)
#                    l_k <- sum((sp_coords_1ix < k * bncol) & (sp_coords_1ix >= (k-1) * bncol + 1))
		    l_k <- sp_coords_1ix_lks[k] #length(sp_coords_1ix_k)

                    if (l_k == 0) next

                    sp_coords_1ix_k <- sp_coords_1ix[(k_ix0 + 1):(k_ix0 + l_k)] - bncol * (k-1)
                    mu_tilde_t[(k_ix0 + 1):(k_ix0 + l_k)] <- mu_tilde_tk_full[sp_coords_1ix_k]

#		    #t(Ftkmat[(t-2)*K+(k-1),] %*%
#    		                                           #Theta_tk_l[, sp_coords_1ix_k])
##		    mu_tilde_t[(k_ix0 + 1):(k_ix0 + l_k)] <- mu_tilde_t[(k_ix0 + 1):(k_ix0 + l_k)] +
##			    Ytkerr[,sp_coords_1ix_k]
#  
#                    Sigmatilde_eta[(k_ix0 + 1):(k_ix0 + l_k), (k_ix0 + 1):(k_ix0 + l_k)] <-
#                            Sigma_eta[sp_coords_1ix_k, sp_coords_1ix_k]
#    
                    k_ix0 <- k_ix0 + l_k
#                }
    
    #            Yterr_full <- crossprod(backsolve(uVt, backsolve(uVt, Yterr_full, transpose=TRUE)), Jt)
    
#                mu_tilde_t <- mu_tilde_t + Yterr_full[, sp_coords_1ix]
            } else {
                #TODO: Only applicable  
                Theta_t_l <- matrix(Theta_L[t,], nrow=p, ncol=S)[,sp_coords_1ix]
                Yterr <- crossprod(backsolve(uVt, backsolve(uVt, Yt_train - Ft_train %*% Theta_t_l, transpose=TRUE)), Jt)[,sp_coords_1ix]
                #TODO: Theta_t_l %*% Ftmat[t-2,] will need to be processed separately, especially for the big data case.
                mu_tilde_t <- t(matrix(t(Ftkmat[((t-2)*K+k):((t-1)*K),]), nrow=1) %*% Theta_t_l) + Yterr
            }
   
            #TODO: Do the below after all k's have had their run

	    if (k_ix0 == Stilde) {
                # construct Bt and bt
                Bt <- duplicate(Sigmatilde_eta, shallow=FALSE) #tau2s[t] * diag(Stilde)
	        Bt <- Bt + tau2s[t] * diag(Stilde)
#       	        .Call("axpy", Stilde*Stilde, tau2s[t], diag(Stilde), Bt)
#       	        .Call("axpy", Stilde, 1.0, rep(tau2s[t], Stilde), diag(Bt))
    
                Bt <- -tau2s[t]^2 * chol2inv(chol(Bt)) 
		Bt <- Bt + tau2s[t] * diag(Stilde)
#                .Call("axpy", Stilde*Stilde, tau2s[t], diag(Stilde), Bt)
#       	        .Call("axpy", Stilde, 1.0, rep(tau2s[t], Stilde), diag(Bt))

		#TODO: debug check
#		if (t == 3) {
#                    print("Bt Fnorm:")
#		    print(sum((tau2s[t] * diag(Stilde) - tau2s[t]^2 * chol2inv(chol(Sigmatilde_eta + tau2s[t] * diag(Stilde))) - Bt)^2))
#		}

    #            invBtSchur <- chol2inv(chol(Sigmatilde_eta + tau2s[t] * diag(Stilde)))
    #            Bt <- diag(Stilde) * tau2s[t] - tau2s[t]^2 * invBtSchur
    
                # Btbt is a placeholder name for bt.
    	        Btbt <- duplicate(mu_tilde_t, shallow=FALSE)

		Btbt <- backsolve(uSigmatilde_eta, backsolve(uSigmatilde_eta, Btbt, transpose=TRUE))
#                .Call("trsm", 'L', 'U', 'T', 'N', Stilde, 1, 1.0, uSigmatilde_eta, Stilde, Btbt, Stilde)
#                .Call("trsm", 'L', 'U', 'N', 'N', Stilde, 1, 1.0, uSigmatilde_eta, Stilde, Btbt, Stilde)
   
                Btbt <- Btbt + t(z[t,] - u[t+1,])/tau2s[t]
		Btbt <- t(Btbt)
#                .Call("axpy", Stilde, 1/tau2s[t], z[t,] - u[t+1,], Btbt)
   
    #            bt <- backsolve(uSigmatilde_eta, backsolve(uSigmatilde_eta, mu_tilde_t, transpose=TRUE)) +
    #                            t(z[t,] - u[t+1,])/tau2s[t]

                Btbt <- Bt %*% Btbt
      
                Btmat[t,] <- c(Bt)
                Btbtmat[t,] <- Btbt
        
                seed_obj_G <- update_seed_newcopy(seed_obj_G, seed_obj_G$seed_step)
                
                # Yt must be sampled from N(Bt * bt, Bt)
                Ytmat[t,] <- rmvnorm(1, Btbt, Bt)
	    }
        }
    }

    return(list(Y = Ytmat, #F = Ftkmat,
		Bt = Btmat, Btbt = Btbtmat, 
		seed_obj_TS = seed_obj_TS, seed_obj_G = seed_obj_G))
}


#TODO: substitute the below with Charlie Geyer's metrop.

#' This function takes in the set of new and old parameters, computes the ratio of probabilities at the Metropolis step, and then takes the probability to determine whether or not to take the Metropolis step. Outputs the Metropolis ratio computed, as well as whether or not the Metropolis step was taken.
#' @returns A list containing the Metropolis ratio and whether or not the Metropolis step is to be taken.

metrop_step <- function(delta, rho_old, rho_new, eta_old, eta_new,
		          logd_rho_eta, logJ_rho_eta,
			  rho_eta_2_Rspace,
			  utmat, tau2s, Uold, Unew,
			  Btbt_lm1, Bt_lm1, Btbt, Bt,
			  yold, ynew,
			  zdiffold, zdiffnew,
			  AR) {

    T <- nrow(zdiffold)
    Stilde <- ncol(zdiffold)

    #TODO: this particular formulation for rho_eta is wrong; we will need rho_eta_2_Rspace and Rspace_2_rho_eta. Rework that Jacobian.
    # Metropolis step for rho and eta begins here.
    # Begin the log-ratio with the log-Jacobian, sum(delta), and d_rho_eta.

    #TODO: What kind of normals for (rho_old, eta_old) and (rho_new, eta_new)???
    #TODO: Would we need to account for the updates in the standard deviations for rho_eta_2_Rspace?? The sd_list would update as well.
    fphi_old <- rho_eta_2_Rspace(rho_old, eta_old)
    fphi_new <- rho_eta_2_Rspace(rho_new, eta_new)
   
    # skip rho if it is the same. 
    l_rho <- 1#length(rho_old)
    if (rho_old == rho_new) {
        l_eta <- length(eta_old)
        fphi_old <- fphi_old[(1+l_rho):(l_rho + l_eta)]
        fphi_new <- fphi_new[(1+l_rho):(l_rho + l_eta)]
    }

    #TODO: relative to p(g(rho, eta)) instead? This is what Xiang did. Try that.
#    logr <- sum(dnorm(fphi_new, mean=0, sd=1, log=TRUE)) - 
#	    sum(dnorm(fphi_old, mean=0, sd=1, log=TRUE)) +
#	    logJ_rho_eta(rho_new, eta_new) - logJ_rho_eta(rho_old, eta_old)

    logr <- logd_rho_eta(rho_new, eta_new) - logd_rho_eta(rho_old, eta_old) -
	    logJ_rho_eta(rho_new, eta_new) + logJ_rho_eta(rho_old, eta_old)

#    logr <- sum(dnorm(fphi_new, mean=0, sd=pi/sqrt(3), log=TRUE)) - 
#	    sum(dnorm(fphi_old, mean=0, sd=pi/sqrt(3), log=TRUE)) +
#	    logJ_rho_eta(rho_new, eta_new) - logJ_rho_eta(rho_old, eta_old)

    #TODO: temporary measure with uniform default densities for (rho, eta). These are constants and can be ignored here, but the code should support more flexible implementations.
    # use uniform
#    logr <- logJ_rho_eta(rho_old, eta_old) - logJ_rho_eta(rho_new, eta_new)
#    logr <- logJ_rho_eta(rho_old, eta_old, delta) 
#    logr <- logJ_rho_eta(rho_old, eta_old, delta) + logd_rho_eta(rho_new, eta_new) - 
#		logd_rho_eta(rho_old, eta_old)

    for (t in 1:T) {
        utm1 <- matrix(data=unlist(utmat[t,]), nrow=1)
        ut <- matrix(data=unlist(utmat[t+1,]), nrow=1)

	logr <- logr + dmvnorm(ut, mean=utm1, sigma=tau2s[t+1]*Unew, log=TRUE) - 
                       dmvnorm(ut, mean=utm1, sigma=tau2s[t+1]*Uold, log=TRUE)

        #TODO: debug
#	    print(t)

        # NOTE: gen_ys initializes the first two time points from the PDE solution directly. Hence, these density terms will be skipped
        #TODO: Is this actually supposed to be the case???
        if (t >= AR) {
		dlogr <- if (!all(Bt[t,] == 0)) {
			          dmvnorm(matrix(data=ynew[t,], nrow=1),
				       mean=matrix(data=Btbt[t,], nrow=1), 
				       sigma=matrix(Bt[t,], nrow=Stilde, ncol=Stilde), log=TRUE)
                         } else {0} - 
			 if (!all(Bt_lm1[t,] == 0)) {
                               dmvnorm(matrix(data=yold[t,], nrow=1),
                                       mean=matrix(data=Btbt_lm1[t,], nrow=1),
                                       sigma=matrix(Bt_lm1[t,], nrow=Stilde, ncol=Stilde),log=TRUE)
                         } else {0}

		 #TODO: debug
#		 print(sprintf("logr pre-sum: %.3f", logr))

            # Note that y is deterministic if Bt is zero. So just add 1 or minus 1 if zero
            logr <- logr + if (!all(Bt[t,] == 0)) {
			          dmvnorm(matrix(data=ynew[t,], nrow=1),
				       mean=matrix(data=Btbt[t,], nrow=1), 
				       sigma=matrix(Bt[t,], nrow=Stilde, ncol=Stilde), log=TRUE)
                           } else {0}

            logr <- logr - if (!all(Bt_lm1[t,] == 0)) {
                               dmvnorm(matrix(data=yold[t,], nrow=1),
                                       mean=matrix(data=Btbt_lm1[t,], nrow=1),
                                       sigma=matrix(Bt_lm1[t,], nrow=Stilde, ncol=Stilde),log=TRUE)
                           } else {0}

	    #TODO: debug
#	    print("y* and y:")
#	    print(ynew[t,])
#	    print(yold[t,])
#	    if (!all(Bt[t-AR,] == 0)) {
#                print(dmvnorm(matrix(data=ynew[t,], nrow=1),
#				       mean=matrix(data=Btbt[t-AR,], nrow=1), 
#				       sigma=matrix(Bt[t-AR,], nrow=Stilde, ncol=Stilde), log=TRUE))
#	    }
#	    if (!all(Bt_lm1[t-AR,] == 0)) {
#                print(dmvnorm(matrix(data=yold[t,], nrow=1),
#                                       mean=matrix(data=Btbt_lm1[t-AR,], nrow=1),
#                                       sigma=matrix(Bt_lm1[t-AR,], nrow=Stilde, ncol=Stilde),log=TRUE))
#	    }

            #TODO: debug
#	    if ((t <= 3) | (t >= T - AR)) {
#            print(sprintf("B_{%d} all zero?", t-AR+1))
#            print(all(Bt_lm1[t-AR+1,] == 0))
#            print(all(Bt[t-AR+1,] == 0))
#	    }

#	    print(dlogr)
#            print(logr)

        }

#	if ((t == AR + 1) | (t == T)) {
#            print("zdiffnew:")
#            print(zdiffnew[t,])
#            print("zdiffold:")
#            print(zdiffold[t,])
#	}
#	print("ut:")
#	print(ut)
#	print("tau2s")
#	print(tau2s)

#	dlogr_sum <- 0
        for (i in 1:Stilde) {
	    #TODO: debug. dlogr term here is just to store the difference.
            dlogr <- dnorm(zdiffnew[t,i], mean=ut[i], sd=sqrt(tau2s[t+1]), log=TRUE) - 
                     dnorm(zdiffold[t,i], mean=ut[i], sd=sqrt(tau2s[t+1]), log=TRUE)
#	    dlogr_sum <- dlogr_sum + dlogr
            logr <- logr + dlogr

	    #TODO: debug. Print only initial and final times.
#	    if ((t == AR + 1) | (t == T)) {
#   	 	   print(sprintf("Scoord %d", i))
#	    print(c(zdiffnew[t,i], zdiffold[t,i]))
#	    print(c(ut[i], sqrt(tau2s[t+1])))
#	    	print(c(dnorm(zdiffnew[t,i], mean=ut[i], sd=sqrt(tau2s[t+1]), log=TRUE),
#                        dnorm(zdiffold[t,i], mean=ut[i], sd=sqrt(tau2s[t+1]), log=TRUE)))
#	    }
#	    print(sprintf("t = %d, spatial index %d, dlogr for z: %.7f", t, i, dlogr))
	}

	#TODO: debug
#	print(dlogr_sum)
#	if ((t == AR + 1) | (t == T)) print(logr)
    }

    metp <- exp(min(logr, 0))

    #TODO: debug
    print(logr)

    # set rho_eta <- rho_eta_new with probability metp (i.e. gen a uniform)
    p_sample <- runif(1)
    return(list(metp = metp, updated = p_sample <= metp))
}


gibbs_1iter <- function(z, y0,
			Y_files, F_files,
			T, N, S, Stilde, p, K, AR,
			n0z, d0z, b,
			uts, m0z, M0z,
			rho, rho_true, eta, eta_list,
	                logd_rho_eta_ud, rho_eta_2_Rspace, Rspace_2_rho_eta, logJ_rho_eta,
			sp_coords_1ix, sp_coords_1ix_lks, sp_dists,
		        pde_fun, grid_dims, n_PDEstates, PDEini, PDEsol_2_Y, gen_F,
			Btbt_lm1, Bt_lm1,
                        Sigma_block,
			uVt, beta,
			seed_obj_TS, seed_obj_G,
			nT = NULL, DT = NULL, dT = NULL, R = NULL,
			Theta_L = NULL,
			st_files = NULL, St_files = NULL,
			eps3 = 5e-2, Upsilon = diag(length(eta) + 1),#length(rho_true)),
			out.format = "csv",
			fst.head = "V",
			save.dir = NULL) {

    zdiffmat <- z - y0
    
    #TODO: rho is supposed to be a scalar in the calibration case we consider. Should I generalize to a vector rho according to the paper?
    U <- exp(-rho * sp_dists)

    seed_obj_G <- update_seed_newcopy(seed_obj_G, seed_obj_G$seed_step)
    tau2s <- sample_tau2s(n0z, d0z, b, uts, zdiffmat, U)

    seed_obj_G <- update_seed_newcopy(seed_obj_G, seed_obj_G$seed_step)
    uts <- FFBS_calibration(zdiffmat, m0z, M0z * tau2s[1], tau2s[2:(T+1)], U, save.dir = save.dir)

    y0_resample <- if (analytic_Y) {
                       gen_ys_analytic(eta, eta_list, T, AR, sp_coords_1ix, sp_coords_1ix_lks, PDEini, pde_fun, n_PDEstates, grid_dims,
                                       PDEsol_2_Y, gen_F, seed_obj_TS = seed_obj_TS, seed_obj_G = seed_obj_G)
                   } else {
                       gen_ys(eta, eta_list, sp_coords_1ix, sp_coords_1ix_lks, PDEini, pde_fun, n_PDEstates, grid_dims,
                              PDEsol_2_Y, gen_F,
                              tau2s[2:(T+1)], z, uts,
                              Y_files, F_files, Sigma_block,
                              uVt, beta,
                              N, S, Stilde, p, K, AR,
                              Theta_L = Theta_L, 
                              st_files = st_files, St_files = St_files,
                              seed_obj_TS = seed_obj_TS,
                              seed_obj_G = seed_obj_G,
                              out.format = out.format)
                   }

    y0 <- y0_resample$Y
    Btbt_lm1 <- y0_resample$Btbt
    Bt_lm1 <- y0_resample$Bt

    seed_obj_G <- y0_resample$seed_obj_G
    seed_obj_TS <- y0_resample$seed_obj_TS

    #TODO: debug time
    t1 <- Sys.time()
    print(sprintf("ut runtime: %.3f", t1 - t0))
    t0 <- t1

    #TODO: Take in a function to translate eta values to infinite support (for all dimensions) and back before doing the update.

    l_rho <- 1#length(rho)
    l_eta <- length(eta)

    phi <- rho_eta_2_Rspace(rho, eta)

    # eta* proposal generation is contained in a while loop to ensure we get one that results in a set of working solutions y_1:T(eta*).
    yf <- NaN
    while (any(is.nan(unlist(yf)))) {
        # update the proposal here.
        seed_obj_G <- update_seed_newcopy(seed_obj_G, seed_obj_G$seed_step)
        delta <- rmvnorm(1, rep(0, l_rho + l_eta), eps3 * Upsilon)
#		rnorm(l_rho + l_eta, 0, 1) 
        if (!is.null(rho_true)) delta[1:l_rho] <- 0
  
        rho_eta_new <- Rspace_2_rho_eta(phi + delta, l_rho, l_eta)
    
        rho_new <- if(!is.null(rho_true)) rho_true else rho_eta_new[1:l_rho]
        eta_new <- rho_eta_new[(l_rho+1):(l_rho + l_eta)]
    
        Unew <- exp(-rho_new * sp_dists)
       
        names(eta_new) <- names(eta)
    
        #TODO: debug
        print(eta)
        print(eta_new)
    
        # iteratively generate the y's and f's.
        yf <- if (analytic_Y) {
                  gen_ys_analytic(eta_new, eta_list, T, sp_coords_1ix, sp_coords_1ix_lks, PDEini, pde_fun, n_PDEstates, grid_dims,
                     PDEsol_2_Y, gen_F, seed_obj_TS = seed_obj_TS, seed_obj_G = seed_obj_G)
              } else {
                  gen_ys(eta_new, eta_list, sp_coords_1ix, sp_coords_1ix_lks, PDEini, pde_fun, n_PDEstates, grid_dims,
                     PDEsol_2_Y, gen_F,
    	             tau2s[2:(T+1)], z, uts,
                     Y_files, F_files, Sigma_block,
    		     uVt, beta,
                     N, S, Stilde, p, K, AR,
                     Theta_L = Theta_L, 
                     st_files = st_files, St_files = St_files,
                     seed_obj_TS = seed_obj_TS,
                     seed_obj_G = seed_obj_G,
                     out.format = out.format)
              }
    
        #TODO: debug time
        t1 <- Sys.time()
        print(sprintf("y runtime: %.3f", t1 - t0))
        t0 <- t1
    }

    zdiffmat <- z - y0
    zdiffmat_star <- z - yf$Y

    # metropolis step here.
    seed_obj_G <- yf$seed_obj_G
    seed_obj_TS <- yf$seed_obj_TS
    seed_obj_G <- update_seed_newcopy(seed_obj_G, seed_obj_G$seed_step)
    M_result <- metrop_step(delta, rho, rho_new, eta, eta_new,
                            logd_rho_eta_ud, logJ_rho_eta,
                            rho_eta_2_Rspace,
                            uts, tau2s, U, Unew,
                            Btbt_lm1, Bt_lm1, yf$Btbt, yf$Bt,
                            y0, yf$Y,
                            zdiffmat, zdiffmat_star,
                            AR)
    
    #TODO: debug
    print(M_result$metp)
    print(M_result$update)

    if (M_result$update) {
        names(eta_new) <- names(eta)
        # return new results
        return(list(rho = rho_new,
		    eta = eta_new,
		    tau2s = tau2s,
		    uts = uts,
		    y = yf$Y,
		    Btbt = yf$Btbt,
		    Bt = yf$Bt,
		    delta = delta,
		    metp = M_result$metp,
		    update = M_result$update,
		    seed_obj_TS = seed_obj_TS,
		    seed_obj_G = seed_obj_G))
    } else {
        # return current results
        return(list(rho = rho,
		    eta = eta,
		    tau2s = tau2s,
		    uts = uts,
		    y = y0,
		    Btbt = Btbt_lm1,
		    Bt = Bt_lm1,
		    delta = NULL,
		    metp = M_result$metp,
		    update = M_result$update,
		    seed_obj_TS = seed_obj_TS,
		    seed_obj_G = seed_obj_G))

    }	    
}

#TODO: Turn Vt back to a function when we can enable it to vary over time for our applications.

#' The calibrator for the code.
#' @param z A T x Stilde matrix containing the field data to be used for the calibrator. Each row of z consists of one instance of zt, measured at Stilde locations.
#' @param Y_files A list containing the file names which contain the Yt's used for emulation. Each row of one of the files contains the flattened (originally N x (S/K)) Yt matrix, where K is the number of episodes per season. For this problem, there should only be one row per file; if there are more, only the row corresponding to eta0 will be used.
#' @param F_files A list containing the file names which contain the Ft's used for emulation. Each row of one of the files contains the flattened (originally N x (p/K)) Ft matrix. For this problem, there should only be one row per file; if there are more, only the first row will be used.
#' @param N The number of different sets of eta used to train the emulator.
#' @param p The number of covariates used to train or evaluate the emulator.
#' @param n0z The shape prior for the tau^2_0 inverse-Gamma distribution.
#' @param d0z The scale prior for the tau^2_0 inverse-Gamma distribution.
#' @param b The discount factor for the tau^2_t inverse-Gamma series of distributions.
#' @param m0z The mean prior for the bias term u_0.
#' @param M0z The variance prior for the bias term u_0, sans tau^2_0.
#' @param rho0 The initial value for the parameter rho, which controls the bias covariance term.
#' @param eta0 The initial value for the parameter eta, which controls the mean terms for y_t. eta0 is a named vector whose names correspond to the parameters of the PDE. e.g. names(eta) == c("eta1", "eta2", "alpha1", "alpha2", "alpha3")
#' @param eta_list An N x d matrix containing the parameter values to train the emulator. Required to construct V_t and J_t.
#' @param logd_rho_eta The log of the joint prior density of rho and eta. logd_rho_eta must take in rho and eta as arguments in that order. 
#' @param rho_eta_2_Rspace The function to transform values of rho and eta to (1 + l(eta))-dimensional space R^{1 + l(eta)}, where l(eta) denotes the length of eta. Takes in rho and eta in that order.
#' @param Rspace_2_rho_eta The function needed to transform values of R^{1 +l(eta)} back to rho and eta space. Takes in iota = c(rho, eta), l_rho the length of rho, and l_eta the length of eta, in that order.
#' @param logJ_rho_eta A function computing the log of the Jacobian between (rho, eta) and its update (rho*, eta*). Takes in (the pre-updated values) rho, eta, and the random update delta. delta is randomly generated by the calibrator, and is as long as the combined length of rho and eta.
#' @param Vt The correlation matrix between the rows of Y_t at the training stage for emulation.
#' @param beta The scalar exponent in Vt sans the distance matrix. Used to construct Jt.
#' @param sp_coords An Stilde x 2 matrix of spatial coordinates used to denote the coordinates of z.
#' @param pde_fun The PDE function to be passed into the PDE solver to generate synthetic data.
#' @param grid_dims A 2-element vector containing the dimensions of the full grid of the PDE function. Note that the product of the entries of grid_dims must be >= Stilde.
#' @param n_PDEstates The number of states in the PDE to keep track of.
#' @param PDEini A grid_dim[1] x (n_PDEstates * grid_dim[2]) matrix containing the initial spatial configurations of each state of the PDE.
#' @param PDEsol_2_Y A function that takes the time index and the raw solutions of the PDE and outputs them to the desired value of Y.
#' @param gen_F A function that takes in the time index, parameter values eta, and Y values to produce a value of F_t as a set of arguments in that order. This function will be used to produce F_t-tilde, corresponding to Y(eta-tilde).
#' @param y0 A T x Stilde initial matrix of y_t at the Stilde spatial points. The rows of y0 contains one of y_t from t = 1 to T. If not specified, one will be chosen from an index closely matching the initial eta0, and then constructed from the Y_t used to train the emulator. Default: NULL.
#' @param Btmat0 The (T - 2) x Stilde^2 initial values of B_t to be taken. Each row contains a flattened (originally Stilde x Stilde) B_t for each value of y_t sampled beyond the second point for the AR2 scenario. Required for non-null y0. If no Btmat0 is specified, then a warning will be printed and calibration will take place with a y0 selected from the rows of Y_files. Default: NULL
#' @param Theta_files A list containing the file names which contain the Theta_t samples generated through emulation. Each row contains one flattened (originally p x (S/K)) sample of the distribution for Theta_t. Default: NULL.
#' @param Sigma_file The file containing the L upper-triangular entries of the bncol x bncol Sigma matrices generated by the emulator. Relevant where the model under consideration is the Matrix-Normal Inverse-Wishart. Default: NULL.
#' @param sigma2_file The file containing the L estimates of sigma^2 generated by the emulator. Relevant where the model under consideration is the Matrix-Normal Inverse-Gamma (MNIG) and sampling sigma2 is not to be done in the calibrator. Default: NULL.
#' @param R The distance scale matrix per episode of the episode-season model. Relevant where the model being considered is MNIG. Default: NULL.
#' @param seed_list_TS A list of random seeds to run. Each seed will control for the next quantity to be randomized. This set of seeds concerns the L samples of Sigma or sigma2 and the L samples of each Theta_tk. The number of such seeds needed is L * T * K if the Sigma's or sigma^2's are saved directly from backwards sampling or L * T * K + L if Sigma or sigma^2 are to be sampled in the calibrator itself. Can be passed as a full set of integers, or as a vector of 1 or 2 integers. If as a vector of 2 integers (a,b), the seed will start at a and be incremented by a step size of b for every instant where Sigma, sigma^2, or Theta_tk needs to be drawn. If as a single integer, the step size will be taken as 1.
#' @param seed_list_G A list of random seeds to run. Each seed will control for the next quantity to be randomized specific to the Gibbs sampler, first the batches of tau^2_t, then u_t, then y_t, and finally rho, eta, and the Metropolis step. The number of such seeds needed is L * (T * K + 3). Can be passed as a full set of integers, or as a vector of 1 or 2 integers. If as a vector of 2 integers (a,b), the seed will start at a and be incremented by a step sizeof b for every instant where one of the quantities tau^2_t, u_t, the individual y_tk's, or rho, eta, and the Metropolis step, need to be sampled. If as a single integer, the step size will be taken as 1. seed_list_G exists to control the randomization of the Gibbs quantities, especially when seed_list_TS has been specified.
#' @param nT The shape parameter for the inverse-Wishart or inverse-Gamma at the end of the FFBS step for the emulation with the MNIW model or the MNIG model respectively. Required for direct sampling of Sigma here. Default: NULL.
#' @param DT The scale matrix parameter for the inverse-Wishart at the end of the FFBS step for each episode at the end of emulation with the MNIW model. Required for direct sampling of Sigma. Default: NULL.
#' @param dT The scale parameter for the inverse-Gamma at the end of the FFBS step for the emulation with the MNIG model. Required for direct sampling of sigma^2. Default: NULL. 
#' @param st_files The array of files containing the T posterior means generated from smoothing. Relevant where L samples each of ThetaSt needs to be generated in the calibrator. Default: NULL.
#' @param St_files The array of files containing the T posterior row-covariance matrices generated from smoothing. Relevant where L samples each of Theta_t needs to be generated in the calibrator. Default: NULL.
#' @param AR The autoregression number used in the data and model. Required to adjust the calibrator to a certain offset. Default: 2.
#' @param L The number of Gibbs iterations to take. Defaults: 1000.
#' @param L_burn The number of Gibbs iterations to serve as burn-in. Default: 5000.
#' @param eps3 The size of the update step to take for rho and eta at each Gibbs step. Defaults to 0.05.
#' @param rho_true The "true value" of rho. If set, rho will be taken as known and will not be estimated. Default: NULL.
#' @param out.head The file output header. Defaults to './calib'.
#' @param save.dir The directory to save processing output, such as the moments of the FFBS at each Gibbs iteration. Defaults to NULL.
#' @param out.format The output format of the files for the calibrator. Supports both 'csv' and 'fst'. Defaults to 'csv'.
#' @param fst.head The header of the name for the columns of both Y_t and F_t. Required to parse the columns from these files. Defaults to 'V' (i.e. the columns of both Y_t and F_t are labeled V1, V2, V3,...)
#' @param max_iter_size The max number of iterations to go through before outputting output. Required because the results of all files may not be able to be stored in their entirety in memory. Will only store the results past burn-in. Default: 1000.
#' @param verbose Notify user for every max_iter_sie iterations of the Gibbs sampler. Default: FALSE.

#TODO: read in K as an option?

#TODO: add options to start the calibrator at certain iterations, e.g. ut0, index0 for file printout.

#TODO: Include Upsilon as an argument for the Metropolis proposal matrix.

calibrator <- function(z, Y_files, F_files, N, p,
		       n0z, d0z, b, m0z, M0z, rho0, eta0, eta_list,
		       logd_rho_eta, rho_eta_2_Rspace, Rspace_2_rho_eta, logJ_rho_eta,
		       Vt, beta,
		       sp_coords,
		       pde_fun, grid_dims, n_PDEstates, PDEini, PDEsol_2_Y, gen_F,
		       y0 = NULL, Btmat0 = NULL,
		       Theta_files = NULL,
		       Sigma_file = NULL, sigma2_file = NULL, R = NULL,
		       seed_list_TS = NULL, seed_list_G = NULL,
		       nT = NULL, DT = NULL, dT = NULL,
		       st_files = NULL, St_files = NULL, AR=2,
		       L = 1000, L_burn = 1000, eps3 = 5e-2,
		       Upsilon = diag(length(eta0) + 1),
		       rho_true = NULL, 
#		       gen_new_samples = FALSE,
		       out.head = "./calib",
		       save.dir = NULL,
		       out.format = "csv",
		       fst.head = "V",
                       max_iter_size = 1000,
		       verbose = FALSE) {

    #TODO: Adjust code to take in NULL sp_coords and grid_dims
    is_spatial <- !(is.null(sp_coords) & is.null(grid_dims))

    # Boolean controllers for the behavior of further code.
    sample_Sigma_here <- is.null(Sigma_file) & (!is.null(nT) & !is.null(DT))
    sample_sigma2_here <- is.null(sigma2_file) & (!is.null(nT) & !is.null(dT))

    is_IW <- !is.null(Sigma_file) | sample_Sigma_here

    is_IG <- (!(is.null(sigma2_file)) | sample_sigma2_here) & !is.null(R)

    #TODO: Rather than write them together (spaghetti code???) write separate functions that permit sampling with moments and preset seeds (for reproducibility) on the spot??? 
    #TODO: split along IW and IG? Or among sample vs no-sample? Better option, write a function to handle both cases

    sample_Theta_here <- (!is.null(st_files) & !is.null(St_files)) & is.null(Theta_files) 
    if (!(is_IW | is_IG)) {
        stop("Either Sigma_file or sigma2_file and R must be specified!")
    }
    if (is_IW & is_IG) {
        print("Warning! Files for both the IW and IG setup are detected. Specify either Sigma_file = NULL or sigma2_file = NULL and R = NULL to resolve this discrepancy. I will default to IW.")
        is_IG <- FALSE
    }
    
#    AR <- 2
    T <- nrow(z)
    Stilde <- ncol(z)
    z <- as.matrix(z) # force z to be a numeric matrix.

    if (is_spatial) {
        if ((grid_dims[1] * grid_dims[2] > 1) & (grid_dims[1] * grid_dims[2] < Stilde))
            stop("Full grid dimensions must be >= Stilde or the underlying ODE must not be spatial (i.e. the grid should be (1,1)).")
    }

    K <- if (!is.null(Theta_files)) {
             length(Theta_files)/T
         } else if (!is.null(st_files)) {
             if (out.format == "csv") {
                 nrow(read.csv(st_files[T])) - 1
	     } else if (out.format == "fst") {
                 nrow(read.fst(st_files[T])) - 1
	     }
         }

    # preprocess seeds here
    # seed list
    # index of the seed list
    seed_ix <- 1
    seed_now <- 0
    seed_step <- 1
    seedarr_explicit_TS <- is.numeric(seed_list_TS) & (length(seed_list_TS) > 2)

    if (!is.null(seed_list_TS)) {
        if (is.numeric(seed_list_TS)) {
            # length 2 or 3
            seed_len <- length(seed_list_TS)
            n_seeds_TS <- L * T * K + L * (sample_Sigma_here | sample_sigma2_here)        

	    if (seed_len <= 2) {
                seed_now <- seed_list_TS[1]
                seed_step <- if (seed_len == 1) 1 else seed_list_TS[2]
	    } else {
                seed_step <- 1
	        if (length(seed_list_TS) < n_seeds_TS) stop(sprintf("The explicit list of seeds you have provided for Sigma or sigma2, and Theta_tk, is too short; there must be %d integers if you want to specify all of them.", n_seeds_TS))
	    }
        }
        #TODO: include implementation for large seed list. Possibly a file.
    }

    # store seeds in this object to be manipulated later.
    seed_obj_TS <- make_seed_obj(seed_list_TS, seed_ix, seed_now, seed_step, seedarr_explicit_TS)

    # iniialize seed_list_G
    seedarr_explicit_G <- is.numeric(seed_list_G) & (length(seed_list_G) > 2)
    
    if (!is.null(seed_list_G)) {
        if (is.numeric(seed_list_G)) {
            # length 2 or 3
            seed_len <- length(seed_list_G)
            n_seeds_G <- L * (T * K + 3)

	    if (seed_len <= 2) {
                seed_now <- seed_list_G[1]
                seed_step <- if (seed_len == 1) 1 else seed_list_G[2]
	    } else {
                seed_step <- 1
	        if (length(seed_list_G) < n_seeds_G) stop(sprintf("The explicit list of seeds you have provided for the Gibbs quantities is too short; there must be %d integers if you want to specify all of them.", n_seeds_G))
	    }
        }
        #TODO: include implementation for large seed list. Possibly a file.
    }

    seed_obj_G <- make_seed_obj(seed_list_G, seed_ix, seed_now, seed_step, seedarr_explicit_G)

    # initialize and pre-allocate arrays.
    l_rho <- if (is.null(rho_true)) 1 else length(rho0)
    l_eta <- ncol(eta_list)
    
    rhos <- matrix(nrow=L+1, ncol=l_rho, data=-1)
    etas <- matrix(nrow=L+1, ncol=l_eta, data=-1)

    # treat rho0 as rho_true if it is set.
    if (!is.null(rho_true)) rho0 <- rho_true
#    eta0 <- unlist(eta_list[eta0_ix,])

    rhos[1,] <- rho0
    etas[1,] <- eta0

    sp_colnames <- paste0(fst.head, seq(1, Stilde))

    #TODO: We need burn in array for L_burn!!
    nrow_out <- min(L_burn, L, max_iter_size)
    if (L_burn == 0) nrow_out <- min(L, max_iter_size)

    uts <- as.data.frame(matrix(nrow = (nrow_out + 1) * (T + 1), ncol=Stilde, data=0))
    names(uts) <- sp_colnames
    uts <- cbind(expand.grid(l = seq(1, (nrow_out + 1)), t = seq(0,T)), uts)


    #TODO: Remove Btbtmat and uptriBtmat and just go with the l - 1 iterate?
    # Init matrices Btbtmat and uBtmat for storage.
#    Btbtmat <- as.data.frame(matrix(nrow=(min(L, max_iter_size)+1) * (T-2), ncol=Stilde, data=0))
#    names(Btbtmat) <- sp_colnames
#    Btbtmat <- cbind(expand.grid(l =seq(0,min(L, max_iter_size)), t=seq(3,T)), Btbtmat)

#    # Store the upper triangular entries of Bt
#    uptriBtmat <- as.data.frame(matrix(nrow=(min(L, max_iter_size)+1)*(T-2), ncol=(Stilde * (Stilde + 1)/2), data=0))
#
#    uptri_sp_colnames <- paste0(fst.head, seq(1, (Stilde * (Stilde+1))/2))
#    names(uptriBtmat) <- uptri_sp_colnames
#    uptriBtmat <- cbind(expand.grid(l=seq(0,min(L, max_iter_size)), t=seq(3,T)), uptriBtmat)

    # Record and aggregate scale to record the update for logd_rho_eta
    logd_rho_eta_ud <- logd_rho_eta
    agg_scale <- 0

    # store metropolis probabilities.
    metps <- rep(-1, L)

    sp_coords_1ix <- if (is_spatial) (sp_coords$x - 1) * grid_dims[2] + sp_coords$y else seq(1,Stilde)

    S <- if (is_spatial) grid_dims[1] * grid_dims[2] else Stilde
    bnrow <- N
    bncol <- S/K
    fbncol <- p/K #TODO: I think we need the f_grid here, but we can leave this here for now owing to its use case.
    
    # What data structure should I use to store the sp_coords_1ix_k? A dict?
    #sp_coords_1ix_k_all <- dict()

    sp_coords_1ix_lks <- rep(0, K)

    for (k in 1:K) {
        sp_coords_1ix_k <- sp_coords_1ix[(sp_coords_1ix <= k * bncol) & 
                                         (sp_coords_1ix >= (k-1) * bncol + 1)] - 
                               bncol * (k-1)
        #sp_coords_1ix_k_all$set(k, sp_coords_1ix_k)

	sp_coords_1ix_lks[k] <- length(sp_coords_1ix_k)
    }

#    sp_coords_uptri_1ix <- rep(1, (Stilde * (Stilde + 1))/2)
#    for (i in 2:Stilde) {
#        sp_coords_uptri_1ix[i:((i*(i+1))/2)] <- seq((i-1)*Stilde+1, (i-1)*Stilde+i)
#    }
   
#    deff <- 0.4 * max(dists)
    uVt <- chol(Vt)

    Jt <- exp(-beta * sqrt(rdist(eta_list, matrix(eta0, nrow=1))))

    #TODO: Since this eta leads to a scale_factor of zero, just use it here.
    scale_factor <- backsolve(uVt, Jt, transpose=TRUE)
    scale_factor <- 1 - sum(scale_factor^2) #dot(scale_factor, scale_factor)

    # Deal with sigma2 and Sigma here.
    if (is_IW) {
        # Record the Stilde * (Stilde + 1)/2 upper-triangular coordinates of the Stilde x Stilde matrix.
        # Subset Sigmamat based on sp_coords_1ix, since we will be reading in individual Sigma that are upper-triangular.
        # From upper-triangular, subset only the columns and rows corresponding to sp_coords_1ix.

       #TODO: check the coordinates here. May not need it.
        sp_coords_1ix_Sigma2Stilde <- rep(sp_coords_1ix[1]^2, (Stilde * (Stilde + 1)/2))
        for (i in 2:Stilde) {
           for (j in 1:i) 
    	       sp_coords_1ix_Sigma2Stilde[(i*(i-1)/2) + j] <- 
    	           sp_coords_1ix[i] * (sp_coords_1ix[i]-1) / 2 + sp_coords_1ix[j]
        }

        # save the upper-triangular 1D indices from 1 to Stilde * Stilde and then subset them.
	if (!is.null(Sigma_file)) {
            Sigmamat <- if (out.format == "csv") {
                            read.csv(Sigma_file)
                        } else if (out.format == "fst") {
                            read.fst(Sigma_file)
                        }
            Sigmamat <- Sigmamat[, sp_coords_1ix_Sigma2Stilde]
	}
    } else if (is_IG) {
        if (!is.null(sigma2_file)) sigma2s <- read.csv(sigma2_file)
        uR <- chol(R)
    }

    #NOTE: that this is one instance of sampling Sigma_block here.
    if (sample_Sigma_here | sample_sigma2_here) seed_obj <- update_seed_newcopy(seed_obj_TS, seed_obj_TS$seed_step)

    # Read or sample Sigma or sigma2
    if (is_IW) {
        # initialize Bt_lm1 and Btbt_lm1 from existing eta_ix0, etc. See gen_ys()
        # Note that each Sigma is per episode. Subdivisions are required.
        if (sample_Sigma_here) {
            Sigma_block <- rinvwishart(1, nT, DT)[,,1]

	} else {
            Sigma_block <- tri.to.sym(Sigmamat[1,], S, lower=FALSE)
	}
        Sigma_eta <- scale_factor * Sigma_block
    } else if (is_IG) {
        if (sample_sigma2_here) {
            sigma2 <- 1/rgamma(1, shape=nT, rate=dT)
	} else {
            sigma2 <- sigma2s[l]
	}
        Sigma_block <- sigma2 * R
        Sigma_eta <- scale_factor * Sigma_block
    }

    sp_dists <- if (is_spatial) rdist(sp_coords, sp_coords) else rdist(sp_coords_1ix, sp_coords_1ix)
    U <- exp(-rho0 * sp_dists)

    # Set priors for tau2s and u_t.
    tau2s <- matrix(ncol=T+1, nrow=L+1, data=-1)
    tau2s[1,1] <- 1/rgamma(1, shape=n0z, rate=d0z)

    uts_prior <- matrix(nrow=T+1, ncol=Stilde, data=0)
    utm1 <- rmvnorm(1, m0z, M0z * tau2s[1,1]) #rep(0, Stilde)
    for (t in 1:T) {
        tau2s[1,t+1] <- 1/rgamma(1, shape=b^t * n0z, rate=b^t * d0z)

        uts_prior[t+1,] <- rmvnorm(1, utm1, tau2s[1,t+1] * U)
	utm1 <- uts_prior[t+1,]
    }

    # With priors for tau2s and u_t, get y from eta0. If eta0 exists in the training set, it will return just that.
    y_out <- gen_ys_prior(eta0, eta_list, sp_coords_1ix, sp_coords_1ix_lks, PDEini, pde_fun, n_PDEstates, grid_dims,
                          PDEsol_2_Y, gen_F, tau2s[1,2:(T+1)], z, uts_prior, Y_files, F_files, Sigma_block,
                          uVt, beta,
                          N, S, Stilde, p, K, AR=AR,
                          Theta_L = Theta_files,
                          st_files = st_files, St_files = St_files,
                          seed_obj_TS = seed_obj_TS, seed_obj_G = seed_obj_G,
                          out.format = out.format,
                          fst.head = fst.head)

    if (any(is.nan(unlist(y_out)))) stop("eta0 yielded bad initial ODE solutions for prior data generation. Select a different eta0.")

    y0 <- y_out$Y
    Bt_lm1 <- y_out$Bt
    Btbt_lm1 <- y_out$Btbt
    seed_obj_TS <- y_out$seed_obj_TS
    seed_obj_G <- y_out$seed_obj_G

    # Generate Btbt and Bt based on y0 and store as l = 0.
    # Initialize Bt_lm1 and Btbt_lm1
    # Do everything here outside of the L-loop so that we don't need the if statement inside the L-loop
    # Generate Bt_lm1 and Btbt_lm1 here for the first iteration.

#    #TODO: How to adjust for arbitrary autoregression? Not just AR2?
#    for (t in 3:T) {
#        Btbt_lm1[t-2,] <- y0[t,] #z[t,] - uts[(uts$l == 1) & (uts$t == t), sp_colnames]
##        Bt_lm1[t-2,] <- c(Bt)
#    }

    #TODO: We can calculate eta0_ix.
    # Subset in calibrator.R to make the most general case.
    # Initialize yt here. 

    dists <- sqrt(rdist(eta_list, eta_list))

    # Define the function for updating logd_rho_eta
    l_rho <- length(rho0)
    l_eta <- length(eta0)
    l_phi <- as.integer(is.null(rho_true)) + length(eta0)
    delta_ixs <- seq((1+!is.null(rho_true)), (l_rho+l_eta))
    deltas <- matrix(nrow=(L_burn + L), ncol=l_phi,
		     data=0)
    accept_ix <- rep(0, L_burn + L)
#    update_logd_rho_eta <- function(logd_rho_eta, deltas) {
#        logd_rho_eta_ud <- function(rho, eta) {
#            logd <- logd_rho_eta(rho, eta)
#            for (j in 1:nrow(deltas)) {
#                logd <- logd + logJ_rho_eta(rho, eta, deltas[j,])
#	    }
#	    return(logd)
#	}
#
#        return(logd_rho_eta_ud)
#    }
    
    rho <- rho0
    eta <- eta0

    n_file_iter <- ceiling(L/nrow_out)

    # Gibbs iterations start here
    for (l in 1:(L + L_burn)) {
        #TODO: debug
        print(sprintf("Iter %d", l))

        Theta_L <- if (!sample_Theta_here) {
                       get_Theta_L(Theta_files, l, p, Stilde, out.format = out.format)
                   } else { NULL }

        #NOTE: that this is one instance of sampling Sigma_block here.
        if (sample_Sigma_here | sample_sigma2_here) seed_obj <- update_seed_newcopy(seed_obj_TS, seed_step)
    
        # Read or sample Sigma or sigma2
        if (is_IW) {
            # initialize Bt_lm1 and Btbt_lm1 from existing eta_ix0, etc. See gen_ys()
            if (sample_Sigma_here) {
                Sigma_block <- rinvwishart(1, nT, DT)[,,1]
    
            } else {
                Sigma_block <- tri.to.sym(Sigmamat[l,], bncol, lower=FALSE)
            }
        } else if (is_IG) {
            if (sample_sigma2_here) {
                sigma2 <- 1/rgamma(1, shape=nT, rate=dT)
            } else {
                sigma2 <- sigma2s[l]
            }
            Sigma_block <- sigma2 * R
        }

	#TODO: Is there something wrong with this kind of indexing? Some sort of problem going on here.
        lm1_fix <- ((l-1) %% nrow_out) + 1

        gibbs_out <- gibbs_1iter(z, y0, Y_files, F_files,
	                         T, N, S, Stilde, p, K, AR,
	                         n0z, d0z, b, 
	                         as.matrix(uts[(uts$l == lm1_fix), sp_colnames]), m0z, M0z,
	                         rho, rho_true, eta, eta_list,
	                         logd_rho_eta_ud, rho_eta_2_Rspace, Rspace_2_rho_eta, logJ_rho_eta,
	                         sp_coords_1ix, sp_coords_1ix_lks, sp_dists,
	                         pde_fun, grid_dims, n_PDEstates, PDEini, PDEsol_2_Y, gen_F,
			         Btbt_lm1, Bt_lm1,
				 Sigma_block,
			         uVt, beta,
	                         seed_obj_TS, seed_obj_G,
	                         nT = nT, DT = DT, dT = dT, R = R,
	                         Theta_L = Theta_L,
	                         st_files = st_files, St_files = St_files,
				 out.format = out.format,
				 fst.head = fst.head,
				 eps3 = eps3, Upsilon = Upsilon,
				 save.dir = save.dir)
	
        # Exterior sampling prior to Gibbs.
        rho <- gibbs_out$rho
	eta <- gibbs_out$eta

	if (!is.null(gibbs_out$delta)) {
            #agg_scale <- agg_scale + gibbs_out$delta #sum(gibbs_out$delta)
            deltas[l,] <- gibbs_out$delta[delta_ixs]
            accept_ix[l] <- l
#            logd_rho_eta_ud <- update_logd_rho_eta(logd_rho_eta,
#						   deltas[accept_ix[accept_ix > 0], delta_ixs])
	}

        uts[uts$l == lm1_fix, sp_colnames] <- gibbs_out$uts
        y0 <- gibbs_out$y
	seed_obj_TS <- gibbs_out$seed_obj_TS
	seed_obj_G <- gibbs_out$seed_obj_G

	# Save the lth iterate.
	if ((l > L_burn)) {
	    tau2s[l - L_burn +1,] <- gibbs_out$tau2s
            rhos[l - L_burn,] <- rho
            etas[l - L_burn,] <- eta
	    metps[l - L_burn] <- gibbs_out$metp
	}

	# progress update
	if (verbose & (l == L_burn)) print(sprintf("Burn-in finished with %d iterations.", L_burn))

	if ((l > L_burn) & ((l - L_burn) %% nrow_out == 0)) {
            if (verbose) print(sprintf("Iteration %d finished.", l - L_burn))

	    group_no <- (l - L_burn) / nrow_out

	    # output files here
            if (out.format == "csv") {
#                write.csv(paste0(out.head, sprintf("_Bt_uptri_%d.csv", group_no)), uptriBtmat)
#                write.csv(paste0(out.head, sprintf("_Btbt_%d.csv", group_no)), Btbtmat)
                write.csv(paste0(out.head, sprintf("_ut_%d.csv", group_no)), uts)
	    } else if (out.format == "fst") {
#                write.fst(uptriBtmat, paste0(out.head, sprintf("_Bt_uptri_%d.fst", group_no)))
#                write.fst(Btbtmat, paste0(out.head, sprintf("_Btbt_%d.fst", group_no)))
                write.fst(uts, paste0(out.head, sprintf("_ut_%d.fst", group_no)))
	    }
	}

	# Update uts after file has been output, because it will loop back to index 1 for uts.
        l_fix <- (l %% nrow_out) + 1
        uts[(uts$l == l_fix), sp_colnames] <- gibbs_out$uts
#	Btbtmat[Btbtmat$l == l_fix, sp_colnames] <- gibbs_out$Btbt
#        uptriBtmat[uptriBtmat$l == l_fix, uptri_sp_colnames] <- gibbs_out$Bt[,sp_coords_uptri_1ix]
	Btbt_lm1 <- gibbs_out$Btbt
	Bt_lm1 <- gibbs_out$Bt
    }

    # output all intermediate files. This includes rhos, etas, Btmat, Btbtmat, uts, tau2s, metps
    if (out.format == "csv") {
       write.csv(paste0(out.head, "_rho.csv"), as.data.frame(rhos))
       write.csv(paste0(out.head, "_eta.csv"), as.data.frame(etas))
       write.csv(paste0(out.head, "_tau2t.csv"), as.data.frame(tau2s))
       write.csv(paste0(out.head, "_metprobs.csv"), as.data.frame(metps))
       write.csv(paste0(out.head, "_metaccept.csv"), as.data.frame(accept_ix))
    } else if (out.format == "fst") {
       write.fst(as.data.frame(rhos), paste0(out.head, "_rho.fst"))
       write.fst(as.data.frame(etas), paste0(out.head, "_eta.fst"))
       write.fst(as.data.frame(tau2s), paste0(out.head, "_tau2t.fst"))
       write.fst(as.data.frame(metps), paste0(out.head, "_metprobs.fst"))
       write.fst(as.data.frame(accept_ix), paste0(out.head, "_metaccept.fst"))
    }
    
    # Return a list of the final estimates
    return(list(rho = rhos[L+1,], eta = etas[L+1,], tau2 = tau2s[L,],
		utmat = uts[uts$l == L,],
		Btbt = Btbt_lm1, Bt = Bt_lm1))#,
#		Btbtmat = Btbtmat[Btbtmat$l == L,], Btmat = uptriBtmat[uptriBtmat$l == L,]))
}

