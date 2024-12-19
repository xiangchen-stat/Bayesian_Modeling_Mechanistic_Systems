library(stringr)
library(fst)

BS_chain_1step <- function(Theta_tpl1, m, G, M, uAt){ #, omega){
    ht <- Theta_tpl1 - G %*% m
    ht <- backsolve(uAt, ht, transpose=TRUE)
    ht <- backsolve(uAt, ht)
    ht <- crossprod(G, ht)
    ht <- M %*% ht
    ht <- m + ht

    Ht <- G %*% M
    Ht <- backsolve(uAt, Ht, transpose=TRUE)
    Ht <- crossprod(Ht)
    Ht <- M - Ht
    Ht <- check.pds(Ht)
    return(list(h=ht, H=Ht))
}

#TODO: code up the small data backwards sampling.
BS_chain <- function(m_list, M_list, uAt_list, Gtmat=NULL) {
    if (is.null(Gtmat)) Gtmat <- matrix(c(diag(p)), nrow=1)
    change_Gt <- nrow(Gtmat) > 1
    Gtplus1 <- NULL
    if (!change_Gt) Gtplus1 <- matrix(Gtmat, nrow=p, ncol=p)
    

    out <- 0
}

#' Generates one set of backwards sampling coefficients for the big data setting given the list of files output by FF.bigdata.
#' @param m_file_list The mt_files in alphabetical order output by FF.bigdata. Each file is a dataframe of posterior means, with each row being the column-concatenated matrix m_t.
#' @param M_file_list The Mt_files in alphabetical order output by FF.bigdata. Each file is a dataframe of upper-triangular entries of posterior left-covariance matrices, with each row being the column-concatenated upper-triangular entries of M_t.
#' @param nT The degrees of freedom n_T for column covariance matrix for the FF. Generated from the FF for either the inverse-Wishart or inverse-Gamma case. 
#' @param DT The scale matrix D_T for the column covariance matrix for the FF. Generated from the FF for the invrese-Wishart case. 
#' @param dT The scale parameter d_T for the column covariance matrix for the FF. Generated from the FF for the inverse-Gamma case.
#' @param R The constant matrix R for the column covariance matrix for the FF. Pre-specified in the inverse Gamma case and should be the same as in the FF.
#' @param Gtmat The matrix of G_t used in FF.bigdata. Defaults to the identity matrix throughout if not specified.
#' @param Wtmat The matrix of W_t used in FF.bigdata. Defaults to the identity matrix throughout if not specified.
#' @param p The number of rows in each block m_tk. If not specified, the program will attempt to figure out the value of p from the file dimensions.
#' @param c The number of columns in each block m_tk. If not specified, the program will attempt to figure out the value of c from the file dimensions.
#' @param save_samples Whether or not to save the Theta_t and Sigma_T samples in files. Defaults to FALSE.
#' @param out.head The file path to output the smoothing coefficients. Defaults to the current directory with the "smooth_out" as the header.
#' @param out.format The format of the files to read and write. Supports 'csv' and 'fst'.
#' @param fst.head The header of the name for the columns of both Y_t and F_t. Required to parse the columns from these files. Defaults to 'V' (i.e. the columns of both Y_t and F_t are labeled V1, V2, V3,...)
#' @returns The list of files in which s_t and S_t are stored.
#' @export

BS_chain.bigdata <- function(m_file_list, M_file_list, nT, DT=NULL, dT=NULL, R=NULL, 
		       Gtmat = NULL, Wtmat = NULL,
		       p = NULL, c = NULL, save_samples=FALSE, 
		       out.head = "./BS_out",
		       out.format = "csv",
	  	       fst.head = "V") {
    # toggle options for IG vs IW like in the FF.
    iw_mode <- !is.null(DT)
    ig_mode <- !(is.null(dT) | is.null(R))

    if (!iw_mode & !ig_mode) stop("Either D0, or d0 and R, must be specified.")

    if (iw_mode & ig_mode) {
        warning("D0 and d0 and R are all specified! Defaulting to the inverse Wishart D0; d0 and R will not be considered.")
        ig_mode <- FALSE
    }

    if (iw_mode) {
        DTinv <- chol2inv(chol(DT))
        Sigma_T <- chol2inv(chol(matrix(rWishart(1, nT, DT),nrow=c)))
    } else if (ig_mode) {
        sigma2_T <- 1/rgamma(1, nT, rate=dT)
        Sigma_T <- sigma2_T * R
    }

    #TODO: try this in a small script in R to test its behavior
    # save Sigma_T. Use fst. We will use this from now on.
    if (save_samples) fst::write_fst(as.data.frame(Sigma_T), paste0(out.head, "_SigmaT.fst"))

    T <- length(m_file_list)

    if (length(M_file_list) != T) stop("The list of files for m_t and M_t must be the same length.")

    # Take in p and c as arguments. If is null, will try to figure out from the files, but this
    #	will take runtime.
    # Read in the columns of M_t. It has p(p+1)/2 columns.

    #TODO: need a way to insert p for fst if it is not given.
    if (is.null(p)) {
        p <- as.integer(system(sprintf("awk -F, '{print NF; exit}' %s", 
					  M_file_list[1]), intern=TRUE))
        # solve quadratic formula
        p <- as.integer((-1 + sqrt(1 + 8 * p))/2)
    }

    if (is.null(c)) c <- as.integer(as.integer(system(sprintf("awk -F, '{print NF; exit}' %s", 
					  m_file_list[1]), intern=TRUE))/p)

    # set default Gtmat and Wtmat if they aren't specified.
    if (is.null(Gtmat)) Gtmat <- matrix(data=c(diag(p)), nrow=1)
    if (is.null(Wtmat)) Wtmat <- matrix(data=c(diag(p)), nrow=1)
    
    # whether or not to transition Gt across time.
    change_Gt <- nrow(Gtmat) > 1
    Gtplus1 <- NULL
    if (!change_Gt) Gtplus1 <- matrix(data=Gtmat[1,], nrow=p, ncol=p)

    # whether or not to transition Wt across time.
    change_Wt <- nrow(Wtmat) > 1
    Wtplus1 <- NULL
    if (!change_Wt) Wtplus1 <- matrix(data=Wtmat[1,], nrow=p, ncol=p)

    ht_out_fname_list <- c()
    Ht_out_fname_list <- c()

    #NOTE: Files will be assumed to be listed in alphabetical order.
    # Get the values at the final block of the final time point
    # The elements are split into multiple lines here. Rejoin them and then split to get all the elements.
    #TODO: do not use system(). Can slow the code down dramatically.
    htplus1 <- system(sprintf("tail -n1 %s", m_file_list[T]), intern=TRUE)
    htplus1 <- paste0(htplus1, collapse='')
    htplus1 <- matrix(as.numeric(str_split_1(htplus1, ',')), nrow=p, ncol=c)

    Htplus1 <- system(sprintf("tail -n1 %s", M_file_list[T]), intern=TRUE)
    Htplus1 <- as.numeric(str_split_1(paste0(Htplus1, collapse=''), ','))
    Htplus1 <- tri.to.sym(Htplus1, p, lower=FALSE)

    Thetatplus1 <- rmn(1, htplus1, Htplus1, Sigma_T)[,,1]

    #TODO: turn this into a mono-loop which we then de-parse later???
    for (t in T:1) {
        # Read in the relevant m_file_list and M_file_list element
        mt_file <- if (out.format == "csv") {
	       	       read.csv(m_file_list[t])
                   } else if (out.format == "fst") {
	       	       read.fst(m_file_list[t])
	           }
        Mt_file <- if (out.format == "csv") {
	       	       read.csv(M_file_list[t])
	           } else if (out.format == "fst") {
	       	       read.fst(M_file_list[t])
	           }

	if (change_Gt) Gtplus1 <- matrix(data=Gtmat[t,], nrow=p, ncol=p)
        if (change_Wt) Wtplus1 <- matrix(data=Wtmat[t,], nrow=p, ncol=p)

	# K+1 is the length of the data frame.
	K <- nrow(mt_file) - 1

        ht_list <- matrix(nrow=K+1, ncol=p*c, data=0)
        Ht_list <- matrix(nrow=K+1, ncol=as.integer(p * (p+1)/2), data=0)

	ht_list[K+1,] <- c(htplus1)
	Ht_list[K+1,] <- Htplus1[upper.tri(Htplus1, diag=TRUE)]

        ht_out_fname <- paste0(out.head, sprintf("_ht_%d.%s", t, out.format))
        Ht_out_fname <- paste0(out.head, sprintf("_Ht_uptri_%d.%s", t, out.format))

        ht_out_fname_list <- c(ht_out_fname_list, ht_out_fname)
        Ht_out_fname_list <- c(Ht_out_fname_list, Ht_out_fname)

        for (k in K:1) {
	    mt <- matrix(unlist(mt_file[k,]), nrow=p, ncol=c)
	    # This is the upper triangle of Mt. Copy the entries to make it lower-triangular.
            Mt_upper <- unlist(Mt_file[k,])
	    Mt <- tri.to.sym(Mt_upper, p, lower=FALSE)
		    
            uMt <- chol(Mt)

            # compute uAt from GMG' + W.
	    uAtplus1 <- uMt %*% t(Gtplus1)
	    uAtplus1 <- chol(crossprod(uAtplus1) + Wtplus1)
	   
            # get Backwards Sampling samples.
            BS.param.list <- BS_chain_1step(Thetatplus1, mt, Gtplus1, Mt, uAtplus1)

	    htplus1 <- BS.param.list$h
	    Htplus1 <- BS.param.list$H

	    # sample Thetatplus1 here
            Thetatplus1 <- rmn(1, htplus1, Htplus1, Sigma_T)[,,1]

	    # save it to file
            if (save_samples) fst::write_fst(as.data.frame(Thetatplus1),
					paste0(out.head, sprintf("_Theta_%d_%d.fst", t, k)))

	    ht_list[k,] <- c(htplus1)
	    Ht_list[k,] <- Htplus1[upper.tri(Htplus1, diag=TRUE)]
	}

	#output files
	if (out.format == "csv") {
	    write.csv(ht_list, file=ht_out_fname, row.names=FALSE)
            write.csv(Ht_list, file=Ht_out_fname, row.names=FALSE)
	} else if (out.format == "fst") {
	    write.fst(ht_list, ht_out_fname, row.names)
            write.fst(Ht_list, Ht_out_fname, row.names)
	}
    }

    # return the list of smoothing files.
    out <- list(ht_files = ht_out_fname_list, Ht_files = Ht_out_fname_list)
    class(out) <- "BS.bigdata"
    return(out)
}


#' A single step of the backwards smoothing algorithm. The parameters are taken
#'   from the forward filter case. This computes the parameters of beta_t | D_T.
#'   beta_{t+1} | D_T is assumed to be normal with mean s and covariance matrix S.
#' 
#' @param s mean of beta_{t+1} | D_T.
#' @param S left-covariance matrix of beta_{t+1} | D_T
#' @param m mean of beta_{t} | D_t, computed in the FF step
#' @param G G_t, the beta transition matrix.
#' @param M left-covariance matrix of beta_{t} | D_t, computed in the FF step
#' @param uAt upper-triangular Cholesky factor of R_t computed in the FF step.
#' @returns The mean and left-covariance matrix of beta_{t} | D_T

smooth_1step <- function(s, S, m, G, M, uAt){ #, omega){
	snew <- s - G %*% m
	snew <- backsolve(uAt, snew, transpose = TRUE)
	snew <- backsolve(uAt, snew)
	snew <- crossprod(G, snew)
	snew <- M %*% snew
	snew <- m + snew

	Snew <- crossprod(uAt)
	Snew <- Snew - S
#	RinvGM <- chol2inv(uAt) %*% G %*% M
	AinvGM <- G %*% M
	AinvGM <- backsolve(uAt, AinvGM, transpose=TRUE)
	AinvGM <- backsolve(uAt, AinvGM)

	Snew <- Snew %*% AinvGM
	Snew <- crossprod(AinvGM, Snew)
	Snew <- M - Snew

	Snew <- check.pds(Snew)
	return(list(s = snew, S = Snew))
}

#' Acquire the smoothing coefficients after having run the Forward Filter. This function takes in the Forard Filter parameters and computes the smoothing coefficients.
#' @param m_list The list of means of beta_{t} | Y_{1:t}
#' @param M_list The list of left-covariance matrices of beta_{t} | Y_{1:t}
#' @param uAt_list The list of upper-Cholesky factors of A_{t}, the covariance matrix of beta_{t+1} | Y_{1:t}
#' @param Gt The transition matrix for beta_{t}: beta_{t} = G_t beta_{t-1} + Gamma_{t}
#' @param discount.factor The discount factor for the covariance matrix. Default: NULL.
#' @returns The list of mean and left-covariance matrices for beta_{t} | Y_{1:T}, for t = 0,...,T.
#' @export

smooth <- function(m_list, M_list, uAt_list, Gt = "I_p"){#, discount.factor = NULL){
    #if (!is.null(discount.factor) & !is.numeric(discount.factor)) {
    #    stop("Discount factor must be numeric if specified.")
    #}
    #if (is.null(discount.factor)) discount.factor <- 1.0

    T <- nrow(m_list) - 1
    p <- sqrt(ncol(M_list))
    m <- ncol(m_list) / p
    if (Gt == "I_p") Gt <- diag(p)

    beta_s_list <- matrix(rep(0, (T + 1) * p * m), nrow = T+1, ncol = p * m)
    M.T <- matrix(M_list[T+1,], p, p)
    s_list <- matrix(rep(0, (T+1) * p * m), nrow = T+1, ncol = p * m)
    S_list <- matrix(rep(0, (T+1) * p * p), nrow = T+1, ncol = p * p)
    s_list[T+1,] <- m_list[T+1,]
    S_list[T+1,] <- M_list[T+1,]

    for (i in T:1){
        mt <- matrix(m_list[i,], nrow=p)  #NOTE: m and M go from 0 to T.
        M <- matrix(M_list[i,],p,p) 
        uAt <- matrix(uAt_list[i,],p,p) # R and G go from 1 to T.
        	
#    	mstar <- time.param.list$GG %*% m_list[i,] 
#    	H <- matrix(S_list[i+1,],p,p)
#           	bs.results.1step <- BS.1step(s_list[i+1,], S, m_list[i,], time.param.list$GG,
#    				     M, uAt, sigma2.T)
        sS.t <- smooth_1step(matrix(s_list[i+1,], nrow=p), 
			     matrix(S_list[i+1,],nrow=p,ncol=p), 
			     mt, Gt, M, uAt)
        s_list[i,] <- c(sS.t$s)
        S_list[i,] <- c(sS.t$S)
    }

    return(list(s = s_list, S = S_list))
}


#' Generates smoothing coefficients for the big data setting given the list of files output by FF.bigdata.
#' @param m_file_list The mt_files in alphabetical order output by FF.bigdata. Each file is a dataframe of posterior means, with each row being the column-concatenated matrix m_t.
#' @param M_file_list The Mt_files in alphabetical order output by FF.bigdata. Each file is a dataframe of upper-triangular entries of posterior left-covariance matrices, with each row being the column-concatenated upper-triangular entries of M_t.
#' @param Gtmat The matrix of G_t used in FF.bigdata. Defaults to the identity matrix throughout if not specified.
#' @param Wtmat The matrix of W_t used in FF.bigdata. Defaults to the identity matrix throughout if not specified.
#' @param p The number of rows in each block m_tk. If not specified, the program will attempt to figure out the value of p from the file dimensions.
#' @param c The number of columns in each block m_tk. If not specified, the program will attempt to figure out the value of c from the file dimensions.
#' @param out.head The file path to output the smoothing coefficients. Defaults to the current directory with the "smooth_out" as the header.
#' @param out.format The format of the files to read and write. Supports 'csv' and 'fst'.
#' @param fst.head The header of the name for the columns of both Y_t and F_t. Required to parse the columns from these files. Defaults to 'V' (i.e. the columns of both Y_t and F_t are labeled V1, V2, V3,...)
#' @returns The list of files in which s_t and S_t are stored.
#' @export

#TODO: Note that uAtplus1 is getting reconstructed here. Load it in as files instead?

smooth.bigdata <- function(m_file_list, M_file_list, Gtmat = NULL, Wtmat = NULL,
			   p = NULL, c = NULL,
			   out.head = "./smooth_out",
			   out.format = "csv",
			   fst.head = "V") {
    T <- length(m_file_list)

    if (length(M_file_list) != T) stop("The list of files for m_t and M_t must be the same length.")

    # If p and c are not supplied as arguments, read in the columns of M_t. It has p(p+1)/2 columns.
    if (is.null(p) | is.null(c)) {
        if (out.format == "csv") {
            p <- as.integer(system(sprintf("awk -F, '{print NF; exit}' %s", 
    					  M_file_list[1]), intern=TRUE))
            p <- as.integer((-1 + sqrt(1 + 8 * p))/2)
            c <- as.integer(as.integer(system(sprintf("awk -F, '{print NF; exit}' %s", 
					  m_file_list[1]), intern=TRUE))/p)
        
            # solve quadratic formula
	} else if (out.format == "fst") {
            p <- ncol(read.fst(M_file_list[1], from=1, to=1))
            p <- as.integer((-1 + sqrt(1 + 8 * p))/2)
	    c <- as.integer(ncol(read.fst(m_file_list[1], from=1, to=1))/p)
	}
    }

    # set default Gtmat and Wtmat if they aren't specified.
    if (is.null(Gtmat)) Gtmat <- matrix(c(diag(p)), nrow=1)
    if (is.null(Wtmat)) Wtmat <- matrix(c(diag(p)), nrow=1)
    
    # whether or not to transition Gt across time.
    change_Gt <- nrow(Gtmat) > 1
    Gtplus1 <- NULL
    if (!change_Gt) Gtplus1 <- matrix(Gtmat[1,], nrow=p, ncol=p)

    # whether or not to transition Wt across time.
    change_Wt <- nrow(Wtmat) > 1
    Wtplus1 <- NULL
    if (!change_Wt) Wtplus1 <- matrix(Wtmat[1,], nrow=p, ncol=p)

    st_out_fname_list <- c()
    St_out_fname_list <- c()

    #NOTE: Files will be assumed to be listed in alphabetical order.
    # Get the values at the final block of the final time point
    # The elements are split into multiple lines here. Rejoin them and then split to get all the elements.
    #TODO: do not use system(). Can slow the code down dramatically.
    if (out.format == "csv") {
        stplus1 <- system(sprintf("tail -n1 %s", m_file_list[T]), intern=TRUE)
        stplus1 <- paste0(stplus1, collapse='')
        stplus1 <- matrix(as.numeric(str_split_1(stplus1, ',')), nrow=p, ncol=c)

        Stplus1 <- system(sprintf("tail -n1 %s", M_file_list[T]), intern=TRUE)
        Stplus1 <- as.numeric(str_split_1(paste0(Stplus1, collapse=''), ','))
    } else if (out.format == "fst") {
        # system won't work with FST. We'll need to read it in.
        stplus1 <- read.fst(m_file_list[T])
        stplus1 <- matrix(unlist(stplus1[nrow(stplus1),]), nrow=p, ncol=c)

	Stplus1 <- read.fst(M_file_list[T])
	Stplus1 <- unlist(Stplus1[nrow(Stplus1),])
    }
    Stplus1 <- tri.to.sym(Stplus1, p, lower=FALSE)


    #TODO: Stplus1 check if it is pds. If not, throw an error.

    for (t in T:1) {
        # Read in the relevant m_file_list and M_file_list element
        mt_file <- if (out.format == "csv") {
		       read.csv(m_file_list[t])
	           } else if (out.format == "fst") {
		       read.fst(m_file_list[t])
	           }
        Mt_file <- if (out.format == "csv") {
		       read.csv(M_file_list[t])
	           } else if (out.format == "fst") {
		       read.fst(M_file_list[t])
	           }

	if (change_Gt) Gtplus1 <- matrix(Gtmat[t,], nrow=p, ncol=p)
        if (change_Wt) Wtplus1 <- matrix(Wtmat[t,], nrow=p, ncol=p)

	# K+1 is the length of the data frame.
	K <- nrow(mt_file) - 1

        st_list <- matrix(nrow=K+1, ncol=p*c, data=0)
        St_list <- matrix(nrow=K+1, ncol=as.integer(p * (p+1)/2), data=0)

	st_list[K+1,] <- c(stplus1)
	St_list[K+1,] <- Stplus1[upper.tri(Stplus1, diag=TRUE)]

        st_out_fname <- paste0(out.head, sprintf("_st_%d.%s", t, out.format))
        St_out_fname <- paste0(out.head, sprintf("_St_uptri_%d.%s", t, out.format))

        st_out_fname_list <- c(st_out_fname_list, st_out_fname)
        St_out_fname_list <- c(St_out_fname_list, St_out_fname)

        for (k in K:1) {
	    mt <- matrix(unlist(mt_file[k,]), nrow=p, ncol=c)
	    # This is the upper triangle of Mt. Copy the entries to make it lower-triangular.
            Mt_upper <- unlist(Mt_file[k,])
	    Mt <- tri.to.sym(Mt_upper, p, lower=FALSE)
		    
            uMt <- chol(Mt)

            # compute uAt from GMG' + W.
	    uAtplus1 <- uMt %*% t(Gtplus1)
	    uAtplus1 <- chol(crossprod(uAtplus1) + Wtplus1)
	    
            smooth.param.list <- smooth_1step(stplus1, Stplus1, mt, Gtplus1, Mt, uAtplus1)

	    stplus1 <- smooth.param.list$s
	    Stplus1 <- smooth.param.list$S

	    st_list[k,] <- c(stplus1)
	    St_list[k,] <- Stplus1[upper.tri(Stplus1, diag=TRUE)]
	}

	#output files
	if (out.format == "csv") {
	    write.csv(st_list, file=st_out_fname, row.names=FALSE)
            write.csv(St_list, file=St_out_fname, row.names=FALSE)
	} else if (out.format == "fst") {
	    write.fst(as.data.frame(st_list), st_out_fname)
            write.fst(as.data.frame(St_list), St_out_fname)
	}
    }

    # return the list of smoothing files.
    out <- list(st_files = st_out_fname_list, St_files = St_out_fname_list)
    class(out) <- "BS.bigdata"
    return(out)
}

#' A single step of the backwards smoothing algorithm. The parameters are taken
#'   from the forward filter case. This computes the parameters of beta_t | D_T
#'   and returns L samples of Theta_t | D_T sampled from MN(s_t, S_t, Sigma).
#'   beta_{t+1} | D_T is assumed to be normal with mean s and covariance matrix S.
#' 
#' @param s mean of beta_{t+1} | D_T.
#' @param S left-covariance matrix of beta_{t+1} | D_T
#' @param m mean of beta_{t} | D_t, computed in the FF step
#' @param G G_t, the beta transition matrix.
#' @param M left-covariance matrix of beta_{t} | D_t, computed in the FF step
#' @param uAt upper-triangular Cholesky factor of R_t computed in the FF step.
#' @param L The number of samples to take for Theta_{t} | D_T. Defaults to 100.
#' @param Sigmas The column-covariance matrix, computed outside the function. Relevant where the column-covariance matrix is treated as inverse-Wishart.
#' @param sigma2 The scalar factor of the inverse-gamma-scale column covariance matrix.
#' @param R The scale matrix where the column-covariance matrix is treated is an inverse Gamma and scalar.
#' @param barebones Whether or not to store the seeds instead of the values of the generated backwards samples in the files. Default: FALSE.
#' @returns The mean and left-covariance matrix of beta_{t} | D_T, and L samples of Theta_t | D_T sampled from MN(s_t, S_t, Sigma).

BS_1step <- function(s, S, m, G, M, uAt, L, Sigmas = NULL, sigma2s = NULL, R = NULL, barebones = FALSE){
    smooth.out <- smooth_1step(s, S, m, G, M, uAt)

    p <- nrow(s)
    c <- ncol(s)

    is.IW <- !is.null(Sigmas)
    is.IG <- !(is.null(sigma2s) | is.null(R))

    if (!is.IW & !is.IG) stop("Either Sigmas, or sigma2s and R, must be specified in BS_1step.")

    if (barebones) {
        seed_ts <- rep(0, L)
        for (l in 1:L) {
	}
    } else {
        Theta_ts <- array(dim=c(p, c, L))
        for (l in 1:L) {
            Theta_ts[,,l] <- if (is.IW) {
    		             rmn(1, s, S, Sigmas[,,l])[,,1]
    	                 } else if (is.IG) {
    		             rmn(1, s, S, sigma2s[l] * R)[,,1]
    	                 }
        }
    }

    return(list(s = smooth.out$s, S = smooth.out$S, Theta_t = Theta_ts))
}


#' Generates backwards-sampling coefficients for the big data setting given the list of files output by FF.bigdata.
#' @param m_file_list The mt_files in alphabetical order output by FF.bigdata. Each file is a dataframe of posterior means, with each row being the column-concatenated matrix m_t.
#' @param M_file_list The Mt_files in alphabetical order output by FF.bigdata. Each file is a dataframe of upper-triangular entries of posterior left-covariance matrices, with each row being the column-concatenated upper-triangular entries of M_t.
#' @param Gtmat The matrix of G_t used in FF.bigdata. Defaults to the identity matrix throughout if not specified.
#' @param Wtmat The matrix of W_t used in FF.bigdata. Defaults to the identity matrix throughout if not specified.
#' @param nT The shape parameter for the inverse-Gamma or inverse-Wishart parameter produced by FF.bigdata.
#' @param DT The scale matrix for the inverse-Wishart produced by FF.bigdata. Relevant to the inverse-Wishart column covariance matrix case.
#' @param dT The scale parameter for the inverse-Gamma produced by FF.bigdata. Relevant to the inverse-Gamma column covariance matrix case.
#' @param p The number of rows in each block m_tk. If not specified, the program will attempt to figure out the value of p from the file dimensions.
#' @param c The number of columns in each block m_tk. If not specified, the program will attempt to figure out the value of c from the file dimensions.
#' @param L The number of backwards samples to draw. Default: 100.
#' @param seed_list A list of random seeds to run. Each seed will control for the next quantity to be randomized, first the L samples of Sigma or sigma2, then L samples of each Theta_tk. Can be passed as a csv or fst file, as a full set of L * (T * K + 1) integers, or as a sequence of 2 or 3 integers. If as a file, the file will contain at least L * (T * K + 1) integers to be used as seeds. If as a sequence of 3 integers (a,b,c), the seeds will be taken from the sequence (a,b) with a step size of c. If as a sequence of 2 integers, the step size c will be treated as 1.
#' @param seed_perchunk The number of seeds to read from file at once if seed_list is passed in as a file. Defaults: 10,000.
#' @param out.head The file path to output the smoothing coefficients. Defaults to the current directory with the "smooth_out" as the header.
#' @param out.format The format of the files to read and write. Supports 'csv' and 'fst'.
#' @param fst.head The header of the name for the columns of both Y_t and F_t. Required to parse the columns from these files. Defaults to 'V' (i.e. the columns of both Y_t and F_t are labeled V1, V2, V3,...)
#' @param verbose Whether or not to print out every 1000 samples of Theta_tk. If less than 1000 samples are sampled, then a message is printed out once all samples are finished. Default: FALSE.
#' @returns The list of files in which s_t and S_t are stored.
#' @export

BS.bigdata <- function(m_file_list, M_file_list, Gtmat = NULL, Wtmat = NULL,
			   nT, DT = NULL, dT = NULL, R = NULL,
			   p = NULL, c = NULL,
			   L = 100, 
			   seed_list = NULL,
			   seed_perchunk = 10000,
			   out.head = "./bs_out",
			   out.format = "csv",
			   fst.head = "V",
			   verbose=FALSE) {

    if (is.null(DT) & is.null(dT)) stop("Either DT or dT must be specified according to the covariance structure of the column covariance matrix Sigma.")
    if (is.null(DT) & (is.null(dT) | is.null(R))) stop("The MNIG structure requires both dT and R to be specified.")

    is.IW <- !is.null(DT)
    is.IG <- !is.null(dT)

    #TODO: copy/paste the code and edit it in
    T <- length(m_file_list)
    if (length(M_file_list) != T) stop("The list of files for m_t and M_t must be the same length.")

    # seed list
    # index of the seed list
    seed_ix <- 1
    seed_now <- 0
    seed_step <- 1

    seedarr_explicit <- is.numeric(seed_list) & (length(seed_list) > 3)
    mt_file <- if (out.format == "csv") {
	           read.csv(m_file_list[1])
               } else if (out.format == "fst") {
	           read.fst(m_file_list[1])
               }
    K <- nrow(mt_file) - 1
    rm(mt_file)
    n_seeds <- L * (T * K + 1)

    if (!is.null(seed_list)) {
        if (is.numeric(seed_list)) {
            # length 2 or 3 
            seed_len <- length(seed_list)
            if (((seed_len == 2) & (seed_list[2] - seed_list[1] + 1 < n_seeds)) | 
	       ((seed_len == 3) & ((seed_list[2] - seed_list[1])/seed_list[3] < n_seeds)))
		    stop(sprintf("seed_list numeric sequence needs to range over numbers that are at least L * (T * K + 1) = %d.", L * (T * K + 1)))
	    seed_now <- seed_list[1]
	    seed_step <- if (seed_len == 2) 1 else seed_list[3]
	}
    }

    # If p and c are not supplied as arguments, read in the columns of M_t. It has p(p+1)/2 columns.
    if (is.null(p) | is.null(c)) {
        if (out.format == "csv") {
            p <- as.integer(system(sprintf("awk -F, '{print NF; exit}' %s", 
    					  M_file_list[1]), intern=TRUE))
            p <- as.integer((-1 + sqrt(1 + 8 * p))/2)
            c <- as.integer(as.integer(system(sprintf("awk -F, '{print NF; exit}' %s", 
					  m_file_list[1]), intern=TRUE))/p)
        
            # solve quadratic formula
	} else if (out.format == "fst") {
            p <- ncol(read.fst(M_file_list[1], from=1, to=1))
            p <- as.integer((-1 + sqrt(1 + 8 * p))/2)
	    c <- as.integer(ncol(read.fst(m_file_list[1], from=1, to=1))/p)
	}
    }

    # set default Gtmat and Wtmat if they aren't specified.
    if (is.null(Gtmat)) Gtmat <- matrix(c(diag(p)), nrow=1)
    if (is.null(Wtmat)) Wtmat <- matrix(c(diag(p)), nrow=1)
    
    # whether or not to transition Gt across time.
    change_Gt <- nrow(Gtmat) > 1
    Gtplus1 <- NULL
    if (!change_Gt) Gtplus1 <- matrix(Gtmat[1,], nrow=p, ncol=p)

    # whether or not to transition Wt across time.
    change_Wt <- nrow(Wtmat) > 1
    Wtplus1 <- NULL
    if (!change_Wt) Wtplus1 <- matrix(Wtmat[1,], nrow=p, ncol=p)

    st_out_fname_list <- c()
    St_out_fname_list <- c()

    #NOTE: Files will be assumed to be listed in alphabetical order.
    # Get the values at the final block of the final time point
    # The elements are split into multiple lines here. Rejoin them and then split to get all the elements.
    #TODO: do not use system(). Can slow the code down dramatically.
    if (out.format == "csv") {
        stplus1 <- system(sprintf("tail -n1 %s", m_file_list[T]), intern=TRUE)
        stplus1 <- paste0(stplus1, collapse='')
        stplus1 <- matrix(as.numeric(str_split_1(stplus1, ',')), nrow=p, ncol=c)

        Stplus1 <- system(sprintf("tail -n1 %s", M_file_list[T]), intern=TRUE)
        Stplus1 <- as.numeric(str_split_1(paste0(Stplus1, collapse=''), ','))
    } else if (out.format == "fst") {
        # system won't work with FST. We'll need to read it in.
        stplus1 <- read.fst(m_file_list[T])
        stplus1 <- matrix(unlist(stplus1[nrow(stplus1),]), nrow=p, ncol=c)

	Stplus1 <- read.fst(M_file_list[T])
	Stplus1 <- unlist(Stplus1[nrow(Stplus1),])
    }
    Stplus1 <- tri.to.sym(Stplus1, p, lower=FALSE)

    Sigma_out_fname <- paste0(out.head, sprintf("_%s_L%d.%s", if (is.IW) "Sigma" else if (is.IG) "sigma2", 
						L, if (is.IW) out.format else "csv"))

    if (is.IW) {
        # TODO: involve seed_list here
	if (is.null(seed_list)) {
            Sigmas <- rinvwishart(L, nT, DT)
	} else {
            Sigmas <- array(dim=c(c, c, L))
            if (is.numeric(seed_list)) {
                if (seedarray_explicit) {
                    for (l in 1:L) {
                        set.seed(seed_list[seed_ix])
                        Sigmas[,,l] <- rinvwishart(1, nT, DT)[,,1]
		        seed_ix <- seed_ix + 1
		    }
		} else {
                    for (l in 1:L) {
                        set.seed(seed_now)
                        Sigmas[,,l] <- rinvwishart(1, nT, DT)[,,1]
                        seed_now <- seed_now + seed_step
    	            }
		}
	    } else {
                #TODO: read in the first seed_perchunk from seed_list
                n_chunks <- ceil(L/seed_perchunk)
	        seed_ix <- 0
	        for (i in 1:n_chunks) {
                    #TODO: Read in that number of chunks and generate Sigma
                    #TODO: 

                    for (j in 1:seed_perchunk) {
                        #TODO: randomize 
                        

                        if (n_chunks * seed_perchunk + j == L) break
		    }
                    seed_ix <- seed_ix + min(seed_perchunk, j)
		}

	    }
	}
        Sigma_mat <- matrix(nrow=L, ncol=(c * (c + 1))/2, data=0)

	for (l in 1:L) Sigma_mat[l,] <- c(Sigmas[,,l][upper.tri(Sigmas[,,l], diag=TRUE)])
	if (out.format == "csv") {
            write.csv(as.data.frame(Sigma_mat), file=Sigma_out_fname, row.names=FALSE)
	} else if (out.format == "fst") {
            write.fst(as.data.frame(Sigma_mat), Sigma_out_fname)
	}
    }

    if (is.IG) {
        # TODO: involve seed_list here
        if (is.null(seed_list)) {
            sigma2s <- 1/rgamma(L, shape=nT, rate=dT)
	} else {
            #TODO: continue here
	}
        data.table::fwrite(as.list(sigma2s), file=Sigma_out_fname)
    }

    #TODO: Stplus1 check if it is pds. If not, throw an error.

    Thetat_out_fname <- paste0(out.head, "_Thetat_t%dk%d.%s")

    for (t in T:1) {

        # Read in the relevant m_file_list and M_file_list element
        mt_file <- if (out.format == "csv") {
		       read.csv(m_file_list[t])
	           } else if (out.format == "fst") {
		       read.fst(m_file_list[t])
	           }
        Mt_file <- if (out.format == "csv") {
		       read.csv(M_file_list[t])
	           } else if (out.format == "fst") {
		       read.fst(M_file_list[t])
	           }

	if (change_Gt) Gtplus1 <- matrix(Gtmat[t,], nrow=p, ncol=p)
        if (change_Wt) Wtplus1 <- matrix(Wtmat[t,], nrow=p, ncol=p)

	# K+1 is the length of the data frame.
	K <- nrow(mt_file) - 1

        st_list <- matrix(nrow=K+1, ncol=p*c, data=0)
        St_list <- matrix(nrow=K+1, ncol=as.integer(p * (p+1)/2), data=0)

	st_list[K+1,] <- c(stplus1)
	St_list[K+1,] <- Stplus1[upper.tri(Stplus1, diag=TRUE)]

        st_out_fname <- paste0(out.head, sprintf("_st_%d.%s", t, out.format))
        St_out_fname <- paste0(out.head, sprintf("_St_uptri_%d.%s", t, out.format))

        st_out_fname_list <- c(st_out_fname_list, st_out_fname)
        St_out_fname_list <- c(St_out_fname_list, St_out_fname)

        for (k in K:1) {
	    mt <- matrix(unlist(mt_file[k,]), nrow=p, ncol=c)
	    # This is the upper triangle of Mt. Copy the entries to make it lower-triangular.
            Mt_upper <- unlist(Mt_file[k,])
	    Mt <- tri.to.sym(Mt_upper, p, lower=FALSE)
		    
            uMt <- chol(Mt)

            # compute uAt from GMG' + W.
	    uAtplus1 <- uMt %*% t(Gtplus1)
	    uAtplus1 <- chol(crossprod(uAtplus1) + Wtplus1)
	   
            BS.param.list <- if (is.IW) {
		                 BS_1step(stplus1, Stplus1, mt, Gtplus1, Mt, uAtplus1, L, Sigmas = Sigmas)
	                     } else if (is.IG) {
		                 BS_1step(stplus1, Stplus1, mt, Gtplus1, Mt, uAtplus1, L, sigma2s = sigma2s, R = R)
	                     }

	    stplus1 <- BS.param.list$s
	    Stplus1 <- BS.param.list$S

	    # Generate L Theta_t samples.
	    Theta_samples <- matrix(nrow=L, ncol=p*c, data=0)
	    for (l in 1:L) {
                Theta_samples[l,] <- c(BS.param.list$Theta_t[,,l])
	        if (verbose & (l %% 1000 == 0)) print(sprintf("%d/%d sample generated for Theta_tk for t = %d, k = %d.", l, L, t, k))
	    }

            if (out.format == "csv") {
                write.csv(as.data.frame(Theta_samples), file=sprintf(Thetat_out_fname, t, k, out.format), row.names=FALSE)
	    } else if (out.format == "fst") {
                write.fst(as.data.frame(Theta_samples), sprintf(Thetat_out_fname, t, k, out.format), compress=100)
	    }

	    st_list[k,] <- c(stplus1)
	    St_list[k,] <- Stplus1[upper.tri(Stplus1, diag=TRUE)]
	}

	#output files
	if (out.format == "csv") {
	    write.csv(st_list, file=st_out_fname, row.names=FALSE)
            write.csv(St_list, file=St_out_fname, row.names=FALSE)
	} else if (out.format == "fst") {
	    write.fst(as.data.frame(st_list), st_out_fname)
            write.fst(as.data.frame(St_list), St_out_fname)
	}
    }

    # return the list of smoothing files.
    out <- list(st_files = st_out_fname_list, St_files = St_out_fname_list)
    class(out) <- "BS.bigdata"
    return(out)
}


