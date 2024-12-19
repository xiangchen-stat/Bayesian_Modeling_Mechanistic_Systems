library(rlang)
#library(Formula)
library(collections)
library(data.table)
library(stringr)
library(fst)

#' one time step for the Forward Filter. Computes the FF parameters given the data at the relevant time step and the relevant parameters from the last time step.
#' @param Yt The data matrix for epoch t.
#' @param Ft The matrix of covariates F_t.
#' @param Gt G_t, the beta transition matrix.
#' @param mtmin1 mean of beta_{t-1} | Y_{1:(t-1)}
#' @param Mtmin1 left-covariance matrix of beta_{t-1} | Y_{1:(t-1)}
#' @param Wt left-covariance matrix of the noise parameter of beta_{t-1} | Y_{1:(t-1)}
#' @param Vt left-covariance matrix of the noise parameter of Y
#' @param ntmin1 the shape parameter, or degrees of freedom, of the right-covariance matrix Sigma | Y_{1:(t-1)}
#' @param Dtmin1 the scale matrix of the right-covariance matrix Sigma | Y_{1:(t-1)}
#' @param N the number of rows of Y
#' @param omega The right-variance matrix discount factor 
#' @returns The mean and covariance matrices m_t and M_t of one filtering step, the updated inverse-Wishart parameters n_t and D_t for the right-covariance matrix, plus other relevant parameters.

FF_1step <- function(Yt, Ft, Gt, mtmin1, Mtmin1, Wt, Vt, ntmin1, Dtmin1, N, omega){#, scale.random=TRUE) {
    if (is.null(omega)) omega <- 1.0

    #TODO: substitute with the .C functions.
    At <- tcrossprod(chol(Mtmin1), Gt)
    At <- crossprod(At) + Wt

    at <- Gt %*% mtmin1
    Qt <- Ft %*% tcrossprod(At, Ft) + Vt
    rownames(Qt) <- NULL
    colnames(Qt) <- NULL
    qt <- Ft %*% at

    uQt <- chol(Qt)
    lQinv_e <- backsolve(uQt, Yt - qt, transpose=TRUE)
    nt <- ntmin1 * omega + N
    Dt <- Dtmin1 * omega + crossprod(lQinv_e)

    lQinvFtA <- backsolve(uQt, Ft %*% At, transpose = TRUE)

    mt <- at + crossprod(lQinvFtA, lQinv_e)
    #	  R %*% t(Ft) %*% Qinv_e
    Mt <- At - crossprod(lQinvFtA)
    Mt <- check.pds(Mt)

    # TODO: limit the parameters that return because not all of them are used for following computations.
    #   Check which ones.
    return(list(at = at, uAt = chol(At), uQt = uQt, qt = qt, nt = nt, Dt = Dt,
    	    mt = mt, Mt = Mt, omega = omega))
}

#TODO: temporary function. Consider replacing or deprecating the above so that we have the pure C implementation
FF_1step_c <- function(Yt, Ft, Gt, mtmin1, Mtmin1, Wt, Vt, ntmin1, Dtmin1, N, omega){
    if (is.null(omega)) omega <- 1.0
    p <- nrow(Gt)
    S <- ncol(Yt)
#    symtol <- 1e-9
    #TODO: Do the above, but replace with .C functions instead of writing one fresh. Try here.

    return(.Call("FF_1step", Yt, Ft, Gt, mtmin1, Mtmin1, Wt, Vt, ntmin1, Dtmin1,
		 N, p, S, omega))
}


#' one time step for the Forward Filter. Computes the FF parameters given the data at the relevant time step and the relevant parameters from the last time step. This one does so where Sigma has a specific structure and only its scale quantity is updated, i.e. Sigma = sigma2 * R_t, for some fixed matrix R_t. sigma2 is assumed to be distributed as an inverse-Gamma.
#' @param Yt The data matrix for epoch t.
#' @param Ft The matrix of covariates F_t.
#' @param Gt G_t, the beta transition matrix.
#' @param mtmin1 mean of beta_{t-1} | Y_{1:(t-1)}
#' @param Mtmin1 left-covariance matrix of beta_{t-1} | Y_{1:(t-1)}
#' @param Wt left-covariance matrix of the noise parameter of beta_{t-1} | Y_{1:(t-1)}
#' @param Vt left-covariance matrix of the noise parameter of Y
#' @param ntmin1 the shape parameter, or degrees of freedom, of the scale for the right-covariance matrix sigma2 | Y_{1:(t-1)}
#' @param dtmin1 the scale matrix of the scale for the right-covariance matrix sigma2 | Y_{1:(t-1)}
#' @param R the structural matrix for Sigma, where Sigma = sigma2 * R.
#' @param N the number of rows of Y
#' @param S the number of columns of Y
#' @param omega The right-variance matrix discount factor 
#' @returns The mean and covariance matrices m_t and M_t of one filtering step, the updated inverse-Gamma parameters n_t and D_t for the right-covariance matrix, plus other relevant parameters.

FF_1step.MNIG <- function(Yt, Ft, Gt, mtmin1, Mtmin1, Wt, Vt, ntmin1, dtmin1, R, N, S, omega){#, scale.random=TRUE) {
    if (is.null(omega)) omega <- 1.0

    #TODO: substitute with the .C functions
    At <- tcrossprod(chol(Mtmin1), Gt)
    At <- crossprod(At) + Wt

    at <- Gt %*% mtmin1
    Qt <- Ft %*% At %*% t(Ft) + Vt
    rownames(Qt) <- NULL
    colnames(Qt) <- NULL
    qt <- Ft %*% at

    uQt <- chol(Qt)
    lQinv_e <- backsolve(uQt, Yt - qt, transpose=TRUE)

    uR <- chol(R)
    lQinv_euRinv <- backsolve(uR, t(lQinv_e), transpose=TRUE)

    nt <- ntmin1 * omega + N * S/2
    dt <- dtmin1 * omega + sum(lQinv_euRinv*lQinv_euRinv)/2

    lQinvFtA <- backsolve(uQt, Ft %*% At, transpose = TRUE)

    mt <- at + crossprod(lQinvFtA, lQinv_e)
    #	  R %*% t(Ft) %*% Qinv_e
    Mt <- At - crossprod(lQinvFtA)
    Mt <- check.pds(Mt)

    return(list(at = at, uAt = chol(At), uQt = uQt, qt = qt, nt = nt, dt = dt,
    	    mt = mt, Mt = Mt, omega = omega))
}



#' The Forward Filtering step. Generates a list of parameters needed to compute the backwards-sampling or smoothing steps.
#' @param formula_y The R formula object denoting the set of variables to treat as the outcome. Variables to be written on right-hand-side only. Excludes an intercept term. (e.g. `~ y1 + y2 + 0`)
#' @param formula_x The R formula object denoting the set of variables to treat as covariates. Variables to be written on right-hand-side only. (e.g. `~ x1 + x2 + x3`)
#' @param data The data frame to be regressed over by the DLM. Must include the variables in formula_y, formula_x, tvar, and indexvar where it is specified.
#' @param n0 The prior degrees of freedom for the right-covariance matrix Sigma.
#' @param D0 The prior scale matrix for the right-covariance matrix Sigma.
#' @param tvar The discrete time variate over which the FF will iterate. Default: "epoch"
#' @param m0 The prior mean matrix for the FF. If not specified, will be initialized to a zero-matrix.
#' @param M0 The pior covariance matrix for the FF. Defaults to diag(p), where p is the number of rows of m0.
#' @param Gt G_t, the beta transition matrix. Assumed known and constant throughout for the purposes of this implementation.
#' @param Wt W_t, the left-covariance matrix of beta_t. Assumed known and constant throughout for the purposes of this implementation.
#' @param Vt The function used to generate the correlation matrix between the Y values. It takes in the size of a chunk of Y and outputs Vt. Defaults to function(chunk) diag(nrow(chunk)).
#' @param indexvar The name of the outcome variable or variables. If set, will sort the data by the variable or set of variables prior to performing any computations.
#' @param discount.factor The right-variance matrix discount factor. If not specified, will default to 1.0.
#' @returns A list of computed values relevant to the smoothing step, as well as recording some input value. Includes m, the list of mt's computed; M for the list of Mt's, nt, Dt, uAt, and uQt follow. N.list is the list of the number of rows in Y at each time step. 
#' @export

#TODO: Substitute Gt and Wt for arrays of matrices.
#TODO: m0 should be a specified vector; M0 should be a specified matrix.
#TODO: Since n_t is expected to be constant, input Vt directly a'la Gt and Wt.

FF <- function(formula_y, formula_x, data, n0, D0, tvar="epoch",
		    m0 = NULL, M0 = "I_p", Gt = "I_p", Wt = "I_p",
		    Vt = function(chunk) diag(nrow(chunk)), indexvar = NULL, discount.factor=NULL) {

    if (!is.null(discount.factor) & !is.numeric(discount.factor)) {
        stop("Discount factor must be numeric if specified.")
    }
    if (!is.null(indexvar)) data <- data %>% arrange(all_of(indexvar))

    form_x_t <- update(formula_x, paste("~. +",tvar))
#    formula_x <- Formula(formula_x)
#    form_x_t <- Formula(form_x_t)

    mfX_t <- stats::model.frame(form_x_t, data)
#    mf_t <- stats::model.frame(form_x_t, data)
#    mt_t <- attr(x = mf_t, which = "terms")
    mtX_t <- attr(x = mfX_t, which = "terms")

    # this gives us X
    #X <- as.data.frame(model.matrix(object = mt, data = mf))
    # we want x with the epoch for processing
    X_aug <- as.data.frame(model.matrix(object = mtX_t, data = mfX_t)) %>%
	    rename(epoch = all_of(tvar))

    x.names <- colnames(X_aug)
    x.names <- x.names[x.names != "epoch"]
    
    # get y
#    Y <- model.response(mf_t, "numeric")
    
    # updates formula_y so that an intercept does not appear here.
    if (!("0" %in% all.vars(formula_y))) formula_y <- update(formula_y, "~. + 0")
    form_y_t <- update(formula_y, paste("~. +",tvar))
    mfY_t <- stats::model.frame(form_y_t, data)
    mtY_t <- attr(x = mfY_t, which = "terms")
    Y_aug <- as.data.frame(model.matrix(object = mtY_t, data = mfY_t)) %>%
	    rename(epoch = all_of(tvar))

    y.names <- colnames(Y_aug)
    y.names <- y.names[y.names != "epoch"]
#    Y_aug <- data %>% select(all_of(c(tvar, y.names))) %>%
#	    rename(epoch = all_of(tvar))
	    #data.frame(Y = model.response(mf_t, "numeric"),
            #            epoch = X_aug[[tvar]])


    #TODO: add support for random effects. We can group them into Vt

    p <- ncol(X_aug) - 1
    m <- ncol(D0)
    if (is.null(m)) m <- 1 # if B0 is a scalar
  
    # sort out the DLM-specific parameters
    # G_t and W_t here. If unspecified default to I(p)
    if (Gt == "I_p") Gt <- diag(p)
    if (Wt == "I_p") Wt <- diag(p)
    if (is.null(m0)) m0 <- matrix(rep(0,p * m), nrow=p)
    if (M0 == "I_p") M0 <- diag(p)
    
    # get all t's.
    t_list <- sort(unique(data[[tvar]])) #%>% 
#	               unique() %>%
#		       arrange(get(tvar)) %>%
#		       pull(all_of(tvar))

    # set an unambiguous prior t
    T <- length(t_list)
    N_t.list <- rep(0, T)
    
    m_list <- matrix(rep(0, (T+1) * m * p), nrow = T+1, ncol = m * p)
    M_list <- matrix(rep(0, (T+1) * p * p), nrow = T+1, ncol = p*p)
    
    n_list <- rep(0, T+1)
    D_list <- matrix(rep(0, (T+1) * m * m), nrow = T+1, ncol = m*m)
    n_list[1] <- n0
    D_list[1,] <- c(D0)
    
    uAt_list <- matrix(rep(0, T * p * p), nrow = T, ncol = p*p)
    uQt_list <- dict()
    
    m_list[1,] <- c(m0)
    M_list[1,] <- c(M0)
    
    #Y.ix.0 <- 0
    for (i in 1:T) {
        t <- t_list[i]
        
        Xt <- as.matrix(X_aug %>% filter(epoch == t) %>% select(-c(epoch)))
        N <- nrow(Xt)
        N_t.list[i] <- N

        Yt <- Y_aug %>% filter(epoch == t) %>% select(-c(epoch))
        
        # select Vt here. If Vt is null, didir.create(outdir, showWarnings=FALSE)agonal by default
        #TODO: Write another function to select a default Vt, e.g. Vt, but with some negative correlation between consecutive days of exercise, based on the string or object input. Use the autocorrelation function for such interfaces. We may need a separate list object to pass in the autocorrelation parameter.
        V_t <- Vt(Xt)
        
        Mt <- matrix(M_list[i,], p, p)

        FF.param.list <- FF_1step(Yt, Xt,
                                  Gt, matrix(m_list[i,], nrow=p), Mt,
                                  Wt, V_t,
                                  n_list[i], matrix(D_list[i,], nrow=m), N, discount.factor)

        m_list[i+1,] <- c(FF.param.list$mt)
        M_list[i+1,] <- c(FF.param.list$Mt)
        n_list[i+1] <- FF.param.list$nt
        D_list[i+1,] <- c(FF.param.list$Dt)
        
        uAt_list[i,] <- c(FF.param.list$uAt)
        uQt_list$set(t_list[i], FF.param.list$uQt)
        
        #Y.ix.0 <- Y.ix.0 + N
    }
    
    # Give m_list the names of Xt as the heads for maximal clarity	
    xy.names <- expand.grid(y.names, x.names)
    xy.names$xy.names <- paste(xy.names$Var1, xy.names$Var2, sep=".")
    colnames(m_list) <- xy.names$xy.names

    #TODO: return the function's call the same way lm() does this
    return(list(formula_y = formula_y,
		formula_x = formula_x,
		m = m_list, M = M_list, n = n_list, D = D_list,
                uAt = uAt_list, uQt = uQt_list, N.list = N_t.list, 
		omega = discount.factor))
}

#TODO: May have to use useMethod() for this. But how to do it with multiple arguments?

#' The Forward Filter, but adapted to an enormous data setting where the size of the output is an enormous X x S matrix, X and S both large, and the file needs to be read into memory in sequence. 
#' @param data_set The big_data_set object containing the list of items to be iterated over by the FF.
#' @param n0 The prior degrees of freedom for the right-covariance matrix Sigma.
#' @param D0 The prior scale matrix for the right-covariance matrix Sigma.
#' @param m0 The prior mean matrix for the FF. If not specified, will be initialized to a zero-matrix.
#' @param M0 The prior covariance matrix for the FF. Defaults to the identity matrix.
#' @param Gtmat G_t, the matrix of beta transition matrices concatenated over time, with each row being a flattened G_t. A single G_t may be set for all time by passing in one covariance matrix flattened into a one-row matrix. Assumed known and constant throughout for the purposes of this implementation. Defaults to the identity matrix if not specified.
#' @param Wtmat W_t, the left-covariance matrices of  of beta_t, with each row being a flattened W_t. A single W_t may be set for all time by passing in one covariance matrix flattened into a one-row matrix. Assumed known and constant throughout for the purposes of this implementation. Defaults to the identity matrix if not specified.
#' @param Vt The function used to generate the correlation matrix between the Y values, given the size of the Y. Defaults to function(chunk) diag(nrow(chunk)).
#' @param omega The right-variance matrix discount factor. If not specified, will default to 1.0.
#' @param out.head The header of the files output by this function.
#' @param out.format The format of the files to read and write. Supports 'csv' and 'fst'. 'fst' is a fast file I/O format developed by Hadley Wickham. Defaults to 'csv'.
#' @param fst.head The header of the name for the columns of both Y_t and F_t. Required to parse the relevant columns from the relevant files. Defaults to 'V' (i.e. the columns of both Y_t and F_t are labeled V1, V2, V3,...).
#' @returns A list of computed values relevant to the smoothing step, as well as recording some input value. Includes mt_list, the list of m's computed; Mt_list for the list of M's, nt_list, and Dt_list.
#' @export

FF.bigdata <- function(data_set, m0, M0, n0, D0 = NULL, 
		       d0 = NULL, R = NULL, Gtmat = NULL, Wtmat = NULL,
		       Vt = function(chunk) diag(nrow(chunk)),
		       omega = NULL, 
#		       tvar="epoch", 
		       out.head = "./ff_out", #grid.traversal.mode = "rowsnake",
		       out.format = "csv",
		       fst.head = "V") {

    if (!(out.format %in% c("csv", "fst"))) stop("Output format of FF must be csv or fst!")

    mtmin1 <- m0
    Mtmin1 <- M0
    ntmin1 <- n0
    Dtmin1 <- D0
    dtmin1 <- d0

    iw_mode <- !is.null(D0)
    ig_mode <- !(is.null(d0) | is.null(R))

    if (!iw_mode & !ig_mode) stop("Either D0, or d0 and R, must be specified.")

    if (iw_mode & ig_mode) {
        warning("D0 and d0 and R are all specified! Defaulting to the inverse Wishart D0; d0 and R will not be considered.")
        ig_mode <- FALSE
    }

    # read in the columns of F
    F_file <- data_set$get("F1")
    file_name_F <- F_file$file_name
    p <- F_file$ncol
#    p <- as.integer(system(sprintf("awk -F, '{print NF; exit}' %s", file_name_F), intern=TRUE))
#    p <- ncol(read_csv(file_name_F, col_names = FALSE, n_max = 1,
#               show_col_types=FALSE))
#	    as.integer(system(sprintf("head -n1 %s | awk -F, '{print NF}'", 
#					  file_name_F), intern=TRUE))
    
    #NOTE: p can also be big
    split_col_F <- data_set$get("split_col_F")
    ptilde <- p
    if (split_col_F) ptilde <- F_file$ncolblock

    # set default Gtmat and Wtmat if they aren't specified.
    if (is.null(Gtmat)) Gtmat <- matrix(c(diag(ptilde)), nrow=1)
    if (is.null(Wtmat)) Wtmat <- matrix(c(diag(ptilde)), nrow=1)

    # whether or not to transition Gt across time.
    change_Gt <- nrow(Gtmat) > 1
    Gt <- NULL
    if (!change_Gt) Gt <- matrix(Gtmat, nrow=ptilde, ncol=ptilde)

    # whether or not to transition Wt across time.
    change_Wt <- nrow(Wtmat) > 1
    Wt <- NULL
    if (!change_Wt) Wt <- matrix(Wtmat, nrow=ptilde, ncol=ptilde)

    # only Y's here
    data_file_list <- unlist(data_set$keys())
    data_file_list <- str_subset(data_file_list, "Y")
    T <- length(data_file_list)

    mt_out_fname_list <- c()
    Mt_out_fname_list <- c()
    nt_out_fname_list <- c()
    Dt_out_fname_list <- c()
    dt_out_fname_list <- c()

    # insert progress bar here.
    #TODO: Use a similar progress bar style as read_csv(). 
#    pbT <- txtProgressBar(min = 0, max = T,
#			  style = 3, char = "=")

    for (data_file_key in data_file_list) {
        t <- as.integer(substring(data_file_key, 2))

	data_file <- data_set$get(data_file_key)
        file_name <- data_file$file_name

	F_file <- data_set$get(paste0("F", t))
	file_name_F <- F_file$file_name
        # get the number of rows and the number of columns in the file without opening it.
#	fnrow <- as.integer(system(sprintf('cat %s | wc -l', 
#					  file_name), intern=TRUE)) - 
#		as.integer(data_file$header)
	fnrow <- data_file$nrow
#	fnrow <- as.integer(str_split_1(system(sprintf("wc -l %s", file_name), intern=TRUE), ' ')[1]) -
#		as.integer(data_file$header)

#        fncol <- as.integer(system(sprintf("head -n1 %s | awk -F, '{print NF}'", 
#					  file_name), intern=TRUE))
	fncol <- data_file$ncol
#	fncol <- as.integer(system(sprintf("awk -F, '{print NF; exit}' %s", file_name), intern=TRUE))

        r <- min(fnrow, data_file$nrowblock)
        c <- min(fncol, data_file$ncolblock)

	# generate grid
        grid <- data_file$grid #generate_grid(fnrow, fncol, r, c, traversal.mode = grid.traversal.mode)

        c_F <- min(p, F_file$ncolblock)
	# only if split_col_F
	F_grid <- NULL
	if (split_col_F) {
            F_grid <- F_file$grid #generate_grid(fnrow, p, r, c_F, traversal.mode = grid.traversal.mode)
	}

	K <- as.integer(fnrow/r) * as.integer(fncol/c)

	# allocate matrices and initialize
        mt_list <- matrix(nrow=K+1, ncol = ptilde * c, data = 0)
	Mt_list <- matrix(nrow=K+1, ncol = as.integer(ptilde * (ptilde+1)/2), data=0)
        nt_list <- rep(0, K+1)
	Dt_list <- matrix(nrow=K+1, ncol = as.integer(c * (c + 1)/2), data=0)
	dt_list <- rep(0, K+1)

	mt_list[1,] <- c(mtmin1)
	Mt_list[1,] <- Mtmin1[upper.tri(Mtmin1, diag=TRUE)]
	nt_list[1] <- ntmin1

	if (iw_mode) Dt_list[1,] <- Dtmin1[upper.tri(Dtmin1, diag=TRUE)]
	if (ig_mode) dt_list[1] <- dtmin1

	#TODO: edit here to allow for fst
	# set output file names
        mt_out_fname <- paste0(out.head, sprintf("_mt_%d.%s", t, out.format))
        Mt_out_fname <- paste0(out.head, sprintf("_Mt_uptri_%d.%s", t, out.format))
        nt_out_fname <- paste0(out.head, sprintf("_nt_%d.csv", t))
        
        if (iw_mode) Dt_out_fname <- paste(out.head, sprintf("Dt_uptri_%d.%s", t, out.format), sep="_")
	if (ig_mode) dt_out_fname <- paste(out.head, sprintf("dt_%d.csv", t), sep="_")

	# store in a file
        mt_out_fname_list <- c(mt_out_fname_list, mt_out_fname)
        Mt_out_fname_list <- c(Mt_out_fname_list, Mt_out_fname)
        nt_out_fname_list <- c(nt_out_fname_list, nt_out_fname)

        if (iw_mode) Dt_out_fname_list <- c(Dt_out_fname_list, Dt_out_fname)
        if (ig_mode) dt_out_fname_list <- c(dt_out_fname_list, dt_out_fname)

	#TODO: do we want this to change for each block?
	if (change_Gt) Gt <- matrix(Gtmat[t,], nrow=ptilde, ncol=ptilde)
	if (change_Wt) Wt <- matrix(Wtmat[t,], nrow=ptilde, ncol=ptilde)

#        pbK <- txtProgressBar(min = 0, max = K,
#				style = 3, char = "-")

	for (k in 1:K) {

	    ud <- unlist(grid[grid$block == k, c("U", "D")])
	    lr <- unlist(grid[grid$block == k, c("L", "R")])

            # read in the Y
            Yt <- if (out.format == "csv") {
		      as.matrix(read_big_csv(file_name, rows = ud,
	  		        cols = unlist(grid[grid$block == k, c("L", "R")]), 
			        has_header=data_file$header,
			        ignore_header=TRUE))
	          } else if (out.format == "fst") {
                      as.matrix(read_big_fst(file_name, rows = ud,
					     cols = paste0(fst.head, seq(lr[1], lr[2]))))
	          }
	    N <- ud[2] - ud[1] + 1

	    # read in the F
            Ft <- NULL
	    if (is.null(F_grid)) {
	        Ft <- if (out.format == "csv") {
		          as.matrix(read_big_csv(file_name_F, rows = ud,
					     has_header=data_file$header,
					     ignore_header=TRUE))
		      } else if (out.format == "fst") {
		          as.matrix(read_big_fst(file_name_F, rows = ud))
		      }
	    } else {
	        ud_F <- unlist(F_grid[F_grid$block == k, c("U", "D")])
	        lr_F <- unlist(F_grid[F_grid$block == k, c("L", "R")])
                Ft <- if (out.format == "csv") {
		          as.matrix(read_big_csv(file_name_F, rows = ud_F,
			            cols = unlist(F_grid[F_grid$block == k, c("L", "R")]), 
			            has_header=F_file$header,
			            ignore_header=TRUE))
		      } else if (out.format == "fst") {
		          as.matrix(read_big_fst(file_name_F, rows = ud_F,
			            cols = paste0(fst.head, seq(lr_F[1], lr_F[2]))))
		      }
	    }

#	    #TODO: useMethod() over here to call the FF_1step default or the MNIG based on the argument structure.
#            #TODO: debug. Print out everything to see that the computation is correct. (It should be.) Manually do the first step yourself to verify the answer (since M1 is apparently not PD).
#	    if ((t == 1) & (k == 1)) {
#                write.csv(as.data.frame(Yt), file="/home/hunanzhou/my_c/Banerjee_GSR/projects/MNIW/Yt.csv", row.names=FALSE)
#                write.csv(as.data.frame(Ft), file="/home/hunanzhou/my_c/Banerjee_GSR/projects/MNIW/Ft.csv", row.names=FALSE)
#                write.csv(as.data.frame(Gt), file="/home/hunanzhou/my_c/Banerjee_GSR/projects/MNIW/Gt.csv", row.names=FALSE)
#                write.csv(as.data.frame(Wt), file="/home/hunanzhou/my_c/Banerjee_GSR/projects/MNIW/Wt.csv", row.names=FALSE)
#                write.csv(as.data.frame(Vt(Ft)), file="/home/hunanzhou/my_c/Banerjee_GSR/projects/MNIW/Vt.csv", row.names=FALSE)
#                write.csv(as.data.frame(mtmin1), file="/home/hunanzhou/my_c/Banerjee_GSR/projects/MNIW/m0.csv", row.names=FALSE)
#                write.csv(as.data.frame(Mtmin1), file="/home/hunanzhou/my_c/Banerjee_GSR/projects/MNIW/M0.csv", row.names=FALSE)
#		write.csv(data.frame(n0 = ntmin1, N = N, is_null_omega = is.null(omega)), file="/home/hunanzhou/my_c/Banerjee_GSR/projects/MNIW/n0_N_omega.csv", row.names=FALSE)
#
#	    }

            # FF_1step over the grid elements, using the output of the next elements as priors for the previous elements.
            if (iw_mode) FF.param.list <- FF_1step(Yt, Ft, Gt, mtmin1, Mtmin1, Wt, Vt(Ft), ntmin1, Dtmin1, N, omega)
	    if (ig_mode) FF.param.list <- FF_1step.MNIG(Yt, Ft, Gt, mtmin1, Mtmin1, Wt, Vt(Ft), ntmin1, dtmin1, R, N, ncol(Yt), omega)

	    # set for next iteration.
	    mtmin1 <- FF.param.list$mt
	    Mtmin1 <- FF.param.list$Mt
	    ntmin1 <- FF.param.list$nt
	    if (iw_mode) Dtmin1 <- FF.param.list$Dt
	    if (ig_mode) dtmin1 <- FF.param.list$dt

	    # Save half the matrix to avoid issues.
	    mt_list[k+1,] <- c(mtmin1)
	    Mt_list[k+1,] <- Mtmin1[upper.tri(Mtmin1, diag=TRUE)]
	    nt_list[k+1] <- ntmin1
	    if (iw_mode) Dt_list[k+1,] <- Dtmin1[upper.tri(Dtmin1, diag=TRUE)]
	    if (ig_mode) dt_list[k+1] <- dtmin1

        }

	# output recorded matrices to file.
	if (out.format == "csv") {
            write.csv(mt_list, file=mt_out_fname, row.names=FALSE)
            write.csv(Mt_list, file=Mt_out_fname, row.names=FALSE)
	} else if (out.format == "fst") {
            write.fst(as.data.frame(mt_list), mt_out_fname)
            write.fst(as.data.frame(Mt_list), Mt_out_fname)
	}
        data.table::fwrite(as.list(nt_list), file=nt_out_fname)
	if (iw_mode) {
	    if (out.format == "csv") write.csv(Dt_list, file=Dt_out_fname, row.names=FALSE)
	    if (out.format == "fst") write.fst(as.data.frame(Dt_list), Dt_out_fname)
	}
	if (ig_mode) data.table::fwrite(as.list(dt_list), file=dt_out_fname)

#	setTxtProgressBar(pbT, t)
    }

#    close(pbT)

    out <- list(nt_files = nt_out_fname_list,
                mt_files = mt_out_fname_list,
                Mt_files = Mt_out_fname_list)

    if (iw_mode) out$Dt_files <- Dt_out_fname_list
    if (ig_mode) out$dt_files <- dt_out_fname_list

    # Return the path to the file streams
    return(out)
}

