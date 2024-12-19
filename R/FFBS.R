library(rlang)

#' The Forward Filter-Backwards Sample algorithm using the big data set structure.
#' @param data_set The big_data_set object containing the list of items to be iterated over by the FF.
#' @param m0 The prior mean matrix for the FF. If not specified, will be initialized to a zero-matrix.
#' @param M0 The prior covariance matrix for the FF. Defaults to the identity matrix.
#' @param n0 The prior degrees of freedom for the right-covariance matrix Sigma.
#' @param D0 The prior scale matrix for the right-covariance matrix Sigma.
#' @param d0 The prior scale for the scale sigma2 where Sigma = sigma2 * R, sigma2 being inverse-gamma.
#' @param Gtmat G_t, the matrix of beta transition matrices concatenated over time, with each row being a flattened G_t. A single G_t may be set for all time by passing in one covariance matrix flattened into a one-row matrix. Assumed known and constant throughout for the purposes of this implementation. Defaults to the identity matrix if not specified.
#' @param Wtmat W_t, the left-covariance matrices of  of beta_t, with each row being a flattened W_t. A single W_t may be set for all time by passing in one covariance matrix flattened into a one-row matrix. Assumed known and constant throughout for the purposes of this implementation. Defaults to the identity matrix if not specified.
#' @param Vt The function used to generate the correlation matrix between the Y values, given the size of the Y. Defaults to function(chunk) diag(nrow(chunk)).
#' @param n_samples The number of backwards samples to draw. Set to NULL to return only the smoothing parameters instead. Defaults to 100.
#' @param seed_list A list of random seeds to run. Each seed will control for the next quantity to be randomized, first the L samples of Sigma or sigma2, then L samples of each Theta_tk. Can be passed as a file, as a full set of L * (T * K + 1) integers, or as a sequence of 2 or 3 integers. If as a file, the file will contain at least L * (T * K + 1) integers to be used as seeds. If as a sequence of 3 integers (a,b,c), the seeds will be taken from the sequence (a,b) with a step size of c. If as a sequence of 2 integers, the step size c will be treated as 1.
#' @param ff.out.head The header of the files output by the forward filter.
#' @param bs.out.head The header of the files output by the backwards sample.
#' @param out.format The format of the files to read and write. Supports 'csv' and 'fst'. 'fst' is a fast file I/O format developed by Hadley Wickham. Defaults to 'csv'.
#' @param fst.head The header of the name for the columns of both Y_t and F_t. Required to parse the relevant columns from the relevant files. Defaults to 'V' (i.e. the columns of both Y_t and F_t are labeled V1, V2, V3,...).
#' @param verbose Whether or not to print every 1000 iterations for backwards sampling. Default: FALSE
#' @returns The list of forward filter and backwards sampling files for each time point.
#' @export

FFBS.bigdata <- function(data_set, m0, M0, n0, D0 = NULL, 
                         d0 = NULL, R = NULL, Gtmat = NULL, Wtmat = NULL,
                         Vt = function(chunk) diag(nrow(chunk)),
                         n_samples = 100, seed_list = NULL,
                         ff.out.head = "./ff_out", #grid.traversal.mode = "rowsnake",
                         bs.out.head = "./bs_out", #grid.traversal.mode = "rowsnake",
                         out.format = "csv",
                         fst.head = "V",
                         verbose=FALSE) {

    is.IW <- !is.null(D0)
    is.IG <- !(is.null(d0) | is.null(R))
    if (!is.IW & !is.IG) stop("FFBS column covariance matrix must be specified to be IW or IG (x scale matrix).")

    ff.file.paths <- FF.bigdata(data_set, m0, M0, n0, D0, 
                                d0, R, Gtmat, Wtmat,
                                Vt, omega = NULL, 
#                                tvar="epoch", 
                                out.head = ff.out.head, #grid.traversal.mode = "rowsnake",
                                out.format = out.format,
                                fst.head = fst.head)

    T <- length(ff.file.paths$nt_files)

    p = nrow(m0)
    c = ncol(m0)

    nT <- unlist(read.csv(ff.file.paths$nt_files[T], header=FALSE))
    names(nT) <- NULL
    nT <- nT[length(nT)]

    DT <- NULL
    dT <- NULL

    if (is.IW) {
        DT <-  if (out.format == "csv") {
                   read.csv(ff.file.paths$Dt_files[T], header=TRUE)
               } else if (out.format == "fst") {
                   read.fst(ff.file.paths$Dt_files[T])
               }
    
        DT <- as.numeric(unlist(DT[nrow(DT),]))
        names(DT) <- NULL
        DT <- tri.to.sym(DT, c, lower=FALSE)
    }
    if (is.IG) {
        dT <- unlist(read.csv(ff.file.paths$dt_files[T], header=FALSE))
        names(dT) <- NULL
        dT <- dT[length(dT)]
    }

    # Backwards sampling here.
    BS.file.paths <- if (is.null(n_samples)) smooth.bigdata(ff.file.paths$mt_files, ff.file.paths$Mt_files, Gtmat, Wtmat,
                                                            p = p, c = c,
                                                            out.head = bs.out.head, out.format = out.format, fst.head = fst.head)
                     else BS.bigdata(ff.file.paths$mt_files, ff.file.paths$Mt_files, Gtmat, Wtmat,
                                                            nT, DT, dT, R,
                                                            p = p, c = c,
                                                            L = n_samples, seed_list = seed_list,
                                                            out.head = bs.out.head, out.format = out.format, fst.head = fst.head,
                                                            verbose=verbose)

    #TODO: return also the Sigma samples for smoothing.
    return(list("FF.files" = ff.file.paths, "BS.files" = BS.file.paths))
}

#'
#' @param tau2s A T-length vector of values of tau2_t, t = 1:T. 

#TODO: document this function precisely

#' The FFBS implemented specifically for calibration.

#TODO: Write an SEXP version of this in C.

FFBS_calibration <- function(zdiffmat, m0z, M0ztau2_0, tau2s, U, save.dir = NULL) {
    #TODO: Split FF_calibration and BS_calibration into functions into FF.R and BS.R respectively?
    T <- nrow(zdiffmat)
    S <- ncol(zdiffmat)

    # initialize Atildemat
#    Atildemat <- matrix(nrow=T, ncol=(S*(S+1)/2), data=0)
#    Qtildemat <- matrix(nrow=T, ncol=(S*(S+1)/2), data=0)

    Atildemat <- matrix(nrow=T, ncol=S*S, data=0)
#    Qtildemat <- matrix(nrow=T, ncol=S*S, data=0)

    mzmat <- matrix(nrow=T+1, ncol=S, data=0)
    mzmat[1,] <- c(m0z)
#    tau2Mzmat <- matrix(nrow=T+1, ncol=(S*(S+1)/2), data=0)
#    tau2Mzmat[1,] <- M0ztau2_0[upper.tri(M0ztau2_0, diag=TRUE)]

    tau2Mzmat <- matrix(nrow=T+1, ncol=S*S, data=0)
    tau2Mzmat[1,] <- c(M0ztau2_0)

    diagS <- diag(S)

    Atildet <- diag(S)
    Qtildet <- diag(S)
    lQinvA <- diag(S)
    tau2tMzt <- diag(S)

    # FF
    for (t in 1:T) {
        #TODO: speed up with .C. BLAS and LAPACK.
        tau2tMztm1 <- matrix(tau2Mzmat[t,], nrow=S) #tri.to.sym(tau2Mzmat[t,], S, lower=FALSE)

        #TODO: Don't shape it yet if we want axpy to work???

        Atildet <- duplicate(tau2tMztm1, shallow=F) #+ tau2s[t] * U
#	.Call("copy", S*S, Atildet, tau2tMztm1)
	#TODO: debug. Why would axpy work for the first one but not for the second??? And why would it all come to the same answer for Qtildet???

	Atildet <- Atildet + tau2s[t] * U
#        .Call("axpy", S*S, tau2s[t], U, Atildet)

#        Atildemat[t,] <- Atildet[upper.tri(Atildet, diag=TRUE)]
        Atildemat[t,] <- c(Atildet)#[upper.tri(Atildet, diag=TRUE)]

#        Qtildet <- Atildet + tau2s[t] * diagS
	Qtildet <- duplicate(Atildet, shallow=F)
#	.Call("copy", S*S, Qtildet, Atildet)
	Qtildet <- Qtildet + tau2s[t] * diagS
#	.Call("axpy", S*S, tau2s[t], diagS, Qtildet)

#        Qtildemat[t,] <- Qtildet[upper.tri(Qtildet, diag=TRUE)]
#        Qtildemat[t,] <- c(Qtildet)

        #TODO: Remove later. We want to be sure our results are correct after we .Call("trsm", ...) for lQinvA below.
        uQtildet <- chol(Qtildet)

        zdifft <- zdiffmat[t,] - mzmat[t,]
        zdifft <- backsolve(uQtildet, backsolve(uQtildet, zdifft, transpose=TRUE))
#        zdifft <- .Call("solve_chol", Qtildet, S, zdifft, S)

#	.Call("trsm", 'L', 'U', 'T', 'N', S, 1, 1.0, uQtildet, S, zdifft, S)
#	.Call("trsm", 'L', 'U', 'N', 'N', S, 1, 1.0, uQtildet, S, zdifft, S)

	#TODO: better just to call solve_chol()
	#TODO: potrs doesn't exist?? I have it defined in lin_alg.c!!! 
#        .Call("potrs", "U", uQtildet, S, zdifft, S)

#        zdifft <- backsolve(uQtildet, zdifft, transpose=TRUE)
#        zdifft <- backsolve(uQtildet, zdifft)

#        lQinvA <- backsolve(uQtildet, Atildet, transpose=TRUE)
        lQinvA <- duplicate(Atildet, shallow=F)

#        .Call("copy", S*S, lQinvA, Atildet)	
	lQinvA <- backsolve(uQtildet, lQinvA, transpose=TRUE)
#	.Call("trsm", "L", "U", "T", "N", S, S, 1.0, uQtildet, S, lQinvA, S)

        mzmat[t+1,] <- mzmat[t,] + zdifft %*% Atildet 

#        .Call("gemv", 'N', S, S, 1.0, Atildet, S, zdifft, 1.0, mzmat[t+1,])

#        tau2tMzt <- duplicate(Atildet, shallow=F)
#        .Call("copy", S*S, tau2tMzt, Atildet)

        tau2tMzt <- Atildet - crossprod(lQinvA)
#        .Call("syrk", 'U', 'T', S, S, -1.0, lQinvA, S, 1.0, tau2tMzt, S)

        tau2Mzmat[t+1,] <- c(tau2tMzt) #tau2tMzt[upper.tri(tau2tMzt, diag=TRUE)]
    }

    htildemat <- matrix(nrow=T+1, ncol=S, data=0)
#    Htildemat <- matrix(nrow=T+1, ncol=(S*(S+1)/2), data=0)
    Htildemat <- matrix(nrow=T+1, ncol=S*S, data=0)
    htildemat[T+1,] <- mzmat[T+1,]
    Htildemat[T+1,] <- tau2Mzmat[T+1,]

    umat <- matrix(nrow=T+1, ncol=S, data=0)
#    umat[T+1,] <- rmvnorm(1, htildemat[T+1], tri.to.sym(Htildemat[T+1,], S, lower=FALSE))
    umat[T+1,] <- rmvnorm(1, htildemat[T+1], matrix(Htildemat[T+1,], nrow=S))

    # BS
    for (t in (T-1):0) {
#        tau2tMzt <- tri.to.sym(tau2Mzmat[t+1,], S, lower=FALSE)
#        Atildetpl1 <- tri.to.sym(Atildemat[t+1,], S, lower=FALSE)
        tau2tMzt <- matrix(tau2Mzmat[t+1,], nrow=S)
        Atildetpl1 <- matrix(Atildemat[t+1,], nrow=S)

        hdiff <- htildemat[t+2,] - mzmat[t+1,]
        Htildetpl1 <- matrix(Htildemat[t+2,], nrow=S)
	Adiff <- duplicate(Atildetpl1, shallow=F)
#        .Call("axpy", S*S, -1.0, Htildetpl1, Adiff)
        Adiff <- Adiff - Htildetpl1

        uAtildetpl1 <- chol(Atildetpl1)

        #TODO: substitute for dpotrs
        hdiff <- backsolve(uAtildetpl1, hdiff, transpose=TRUE)
        hdiff <- backsolve(uAtildetpl1, hdiff)
        
#	.Call("trsm", 'L', 'U', 'T', 'N', S, 1, 1.0, uAtildetpl1, S, hdiff, S)
#	.Call("trsm", 'L', 'U', 'N', 'N', S, 1, 1.0, uAtildetpl1, S, hdiff, S)

        htildet <- mzmat[t+1,] + hdiff %*% tau2tMzt
        htildemat[t+1,] <- c(htildet)

#        Htildetpl1 <- tri.to.sym(Htildemat[t+2,], S, lower=FALSE)
        lAdiff <- t(chol(Adiff))
        lAdiff <- backsolve(uAtildetpl1, lAdiff, transpose=TRUE)
        lAdiff <- backsolve(uAtildetpl1, lAdiff)

#	.Call("trsm", 'L', 'U', 'T', 'N', S, S, 1.0, uAtildetpl1, S, lAdiff, S)
#	.Call("trsm", 'L', 'U', 'N', 'N', S, S, 1.0, uAtildetpl1, S, lAdiff, S)

        lAdiff <- tau2tMzt %*% lAdiff
        Htildet <- tau2tMzt - tcrossprod(lAdiff)
       
#        Htildet <- tau2tMzt
#        .Call("syrk", 'U', 'T', S, S, -1.0, lAdiff, S, 1.0, Htildet, S)

	#TODO: syrk didn't work above for the FF, but it works here??? 
#        .C("copytri", A=Htildet, n_=S, uplo='U', PACKAGE="MNIW")$A

#        Htildemat[t+1,] <- Htildet[upper.tri(Htildet, diag=TRUE)]
        Htildemat[t+1,] <- c(Htildet)

        umat[t+1,] <- c(rmvnorm(1, htildet, Htildet))
    }

    if (!is.null(save.dir)) {
        #TODO: save FF and BS moment matrices.
    }
    
    return(umat)
}

