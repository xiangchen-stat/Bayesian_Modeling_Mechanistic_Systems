library(matrixsampling)
library(msos)
library(HypergeoMat)

#' Sample from the Hyper-T distribution.
#' @param n Number of samples
#' @param m The mean of the matrix-T random variable.
#' @param M The p x p row covariance matrix of the matrix-T random variable.
#' @param nu Degrees of freedom of the inverse-Wishart.
#' @param D The S x S scale matrix of inverse-Wishart
#' @returns A 3D array of matrix-T-distributed samples, with the third dimension being the index for each individual sample.
#' @export

#TODO: offload this to C, especially since the for loop is really bad for its runtime.
rmatrixt <- function(n, m, M, nu, D){
    p <- nrow(m)
    S <- ncol(m)

    Sigmas <- rinvwishart(n, nu=nu, Omega=D, epsilon=0, checkSymmetry=F)

    B <- array(dim = c(p, S, n))
    for (i in 1:n) {
        B[,,i] <- rmn(1, m, M, Sigmas[,,i])[,,1]
    }
    return(B)
}


#' Computes the density of a matrix-T variable. Defaults to the log density.
#' @param q The matrix to compute the matrix-T.
#' @param m The mean of the matrix-T random variable.
#' @param M The p x p row covariance matrix of the matrix-T random variable.
#' @param nu Degrees of freedom of the inverse-Wishart.
#' @param D The S x S scale matrix of inverse-Wishart
#' @param log Whether to output the log-density or the density. Defaults to TRUE (i.e., outputs log-density).
#' @returns The log-density or density.

#TODO: How to do this efficiently for a stack of matrices q?

dmatrixt <- function(q, m, M, nu, D, log=TRUE){
    p <- nrow(m)
    S <- ncol(m)

    uM <- chol(M)
    uD <- chol(D)
    Qbig <- q - m
    Qbig <- backsolve(uM, Qbig, transpose=TRUE)
    Qbig <- crossprod(Qbig)
    Qbig <- backsolve(uD, backsolve(uD, Qbig, transpose=TRUE))
    Qbig <- Qbig + diag(S)

    logp <- -p*S/2 * log(pi) + lmvgamma((nu + p)/2, S) - lmvgamma(nu/2, S) -
	    S/2 * logdet(M) - p/2 * logdet(D) - (nu + p)/2 * logdet(Qbig)

    if (!log) return(exp(logp))
    return(logp)
}


#' Sample from the Hyper-T distribution where Sigma = sigma2 * R, with sigma2 being distributed as an inverse-Gamma with paramaters (nu, d) and R being a scale matrix.
#' @param n Number of samples
#' @param m The mean of the matrix-T random variable.
#' @param M The p x p row covariance matrix of the matrix-T random variable.
#' @param nu Degrees of freedom of the inverse-Wishart.
#' @param d The scale of the inverse gamma sigma2.
#' @param R The scale column-covariance matrix.
#' @returns List of mean and variance samples from the matrix-normal inverse-Wishart distribution.
#' @export

#TODO: rename this???
rmatrixtscS <- function(n, m, M, nu, d, R){
    p <- nrow(m)
    S <- ncol(m)

    sigma2s <- 1/rgamma(n, nu, rate=d)

    B <- array(dim = c(p, S, n))
    # Y <- array(dim = c(n, q, nsam))
    for (i in 1:n) {
        B[,,i] <- rmn(1, m, M, sigma2s[i] * R)[,,1]
    }
    return(B)
}

#' Computes the density of a matrix-T variable where Sigma = sigma2 * R, with sigma2 being distributed as an inverse-Gamma with parameters (nu, d) and R being a scale matrix. Defaults to the log density.
#' @param q The matrix to compute the matrix-T.
#' @param m The mean of the matrix-T random variable.
#' @param M The p x p row covariance matrix of the matrix-T random variable.
#' @param nu Degrees of freedom of the inverse-Wishart.
#' @param d The inverse-gamma scale parameter. Corresponds to the rate in the gamma.
#' @param R The S x S scale matrix.
#' @param log Whether to output the log-density or the density. Defaults to TRUE (i.e., outputs log-density).
#' @returns The log-density or density.

dmatrixtscS <- function(q, m, M, nu, d, R, log=TRUE) {
    p <- nrow(m)
    S <- ncol(m)

    uM <- chol(M)
    uR <- chol(R)
    Qbig <- q - m
    Qbig <- backsolve(uM, Qbig, transpose=TRUE)
    Qbig <- backsolve(uR, t(Qbig), transpose=TRUE)

    logp <- -(p * S/2 + nu) * log(1 + sum(Qbig * Qbig)/(2 * d)) + lmvgamma(p*S/2 + nu, 1) - lmvgamma(nu, 1) - 
	    p*S/2 * log(2 * pi * d) - S/2 * logdet(M) - p/2 * logdet(R)

    if (!log) return(exp(logp))
    return(logp)
}

