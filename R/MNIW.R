library(msos)

#' samples a matrix normal variable from MN(mu, V, S)
#' @param N The number of random matrix normals to generate, stacked along the third dimension of the array.
#' @param mu The mean matrix of the matrix-normal
#' @param V The covariance matrix for the rows
#' @param S The covariance matrix for the columns
#' @returns A random matrix normal from MN(mu, V, S)
#' @export
rmn <- function(N, mu, V=diag(nrow(mu)), S=diag(ncol(mu))) {
    n <- nrow(mu)
    q <- ncol(mu)

    # preprocessing
    if (is.null(n)) {
        n <- length(mu)
        if (is.null(V)) V <- diag(n)
    }
    if (is.null(q)) {
        q <- 1
        if (is.null(S)) S <- as.matrix(1)
    }

    u_V <- chol(V)
    u_S <- chol(S)

    #TODO: How are the slices of this array populated? In what order?
    X <- array(rnorm(n*q*N), dim=c(n,q,N))

#    X <- X %*% u_S

    for (i in 1:N) {
        X[,,i] <- X[,,i] %*% u_S
#        X[,,i] <- .Call("trmm", "R", "U", "N", "N", n, q, 1.0, u_S, q, X[,,i], n)
#    .C("rmn", X=X, n=n, q=q, mu=mu, u_V=u_V, u_S=u_S, PACKAGE="MNIW")$X

        X[,,i] <- crossprod(u_V, X[,,i])
#        X[,,i] <- .Call("trmm", "L", "U", "T", "N", n, q, 1.0, u_V, n, X[,,i], n)
	X[,,i] <- X[,,i] + mu

    }
    return(X)
}

#' The density of the matrix normal.
#' @param y The n x p matrix over which to compute the value of the density.
#' @param mu The n x p mean matrix. Defaults to the zero matrix when a default is absent.
#' @param V The p x p row-covariance matrix. Defaults to the identity.
#' @param S The n x n column-covariance matrix. Defaults to the identity.
#' @param log Whether to return the log-density at y or the density. Defaults to TRUE.
#' @returns The log-density of the matrix normal at y, or the density if log is set to false.

dmn <- function(y, mu = matrix(0, nrow=nrow(y), ncol=ncol(y)), 
		V = diag(nrow(y)), S = diag(ncol(y)), log=TRUE) {
    p <- ncol(mu)
    n <- nrow(mu)
    
    uV <- chol(V)
    uS <- chol(S)
    qdiff <- y - mu

    qdiff <- .Call("trsm", "L", "U", "T", "N", n, p, 1.0, uV, n, qdiff, n)
    qdiff <- .Call("trsm", "R", "U", "N", "N", n, p, 1.0, uS, p, qdiff, n)

    # qdiff * qdiff is pointwise product. It is the same as computing `tr(t(qdiff) %*% qdiff)`.
    logp <- -(p*n*log(2 * pi) + p * logdet(V) + n * logdet(S) +
  	      sum(qdiff*qdiff))/2

    if (!log) return(exp(logp))
    return(logp)
}


#' matrix-normal posterior update of Theta with data "point" matrix Y.
#' @param m0 The p x q mean prior of Theta.
#' @param M0 The p x p row-covariance matrix prior of Theta.
#' @param Sigma The q x q column-covariance matrix of Theta. Sigma is sampled via an Inverse Wishart distribution, but is included as an input here.
#' @param Y An n x q matrix which is assumed to follow from the linear equation Y = X * Theta + E, where E ~ MN(0, V, Sigma).
#' @param X An n x p matrix of covariates.
#' @param V The p x p row-covariance matrix of E.
#' @returns The matrix normal mean and row-covariance matrix for Theta | Sigma, Y: (m, M)
#' @export
Theta_post <- function(m0, M0, Sigma, Y, X, V) {
    u_M0 <- chol(M0)
    M <- u_M0 %*% t(X)

    #TODO: syrk into a storage matrix?
    M <- crossprod(M) + V
    u_Mint <- chol(M)
    M <- X %*% M0
    
    n <- nrow(X)
    p <- ncol(X)

    M <- .Call("trsm", "L", "U", "T", "N", n, p, 1.0, u_Mint, n, M, n)
#    M <- backsolve(u_Mint, M, transpose=TRUE)
    M <- M0 - crossprod(M)

    u_V <- chol(V)
    m <- crossprod(X, backsolve(u_V, backsolve(u_V, Y, transpose=TRUE), transpose=FALSE)) +
        backsolve(u_M0, backsolve(u_M0, m0, transpose=TRUE), transpose=FALSE)

    m <- M %*% m

    return(list(m = m, M = M))
}

#' For the model Y = X * Theta + E, where E ~ MN(0, V, Sigma), Sigma is assumed to be inverse-Wishart distributed, with degrees of freedom nu0 and scale matrix Psi0.
#' @param nu0 The degrees of freedom of the prior for Sigma.
#' @param Psi0 The q x q prior scale matrix for Sigma.
#' @param Y The n x q data "point" matrix.
#' @param X The n x p matrix of covariates.
#' @param V The p x p row-covariance matrix of E.
#' @param m0 The p x q mean prior of Theta.
#' @param M0 The p x p row-covariance matrix prior of Theta.
#' @returns The parameters - the degrees of freedom nu and scale matrix Psi - of Sigma | Y.
#' @export
Sigma_post <- function(nu0, Psi0, Y, X, V, m0, M0) {
    p <- ncol(X)
    if (nu0 <= p - 1) stop("Degrees of freedom nu0 is smaller than p.")
    n <- nrow(X)

    theta_uds <- Theta_post(m0, M0, Sigma, Y, X, V)
    m <- theta_uds$m

    u_V <- chol(V)
    u_M0 <- chol(M0)

    m1 <- crossprod(X, backsolve(u_V, backsolve(u_V, Y, transpose=TRUE), transpose=FALSE)) +
        backsolve(u_M0, backsolve(u_M0, m0, transpose=TRUE), transpose=FALSE)

    Ysc <- backsolve(u_V, Y, transpose=TRUE)
    m0sc <- backsolve(u_M0, m0, transpose=TRUE)

    Psi <- crossprod(Ysc) + crossprod(m0sc) - crossprod(m1, m) + Psi0

    q <- nrow(Psi0)
    nu <- nu0 + n 
    
    return(list(nu = nu, Psi = Psi))
}

