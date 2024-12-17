library(matrixcalc)

#' Checks whether a matrix is positive definite symmetric, within a given tolerance.
#' @param C The matrix to be evaluated
#' @param tol The tolerance to determine if C is symmetric. If C - t(C) is greater than or equal to this tolerance, an error is thrown. Otherwise, C is set to the average of C and t(C).
#' @returns The matrix after forcing its symmetry.
#' @export
check.pds <- function(C, tol = 10^-9) {
  if (norm(C - t(C), type="F") < tol) {
	  C <- (C + t(C))/2
          if (!is.positive.definite(C)) {
		  cat("Warning: Matrix is not positive definite. Listing eigenvalues:", eigen(C)$values, "\n")
	  }
	  return(C)
  } else {
	  throw("Matrix symmetry is outside of tolerance: ",tol)
  }
}

#' Turns the array of the entries of a triangular matrix into a symmetric matrix and returns said matrix.
#' @param M.tri The array of triangular values, including the diagonal entries.
#' @param p The dimension of the full matrix.
#' @param lower Whether M.tri represents the lower- or upper-triangular values. Defaults to FALSE.
#' @returns M, the full symmetric matrix comprised of the entries of M.tri.

tri.to.sym <- function(M.tri, p, lower=FALSE) {
    M <- matrix(0, p, p)
    if (lower) {
        M[lower.tri(M, diag=TRUE)] <- M.tri
        M[upper.tri(M)] <- t(M)[upper.tri(M)]
    } else {
        M[upper.tri(M, diag=TRUE)] <- M.tri
        M[lower.tri(M)] <- t(M)[lower.tri(M)]
    }
    return(M)
}

