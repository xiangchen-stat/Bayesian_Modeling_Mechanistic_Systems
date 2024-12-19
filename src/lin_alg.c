#define USE_FC_LEN_T

#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>
#include "lin_alg.h"

#ifndef FCONE
# define FCONE
#endif

//Generates a random number in N(0,1).
//  Taken from: https://kcru.lawsonresearch.ca/research/srk/normalDBN_random.html#:~:text=The%20rand%20%28%29%20call%20returns%20a%20random%20number,to%20double%20to%20make%20the%20division%20floating%20point.
double randn ()
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return ((double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  //X2 = U2 * mult;

  call = !call;

  return ((double) X1);
}

// Copies X_ to Y_
SEXP copy(SEXP n_, SEXP Y_, SEXP X_) {
    int one = 1;
    int n = asInteger(n_);
    double* X = REAL(X_);
    double* Y = REAL(Y_);

    F77_CALL(dcopy)(&n, X, &one, Y, &one);

    return(Y_);
}


SEXP axpy(SEXP n_, SEXP A_, SEXP X_, SEXP Y_) {

    int n = asInteger(n_);
    double A = asReal(A_);
    double* X = REAL(X_);
    double* Y = REAL(Y_);

    int one = 1;

    F77_CALL(daxpy)(&n, &A, X, &one, Y, &one);

    return(Y_);
}

SEXP gemv(SEXP transa_, SEXP m_, SEXP n_, SEXP alpha_,
		SEXP A_, SEXP ldA_, SEXP x_, SEXP beta_, SEXP y_) {

    const char * transa = CHAR(asChar(transa_));
    
    int m = asInteger(m_);
    int n = asInteger(n_);
    double alpha = asReal(alpha_);
    double beta = asReal(beta_);

    double* A = REAL(A_);
    int ldA = asInteger(ldA_);

    double* x = REAL(x_);
    double* y = REAL(y_);

    int one = 1;

    F77_CALL(dgemv)(transa, &m, &n, &alpha, A, &ldA,
		    x, &one, &beta, y, &one FCONE);

    return(y_);
}

SEXP gemm(SEXP transa_, SEXP transb_, SEXP m_, SEXP n_, SEXP k_,
		SEXP alpha_, SEXP A_, SEXP ldA_, SEXP B_, SEXP ldB_,
		SEXP beta_, SEXP C_, SEXP ldC_) {

    const char * transa = CHAR(asChar(transa_));
    const char * transb = CHAR(asChar(transb_));
    
    int m = asInteger(m_);
    int n = asInteger(n_);
    int k = asInteger(k_);
    double alpha = asReal(alpha_);
    double beta = asReal(beta_);

    double* A = REAL(A_);
    int ldA = asInteger(ldA_);

    double* B = REAL(B_);
    int ldB = asInteger(ldB_);

    double* C = REAL(C_);
    int ldC = asInteger(ldC_);

    F77_CALL(dgemm)(transa, transb, &m, &n, &k, &alpha, A, &ldA, B, &ldB,
		    &beta, C, &ldC FCONE FCONE);

    return(C_);
}


SEXP trmm(SEXP side_, SEXP uplo_, SEXP transa_, SEXP diag_, SEXP m_, SEXP n_,
		SEXP alpha_, SEXP A_, SEXP ldA_, SEXP B_, SEXP ldB_) {
    const char * side = CHAR(asChar(side_));
    const char * uplo = CHAR(asChar(uplo_));
    const char * transa = CHAR(asChar(transa_));
    const char * diag = CHAR(asChar(diag_));

    int m = asInteger(m_);
    int n = asInteger(n_);
    double alpha = asReal(alpha_);

    double* A = REAL(A_);
    int ldA = asInteger(ldA_);

    double* B = REAL(B_);
    int ldB = asInteger(ldB_);

    F77_CALL(dtrmm)(side, uplo, transa, diag, &m, &n, &alpha, A, &ldA, B, &ldB
		    FCONE FCONE FCONE FCONE);

    return(B_);
}

SEXP trsm(SEXP side_, SEXP uplo_, SEXP transa_, SEXP diag_, SEXP m_, SEXP n_,
		SEXP alpha_, SEXP A_, SEXP ldA_, SEXP B_, SEXP ldB_) {

    const char * side = CHAR(asChar(side_));
    const char * uplo = CHAR(asChar(uplo_));
    const char * transa = CHAR(asChar(transa_));
    const char * diag = CHAR(asChar(diag_));

    int m = asInteger(m_);
    int n = asInteger(n_);
    double alpha = asReal(alpha_);

    double* A = REAL(A_);
    int ldA = asInteger(ldA_);

    double* B = REAL(B_);
    int ldB = asInteger(ldB_);
    
    F77_CALL(dtrsm)(side, uplo, transa, diag, &m, &n, &alpha, A, &ldA, B, &ldB
		    FCONE FCONE FCONE FCONE);

    return(B_);
}

SEXP syrk(SEXP uplo_, SEXP trans_, SEXP n_, SEXP k_, 
		SEXP alpha_, SEXP A_, SEXP ldA_, SEXP beta_, SEXP C_, SEXP ldC_) {
    const char * uplo = CHAR(asChar(uplo_));
    const char * trans = CHAR(asChar(trans_));

    int n = asInteger(n_);
    int k = asInteger(k_);

    double alpha = asReal(alpha_);
    double* A = REAL(A_);
    int ldA = asInteger(ldA_);

    double beta = asReal(beta_);
    double* C = REAL(C_);
    int ldC = asInteger(ldC_);

    F77_CALL(dsyrk)(uplo, trans, &n, &k, &alpha, A, &ldA,
		    &beta, C, &ldC FCONE FCONE);

    return(C_);
}

// Solve AX = B, where A is factored using the Cholesky decomposition outside of the loop.
SEXP potrs(SEXP uplo_, SEXP A_, SEXP n_, SEXP B_, SEXP q_)
{
    const char * uplo = CHAR(asChar(uplo_));
    int info;
    int n = asInteger(n_);
    int q = asInteger(q_);
    double* A = REAL(A_);
    double* B = REAL(B_);

    F77_CALL(dpotrs)(uplo, &n, &q, A, &n, B, &n, &info FCONE);

    return(B_);
}

// Generates a matrix random normal variable in-place.
void rmn(double *X, int *n, int *q, double *mu, double *u_V, double *u_S) {
    int in = n[0];
    int iq = q[0];

    for(int i = 0; i < in * iq; i++) {
        // in-place samples N(0,1) variates for X
        X[i] = randn();
    }

    // Left- and right-multiply.
    double one = 1.0;
    
    char* right = "R";
    char* left = "L";
    char* up = "U";
    char* trans = "T";
    char* no = "N";

    F77_CALL(dtrmm)(right, up, no, no, n, q, &one, u_S, q, X, n
		    FCONE FCONE FCONE FCONE);
    F77_CALL(dtrmm)(left, up, trans, no, n, q, &one, u_V, n, X, n
		    FCONE FCONE FCONE FCONE);

    for(int i = 0; i < in * iq; i++) {
        X[i] += mu[i];
    }
}

// Copies the upper-triangular matrix to a lower-triangular matrix or vice-versa,
//   where the matrix A is the flattened column-stacked array.
// uplo determines which segment of the matrix to copy across the diagonal. If 
//   uplo == 'U', then copytri copies the upper-triangular part of the matrix to the
//   lower-triangle. If uplo == 'D', copytri copies the lower-triangular part of the
//   matrix to the upper-triangular part.
   
void copytri(double *A, int *n_, char uplo) {
    int n = n_[0];
    // check that A is n^2
    size_t nA = sizeof(A)/sizeof(A[0]);
    if (nA != n*n) {
        printf("Error: argument matrix A to copytri() is not %d x %d.\n", n, n);
        exit(1);
    }
    
    //TODO: Does the copy actually do what it should be doing??? Test it (you can use another language)!
    if (uplo == 'U') {
	for (int i=0; i<n; i++) {
            for (int j=i+1; j<n; j++) {
                A[i*n+j] = A[j*n+i];
	    }
	}
    } else if (uplo == 'D') {
	for (int i=0; i<n; i++) {
            for (int j=i+1; j<n; j++) {
                A[j*n+i] = A[i*n+j];
	    }
	}	
    } else {
        printf("Error: argument uplo to copytri() must be either 'U' or 'D'.");
	exit(1);
    }
}


//// SEXP version of copytri
//
//SEXP copytri(SEXP A_, SEXP n_, SEXP uplo_) {
//    const char * uplo = CHAR(asChar(uplo_));
//
//    int n = asInteger(n_);
//
//
//}


// Replaces the upper-triangular entries of an n x n PDS matrix A with its
//   Cholesky factor. 
void chol(double *A, int *n) {
    int info;
    F77_CALL(dpotrf)("U", n, A, n, &info FCONE);
    if (info != 0)
        error("Cholesky decomposition failed.");
}


// Replace the inverse of an n x n PDS matrix V via Cholesky decomposition.
// Code sourced from: https://github.com/cjgeyer/mat/blob/master/package/baz/src/i.c

void chol2inv(double *A, int *n)
{
    int info;

    F77_CALL(dpotrf)("L", n, A, n, &info FCONE);
    if (info != 0)
        error("Cholesky decomposition failed.");

    F77_CALL(dpotri)("L", n, A, n, &info FCONE);
    if (info != 0)
        error("Cholesky inversion failed.");

    // Copy the lower-triangle to the upper-triangle.
    // Note that the entries of a matrix are stored row-wise.
    int in = n[0];
    for(int i = 0; i < in; i++) {
        for (int j = 0; j < i; j++) {
            A[j + in * i] = A[i + in * j];
        }
    }
}

// Solve AX = B, where A is Positive-Symmetric Definite, given the Cholesky decomposition of A.
// Essentially the same as computing potrf(), and then potrs().
SEXP solve_chol(SEXP A_, SEXP n_, SEXP B_, SEXP q_)
{
    int info;
    int n = asInteger(n_);
    int q = asInteger(q_);
    double* A = REAL(A_);
    double* B = REAL(B_);

    F77_CALL(dpotrf)("L", &n, A, &n, &info FCONE);
    if (info != 0)
        error("Cholesky decomposition failed.");

    F77_CALL(dpotrs)("L", &n, &q, A, &n, B, &n, &info FCONE);

    return(B_);
}

