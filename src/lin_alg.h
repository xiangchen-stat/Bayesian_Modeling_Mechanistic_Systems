#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

#ifdef __cplusplus
extern "C" {
#endif
    SEXP copy(SEXP n_, SEXP Y_, SEXP X_);
    
    SEXP axpy(SEXP n_, SEXP A_, SEXP X_, SEXP Y_);
    
    SEXP gemv(SEXP transa_, SEXP m_, SEXP n_, SEXP alpha_, SEXP A_, SEXP ldA_, SEXP x_, SEXP beta_, SEXP y_);
    
    SEXP gemm(SEXP transa_, SEXP transb_, SEXP m_, SEXP n_, SEXP k_,
              SEXP alpha_, SEXP A_, SEXP ldA_, SEXP B_, SEXP ldB_,
              SEXP beta_, SEXP C_, SEXP ldC_);
    
    SEXP trmm(SEXP side_, SEXP uplo_, SEXP transa_, SEXP diag_, SEXP m_, SEXP n_,
              SEXP alpha_, SEXP A_, SEXP ldA_, SEXP B_, SEXP ldB_);
    
    SEXP trsm(SEXP side_, SEXP uplo_, SEXP transa_, SEXP diag_, SEXP m_, SEXP n_,
              SEXP alpha_, SEXP A_, SEXP ldA_, SEXP B_, SEXP ldB_);
    
    SEXP syrk(SEXP uplo_, SEXP trans_, SEXP n_, SEXP k_, 
              SEXP alpha_, SEXP A_, SEXP ldA_, SEXP beta_, SEXP C_, SEXP ldC_);
    
    SEXP potrs(SEXP uplo_, SEXP A_, SEXP n_, SEXP B_, SEXP q_);
    
    SEXP solve_chol(SEXP A_, SEXP n_, SEXP B_, SEXP q_);

#ifdef __cplusplus
}
#endif

