#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;

// Function to generate a matrix normal distribution
// [[Rcpp::export]]
Eigen::MatrixXd rmn_cpp(const Eigen::MatrixXd& m, const Eigen::MatrixXd& U, const Eigen::MatrixXd& V) {
    int p = m.rows(); // Number of rows in m
    int S = m.cols(); // Number of cols in m
    Eigen::MatrixXd Z = Eigen::MatrixXd::NullaryExpr(p, S, [](){ return R::rnorm(0, 1); });
    Eigen::LLT<Eigen::MatrixXd> cholU(U);  // Cholesky decomposition of U
    Eigen::LLT<Eigen::MatrixXd> cholV(V);  // Cholesky decomposition of V
    return m + cholU.matrixL() * Z * cholV.matrixL().transpose();  // Construct the matrix normal distribution
}

// [[Rcpp::export]]
Eigen::MatrixXd rmn_chol_cpp(const Eigen::MatrixXd& m, const Eigen::MatrixXd& RM, const Eigen::MatrixXd& RSigma) {
    int p = m.rows(); // Number of rows in m
    int S = m.cols(); // Number of cols in m

    // Generate a matrix of standard normal random variables
    Eigen::MatrixXd Z = Eigen::MatrixXd::NullaryExpr(p, S, [](){ return R::rnorm(0, 1); });

    // Calculate the result
    Eigen::MatrixXd theta = m + RM.transpose() * Z * RSigma;

    return theta;
}