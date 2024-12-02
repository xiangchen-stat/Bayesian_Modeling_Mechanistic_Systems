#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List BS_1step_cpp(const Eigen::MatrixXd& mt, 
                      const Eigen::MatrixXd& Mt, 
                      const Eigen::MatrixXd& st1, 
                      const Eigen::MatrixXd& St1, 
                      const Eigen::MatrixXd& at1, 
                      const Eigen::MatrixXd& At1, 
                      const Eigen::MatrixXd& Gt1, 
                      double delta = 1.0) {
  
  Eigen::MatrixXd At1inv = At1.inverse();
  Eigen::MatrixXd MttGt1At1inv = Mt * Gt1.transpose() * At1inv;
  Eigen::MatrixXd st = mt + MttGt1At1inv * (st1 - at1);
  Eigen::MatrixXd St = Mt - MttGt1At1inv * (At1 - St1) * MttGt1At1inv.transpose();
  
  return Rcpp::List::create(Rcpp::Named("st") = st,
                            Rcpp::Named("St") = St);
}

