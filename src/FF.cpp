#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List FF_1step_cpp(const Eigen::MatrixXd& Yt, 
                        const Eigen::MatrixXd& Ft, 
                        const Eigen::MatrixXd& Gt, 
                        const Eigen::MatrixXd& Wt, 
                        const Eigen::MatrixXd& Vt, 
                        const Eigen::MatrixXd& mt_1, 
                        const Eigen::MatrixXd& Mt_1, 
                        double nt_1, 
                        const Eigen::MatrixXd& Dt_1, 
                        double delta = 1.0) {
  
  int N = Yt.rows();
  
  Eigen::MatrixXd at = Gt * mt_1;
  Eigen::MatrixXd At = Gt * Mt_1 * Gt.transpose() + Wt;
  Eigen::MatrixXd qt = Ft * at;
  Eigen::MatrixXd FtAt = Ft * At;
  Eigen::MatrixXd Qt = FtAt * Ft.transpose() + Vt;
  Eigen::MatrixXd Qtinv = Qt.inverse();
  Eigen::MatrixXd AttFtQtinv = At * Ft.transpose() * Qtinv;
  Eigen::MatrixXd Yt_qt = Yt - qt;
  Eigen::MatrixXd mt = at + AttFtQtinv * Yt_qt;
  Eigen::MatrixXd Mt = At - AttFtQtinv * FtAt;
  double nt = nt_1 + N;
  Eigen::MatrixXd Dt = Dt_1 + Yt_qt.transpose() * Qtinv * Yt_qt;
  
  return Rcpp::List::create(Rcpp::Named("nt") = nt,
                            Rcpp::Named("Dt") = Dt,                            
                            Rcpp::Named("at") = at,
                            Rcpp::Named("At") = At,
                            Rcpp::Named("mt") = mt,
                            Rcpp::Named("Mt") = Mt);
}
