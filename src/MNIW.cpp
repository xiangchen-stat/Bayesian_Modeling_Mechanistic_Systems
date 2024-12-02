#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers


// // [[Rcpp::export]]
// Rcpp::List MNIW_naiive_cpp(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, 
//                       const Eigen::MatrixXd& H, double v,                        
//                       const Eigen::MatrixXd& S, const Eigen::MatrixXd& C, 
//                       const Eigen::MatrixXd& Vb) {
  
//   int n = Y.rows();
//     // int q = Y.cols(); // q is calculated but not used in the given function
//     // int p = X.cols(); // p is calculated but not used in the given function

//     Eigen::MatrixXd Vbinv = Vb.inverse();
//     Eigen::MatrixXd Hinv = H.inverse();
//     Eigen::MatrixXd XtHinv = X.transpose() * Hinv;
//     Eigen::MatrixXd Oinv = Vbinv + XtHinv * X;
//     Eigen::MatrixXd O = Oinv.inverse();
//     Eigen::MatrixXd VbinvC = Vbinv * C;
//     Eigen::MatrixXd M = VbinvC + XtHinv * Y;
//     Eigen::MatrixXd Snew = S + C.transpose() * VbinvC + Y.transpose() * Hinv * Y - M.transpose() * O * M;
//     double vnew = v + n;
//     Eigen::MatrixXd Cnew = O * M;
//     Eigen::MatrixXd Vbnew = O;
  
//   // Return the results
//   Rcpp::List output = Rcpp::List::create(
//     Rcpp::Named("vnew") = vnew, 
//     Rcpp::Named("Snew") = Snew,
//     Rcpp::Named("Cnew") = Cnew, 
//     Rcpp::Named("Vbnew") = Vbnew);
//   return output;
// }



// [[Rcpp::export]]
Rcpp::List MNIW_cpp(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, 
                           const Eigen::MatrixXd& H, double v,                        
                           const Eigen::MatrixXd& S, const Eigen::MatrixXd& C, 
                           const Eigen::MatrixXd& Vb) {
  
  int n = Y.rows();
  int q = Y.cols();
  int p = X.cols();
  
  Eigen::MatrixXd Vbinv = Vb.inverse();
  // Eigen::MatrixXd Vbinv = Vb.llt().solve(Eigen::MatrixXd::Identity(p, p));
  Eigen::MatrixXd Hinv = H.inverse();
  // Eigen::MatrixXd Hinv = H.llt().solve(Eigen::MatrixXd::Identity(n, n));
  Eigen::MatrixXd XtHinv = X.transpose() * Hinv;
  Eigen::MatrixXd Oinv = Vbinv + XtHinv * X;
  Eigen::MatrixXd O = Oinv.inverse();
  // Eigen::MatrixXd O = Oinv.llt().solve(Eigen::MatrixXd::Identity(p, p));
  Eigen::MatrixXd VbinvC = Vbinv * C;
  Eigen::MatrixXd M = VbinvC + XtHinv * Y;
  Eigen::MatrixXd Snew = S + C.transpose() * VbinvC + Y.transpose() * Hinv * Y - M.transpose() * O * M;
  double vnew = v + n;
  Eigen::MatrixXd Cnew = O * M;
  Eigen::MatrixXd Vbnew = O;
  
  // Return the results
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("vnew") = vnew, 
    Rcpp::Named("Snew") = Snew,
    Rcpp::Named("Cnew") = Cnew, 
    Rcpp::Named("Vbnew") = Vbnew);
  return output;
}



