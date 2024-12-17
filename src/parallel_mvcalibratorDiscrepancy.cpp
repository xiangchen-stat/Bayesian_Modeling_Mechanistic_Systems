// #define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#include <RcppDist.h>
#include <RcppThread.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppThread)]]

using namespace arma;

vec backwardVolatility(const vec& n, const vec& d, const double delta_v){
  uword T = n.n_elem;
  vec v = n;
  double vinv, v_star;
  vinv = R::rgamma(n(T-1), 1.0/d(T-1));
  v(T-1) = 1.0/vinv;
  for(int t=T-2; t>=0; t--){
    v_star = R::rgamma((1-delta_v)*n(t+1), 1.0/d(t));
    vinv = delta_v*vinv + v_star;
    v(t) = 1.0/vinv;
  }
  return(v);
}

mat backwardState(const mat& a, 
                  const mat& m, 
                  const cube& CC, 
                  const cube& RR,
                  const vec v,
                  const double nugget){
  uword T = a.n_cols;
  mat theta = m; 
  double vinv;
  cube C = CC, R = RR;
  for(int i=0; i<C.n_slices; i++){
    R.slice(i) = (v(i))*R.slice(i);
    C.slice(i) = (v(i))*C.slice(i);
  }
  mat s = m.col(T-1);
  mat S = C.slice(T-1);
  theta.col(T-1) = rmvnorm(1, s, S+nugget*eye(S.n_rows,S.n_cols)).t();
  mat Rinv = R.slice(0);
  for(int t=T-2; t>=0; t--){
    Rinv = inv(R.slice(t+1));
    s = m.col(t) + C.slice(t)*Rinv*(s-a.col(t+1));
    S = C.slice(t) + C.slice(t)*Rinv*(R.slice(t+1)-S)*Rinv*C.slice(t);
    theta.col(t) = rmvnorm(1, s, S+nugget*eye(S.n_rows,S.n_cols)).t();
  }
  return(theta);
}

std::tuple<mat, vec> forwardBackward(const mat& y, const vec& m0, const mat& C0,
                                     double a0, double b0, double c0,
                                     const cube& FF, const mat& G,
                                     const mat& Sigma,
                                     const vec& n0, const vec& d0, 
                                     const double& delta_w, const double& delta_v)
{
  int T = y.n_cols;
  int p = y.n_rows;
  int p2 = G.n_cols;
  
  mat theta(p2, T, fill::zeros);
  mat a(p2, T, fill::zeros);
  mat m(p2, T, fill::zeros);
  vec f(p, fill::zeros);
  mat error(p, T, fill::zeros);
  mat H(p2, p2,fill::eye);
  mat h(p2, 1, fill::ones);
  
  cube R(p2, p2, T, fill::ones);
  mat Q(p, p, fill::ones);
  cube C(p2, p2, T, fill::zeros);
  vec v_draw(T, fill::zeros);
  vec n(T, fill::zeros), d(T, fill::zeros), s(T, fill::zeros);
  // 
  a.col(0) =  m0.col(0);
  R.slice(0) = C0/delta_w;
  f = FF.slice(0) * a.col(0);
  Q = FF.slice(0) * R.slice(0) * FF.slice(0).t() + Sigma;
  mat Qinv = inv(Q);
  m.col(0) = a.col(0) + R.slice(0) * FF.slice(0).t() * Qinv * (y.col(0) - f.col(0));
  C.slice(0) = R.slice(0) - R.slice(0) * FF.slice(0).t() * Qinv * FF.slice(0) * R.slice(0);
  // modify alpha and beta for discount factor
  n(0) = as_scalar(delta_v*n0) + 0.5*y.n_rows;
  d(0) = as_scalar(delta_v*d0) + 0.5*dot(y.col(0) - f.col(0), Qinv * (y.col(0) - f.col(0)));
  for(uword t=1; t<T; t++){
    // # // predict
    a.col(t) =  m.col(t-1);
    R.slice(t) =  C.slice(t-1)/delta_w;
    // # // marginalize
    f = FF.slice(t) * a.col(t);
    Q = FF.slice(t) * R.slice(t) * FF.slice(t).t() + Sigma;
    Qinv = inv(Q);
    // # // discount factors
    n(t) = delta_v*n(t-1) + 0.5*y.n_rows;
    d(t) = delta_v*d(t-1) + 0.5*dot(y.col(t) - f, Qinv * (y.col(t) - f));
    // # // filter
    m.col(t) = a.col(t) + R.slice(t) * FF.slice(t).t() * Qinv * (y.col(t) - f);
    C.slice(t) = (R.slice(t) - R.slice(t) * FF.slice(t).t() * Qinv * FF.slice(t) * R.slice(t));
  }
  
  // mat m2 = join_rows(m0,m);
  // cube C2 = join_slices(C0,C);
  // // expressions from Bayesian Filtering and Smoothing book
  // theta.col(T) = rmvnorm(1, m2.col(T), C2.slice(T)).t();
  // double vinv;
  // for(int t=T-1; t>=0; t--){
  //   theta.col(t) = rmvnorm(1,
  //             (1-delta_w)*m2.col(t) + delta_w*theta.col(t+1),
  //             (1-delta_w)*C2.slice(t)).t();
  // }
  // for(int t=T-1; t>=0; t--){
  //   //vinv = R::rgamma(n(t)/2.0, 2.0/beta(t));
  //   v_draw(t) = 1/R::rgamma(n(t)/2.0, 2.0/d(t));
  // }
  v_draw = backwardVolatility(n, d, delta_v);
  theta = backwardState(a, m, C, R, v_draw, 0);
  // mat SSs;
  // double ssy;
  // uword k=0;
  // mat theta2 = theta;
  // for(int t=0; t<T; t++){
  //   k = t+1;
  //   ssy = dot(y.col(t) - FF.slice(t)*theta.col(t), y.col(t) - FF.slice(t)*theta.col(t));
  //   v_draw(t) = 1/R::rgamma(a0,1/ssy);
  //   SSs = (G*theta2.col(k) - G*theta2.col(k-1)) * (G*theta2.col(k) - G*theta2.col(k-1)).t()/2;
  //   W_draw.slice(t) = inv(rwish(b0+t/2.0, inv(diagmat(c0*ones(p2)) + SSs)));
  // }
  return std::make_tuple(theta, v_draw);
}

// eventually change this to suggestion by Higdon
double priorBeta(vec beta)
{
  vec out(beta.n_elem, fill::zeros);
  for(uword i=0; i < out.n_elem; i++){
    out(i) = R::dgamma(beta(i),1,.1,true);
  }
  return(sum(out));
}

double calPrior(rowvec calibrate){
  rowvec out(calibrate.n_elem, fill::zeros);
  for(uword i=0; i < out.n_elem; i++){
    out(i) = R::dbeta(calibrate(i),4,4,true);
  }
  return(sum(out));
}

double priorBeta2(vec beta){
  rowvec out(beta.n_elem, fill::zeros);
  rowvec rho(beta.n_elem, fill::zeros);
  rowvec post(beta.n_elem, fill::zeros);
  for(uword i=0; i < out.n_elem; i++){
    rho(i) = R::rbeta(1,.1);
    out(i) = -4*log(rho(i));
    post(i) = pow(1-exp(out(i)/4),.9)*exp(-out(i)/4);
  }
  return(log(prod(post)));
}


void kernel(const cube& y, const mat& params, const cube& theta, vec& beta,
            const field<cube>& FF, mat& K, const vec tune, const mat& v, 
            const int S, const double nugget)
{
  uword T = y.n_cols;
  //compute proposed kernel
  vec eta = log(beta);
  vec beta_star = exp(eta+rmvnorm(1, zeros(beta.n_rows,1), diagmat(tune)).t());
  mat K_star = eye(params.n_rows, params.n_rows);
  for(uword i=0; i < params.n_rows; i++) {
    for (uword j = (i+1); j< params.n_rows; j++){
      // squared exp
      K_star(i, j) = exp(-dot(beta_star.t() % (params.row(i) - params.row(j)),(params.row(i) - params.row(j))));
      // abs exponential
      // K_star(i, j) = exp(-sum(beta_star.t() % abs(params.row(i) - params.row(j))));
      K_star(j, i) = K_star(i, j);
    }
  }
  // compute likelihoods
  mat L(S, T, fill::zeros), L_star(S, T, fill::zeros);
  for(uword i=0; i < T; i++){
    for(uword sp=0; sp < S; sp++){
      L_star(sp,i) = dmvnorm(y.slice(sp).col(i).t(), FF(sp).slice(i)*theta.slice(sp).col(i), v(i,sp)*K_star, true)(0,0);
      L(sp, i) = dmvnorm(y.slice(sp).col(i).t(), FF(sp).slice(i)*theta.slice(sp).col(i), v(i, sp)*K, true)(0,0);
    }
  }
  // // compute MH ratio
  double logr = accu(L_star) + priorBeta(beta_star) + sum(log(1/beta))
    - accu(L) - priorBeta(beta) - sum(log(1/beta_star));
  if(log(R::runif(0,1)) < logr){
    beta = beta_star;
    K=K_star;
  }
  K = K + nugget*eye(K.n_rows, K.n_cols);
}

void kernel2(const mat& y, const mat& params, const mat& theta, vec& beta, mat& K, const vec tune, const double& v, 
             const int S, const double nugget)
{
  uword T = y.n_cols;
  //compute proposed kernel
  vec eta = log(beta);
  vec beta_star = exp(eta+rmvnorm(1, zeros(beta.n_rows,1), diagmat(tune)).t());
  mat K_star = eye(params.n_rows, params.n_rows);
  for(uword i=0; i < params.n_rows; i++) {
    for (uword j = (i+1); j< params.n_rows; j++){
      K_star(i, j) = exp(-dot(beta_star.t() % (params.row(i) - params.row(j)),(params.row(i) - params.row(j))));
      K_star(j, i) = K_star(i, j);
    }
  }
  // compute likelihoods
  vec L(T, fill::zeros), L_star(T, fill::zeros);
  for(uword i=0; i < T; i++){
    L_star(i) = dmvnorm(y.col(i).t(), theta.col(i), v*K_star, true)(0,0);
    L(i) = dmvnorm(y.col(i).t(), theta.col(i), v*K, true)(0,0);
  }
  // // compute MH ratio
  double logr = accu(L_star) + priorBeta(beta_star) + sum(log(1/beta))
    - accu(L) - priorBeta(beta) - sum(log(1/beta_star));
  //std::cout << logr << std::endl;
  if(log(R::runif(0,1)) < logr){
    beta = beta_star;
    K=K_star;
  }
  K = K + nugget*eye(K.n_rows, K.n_cols);
}


void calibrateParams(const mat& z, const cube& y, const cube& theta, const mat& params,
                     const mat& v, const vec tune, const double tune2, 
                     const vec& beta, const mat& K, const field<cube>& FF, rowvec& calibrate, 
                     double& sigma2_z, double a0, double b0, const int Sp){
  // compute S_t and M_t where M_t=phi[n,]*z[n-1]+r'*L*(y[,n]-to_matrix(F[,,n-lags]) * phi[n,])
  // S_t = sig_z+1-r'*L*r
  int T = y.n_cols;
  vec r(params.n_rows, fill::zeros), r_star(params.n_rows, fill::zeros);
  mat M(T, Sp, fill::zeros), M_star(T, Sp, fill::zeros), 
  S(T, Sp, fill::ones), S_star(T, Sp, fill::ones);
  mat Kinv = inv(K);
  rowvec calibrate_star = calibrate;
  //logit transform
  rowvec calibrate_star_logit = log(calibrate)-log(1-calibrate)+rmvnorm(1, zeros(calibrate.n_elem), diagmat(tune));
  calibrate_star = exp(calibrate_star_logit)/(1+exp(calibrate_star_logit));
  
  calibrate_star = calibrate + rmvnorm(1, zeros(calibrate.n_elem), diagmat(tune));
  
  for(int j=0; j < params.n_rows; j++){
    r(j) = exp(-dot(beta.t() % (calibrate - params.row(j)),(calibrate - params.row(j))));
    r_star(j) = exp(-dot(beta.t() % (calibrate_star - params.row(j)),(calibrate_star - params.row(j))));
    // exponential kernel rather than squared exponential
    // r(j) = exp(-sum(beta.t() % abs(calibrate - params.row(j))));
    // r_star(j) = exp(-sum(beta.t() % abs(calibrate_star - params.row(j))));
  }
  /// METROPOLIS STEP FOR CALIBRATION PARAMETERS
  // check if t is starting correctly, ie starting at 0 vs. 1
  int p = theta.slice(0).n_rows;
  vec regPars(p, fill::zeros);
  for(int sp=0; sp<Sp; sp++){
    for(int t=1; t<T; t++){
      regPars(0) = z(t-1,sp);
      for(int k=1; k < p; k++){
        regPars(k) = 1;
        regPars(k) = params.row(sp)(k-1);
      }
      M(t,sp) = as_scalar(dot(theta.slice(sp).col(t), regPars) + r.t() * Kinv * (y.slice(sp).col(t) - FF(sp).slice(t) * theta.slice(sp).col(t)));
      M_star(t,sp) = as_scalar(dot(theta.slice(sp).col(t), regPars) + r_star.t() * Kinv * (y.slice(sp).col(t) - FF(sp).slice(t) * theta.slice(sp).col(t)));
      S(t,sp) = as_scalar(sigma2_z + v(t,sp) - v(t) * r.t() * Kinv*r);
      S_star(t,sp) = as_scalar(sigma2_z + v(t, sp) - v(t,sp)*r_star.t()*Kinv*r_star);
    }
  }
  mat ll(1,1,fill::zeros), ll_star(1,1,fill::zeros);
  ll =  -0.5*sum(log(S)) - 0.5*sum(pow(z-M,2) / S) + calPrior(calibrate);
  ll_star = -0.5*sum(log(S_star)) - 0.5*sum(pow(z-M_star,2) / S_star) + calPrior(calibrate_star);
  
  // // compute MH ratio;
  mat logr = ll_star - ll 
    + log(prod(calibrate))+log(prod(1-calibrate))
    -log(prod(calibrate_star))-log(prod(1-calibrate_star));
    if(log(R::runif(0,1)) < logr(0,0)){
      calibrate = calibrate_star;
      M = M_star;
    }
    
    //// METROPOLIS STEP FOR OBSERVATIONAL VARIANCE
    double sigma2_star = exp(log(sigma2_z) + Rcpp::rnorm(1,0,tune2)(0));
    for(int j=0; j < params.n_rows; j++){
      r(j) = exp(-dot(beta.t() % (calibrate - params.row(j)),(calibrate - params.row(j))));
      /// exponential
      // r(j) = exp(-sum(beta.t() % abs(calibrate - params.row(j))));
    }
    for(int sp=0; sp<Sp; sp++){
      for(int t=1; t<T; t++){
        // regPars(0) = z(t-1,sp);
        // for(int k=1; k < p; k++){
        //   regPars(k) = 1;
        //   regPars(k) = params.row(sp)(k);
        // }
        // M(t, sp) = as_scalar(dot(theta.slice(sp).col(t), regPars) + r.t() * Kinv * (y.slice(sp).col(t) - FF(sp).slice(t) * theta.slice(sp).col(t)));
        S(t, sp) = as_scalar(sigma2_z + v(t) - v(t) * r.t() * Kinv * r);
        S_star(t, sp) = as_scalar(sigma2_star + v(t) - v(t) * r.t() * Kinv * r);
      }
    }
    ll =  -0.5*sum(log(S)) - 0.5*sum(pow(z-M,2) / S) + R::dgamma(sigma2_z,1,1,true);
    ll_star = -0.5*sum(log(S_star)) - 0.5*sum(pow(z-M,2) / S_star) +  R::dgamma(sigma2_star,1,1,true);
    
    logr = ll_star - ll + log(sigma2_star) - log(sigma2_z);
    if(log(R::runif(0,1)) < logr(0,0)){
      sigma2_z = sigma2_star;
    }
    
    // Gibbs step for time-varying observational variance
    // for(int t=0; t<T; t++){
    //   sigma2_z(t) = 1/Rcpp::rgamma(1.0,a0,
    //            1.0/(b0 + 0.5*dot(z(t) - M.row(t),
    //                              z(t) - M.row(t))))(0);
    //   // sigma2_z(t) = 0;
    // }
}

void sampleVar(mat z, mat emulated, double& sigma2){
  double SSR = accu(pow(z - emulated,2));
  double vinv = Rcpp::rgamma(1,z.n_elem/2.0,2.0/SSR)(0);
  sigma2 = 1.0 / vinv;
}


vec emulate(const mat& y, const rowvec& pred_params, const mat& params, 
            const vec& beta, const mat& theta,
            const mat& K, const cube& F, double yinit, const vec v){
  // this function will be built so that it gives one posterior predictive draw
  // thus, it must be called multiple times to build a posterior predictive distribution
  // the input arguments are one draw from the respective posterior parameters
  int T = y.n_cols;
  vec gamma_z(params.n_rows, fill::zeros);
  vec mu_z(T,fill::zeros);
  vec sigma_z(T,fill::ones);
  vec ypred(T, fill::ones);
  //vec v(T,fill::ones);
  
  mat Sigma_inv = inv(K); // covariance matrix
  mat errors(y.n_rows, T, fill::zeros);
  int p = F.slice(0).n_cols;
  for(int j = 0; j < params.n_rows; j++){
    gamma_z(j) = exp(-dot(beta.t() % (pred_params - params.row(j)),(pred_params - params.row(j))));
    // gamma_z(j) = exp(-sum(beta.t() % abs(pred_params - params.row(j))));
  }
  vec yinit2(1, fill::zeros);
  yinit2(0) = yinit;
  for(int t=0; t<T; t++){
    if(t==0){
      mu_z(t) = dot(theta.col(t),yinit2)+dot(gamma_z,Sigma_inv*(y.col(t) - F.slice(t)*theta.col(t)));
      sigma_z(t) = v(t)*(1-dot(gamma_z, Sigma_inv*gamma_z));
      ypred(t) = Rcpp::rnorm(1,mu_z(t), sqrt(sigma_z(t)))(0);
    }else{
      yinit2(0) = ypred(t-1);
      mu_z(t) = dot(theta.col(t),yinit2)+dot(gamma_z,Sigma_inv*(y.col(t) - F.slice(t)*theta.col(t)));
      sigma_z(t) = v(t)*(1-dot(gamma_z, Sigma_inv*gamma_z));
      ypred(t) = Rcpp::rnorm(1,mu_z(t), sqrt(sigma_z(t)))(0);
    }
  }
  return(ypred);
}



void forward2(const mat& y, const vec& m0, const mat& C0,
              const cube& FF, const mat& G,
              const mat& V, const mat& W,
              mat& a, mat& m, mat& f, mat& error, 
              cube& C, cube& Q, cube& R)
{
  int T = y.n_cols;
  a.col(0) = G * m0.col(0);
  R.slice(0) = G * C0 * G.t() + W;
  f.col(0) = FF.slice(0) * a.col(0);
  Q.slice(0) = FF.slice(0) * R.slice(0) * FF.slice(0).t() + V;
  error.col(0) = y.col(0) - f.col(0);
  m.col(0) = a.col(0) + R.slice(0) * FF.slice(0).t() * inv(Q.slice(0)) * error.col(0);
  C.slice(0) = R.slice(0) - R.slice(0) * FF.slice(0).t() * inv(Q.slice(0)) * FF.slice(0) * R.slice(0);
  for(int t=1; t<T; t++){
    // # predict
    a.col(t) = G * m.col(t-1);
    R.slice(t) = G * C.slice(t-1) * G.t() + W;
    // # // marginalize
    f.col(t) = FF.slice(t) * a.col(t);
    Q.slice(t) = FF.slice(t) * R.slice(t) * FF.slice(t).t() + V;
    // # // filter
    error.col(t) = y.col(t) - f.col(t);
    m.col(t) = a.col(t) + R.slice(t) * FF.slice(t).t() * inv(Q.slice(t)) * error.col(t);
    C.slice(t) = R.slice(t) - R.slice(t) * FF.slice(t).t() * inv(Q.slice(t)) * FF.slice(t) * R.slice(t);
  }
}

void backward2(mat& theta, mat& h, mat& H, const mat& a, const mat& m, const cube& C, 
               const cube& R, const mat& G)
{
  int T = a.n_cols;
  int p = a.n_rows;
  int p2 = G.n_cols;
  theta.col(T) = rmvnorm(1, m.col(T), C.slice(T)).t();
  //// expressions from Bayesian Filtering and Smoothing book
  for(int t=T-1; t>=0; t--){
    h = m.col(t) + C.slice(t)*G.t()*inv(R.slice(t))*(theta.col(t+1)-a.col(t));
    H = C.slice(t) - C.slice(t)*G.t()*inv(R.slice(t))*G*C.slice(t);
    theta.col(t) = rmvnorm(1, h, H).t();
  }
}

void emission(const mat& y, const mat& theta,
              const cube& FF, const mat&G, 
              const vec& v0, const mat& S0v, const mat& S0w, mat& V, mat& W)
{
  int T = y.n_cols, p = y.n_rows, p2 = G.n_cols, k, t;
  mat theta_ = theta;
  theta_.shed_col(0);
  mat SSy(p,p,fill::zeros), SSs(p2,p2,fill::zeros);
  for(t=0; t<T; t++){
    k = t+1;
    SSy += (y.col(t) - FF.slice(t)*theta_.col(t)) * (y.col(t) - FF.slice(t)*theta_.col(t)).t()/2;
    SSs += (G*theta.col(k) - G*theta.col(k-1)) * (G*theta.col(k) - G*theta.col(k-1)).t()/2;
  }
  V = inv(rwish(v0(0)+T/2.0, inv(S0v + SSy)));
  W = inv(rwish(v0(1)+T/2.0, inv(S0w + SSs)));
}


// [[Rcpp::export]]
Rcpp::List parallel_mvcalibrator_discrepancy(const mat& z, const cube& y, const uword& niter, 
                                             const int& burnin,
                                             mat& m0, mat& C0,
                                             const field<cube>& FF, const mat& G, 
                                             const mat& params,
                                             const double a0, const double b0, const double c0,
                                             vec& tune, const vec tune2, 
                                             const double tune3, const vec sptune,
                                             const double nugget, const int S,
                                             const vec& alpha0, const vec& beta0, 
                                             const double delta_v, const double delta_w,
                                             const double delta_vw, const double delta_ww,
                                             double a0z, double b0z,
                                             const vec yinit,
                                             const int discrep, double sigma_z)
{
  int T = y.n_cols;
  int p = y.n_rows;
  int p2 = G.n_cols;
  int p3 = 1;
  
  field<cube> FFH(S);
  cube ffh_blank(params.n_rows, p2, T, fill::zeros);
  FFH.fill(ffh_blank);
  
  cube mc_W(p2, p2,niter-burnin,fill::ones),
  mc_K(p, p, niter-burnin,fill::ones);
  
  mat mc_beta(params.n_cols, niter-burnin, fill::zeros);
  cube mc_v(T, S, niter-burnin, fill::zeros);
  
  field<cube> mc_theta(niter-burnin);
  cube blank(p2,p2,T,fill::zeros), blank2(p2, T, S,fill::zeros);
  // W.fill(blank);
  mc_theta.fill(blank2);
  cube theta(p2,T,S,fill::zeros), mc_w(T, S, niter-burnin, fill::zeros);
  // for(int i=0;i<S;i++){
  //   for(int j=0; j<T; j++){
  //     W(i).slice(j) =  eye(p2,p2);
  //   }
  // }
  
  mat K(p,p,fill::eye);
  mat v(T, S, fill::ones);
  vec gpbeta = pow(Rcpp::rnorm(params.n_cols, 5, 1),2), gpbeta2 = pow(Rcpp::rnorm(2, 5, 1), 2);
  for(uword i=0; i < p; i++) {
    for (uword j = (i+1); j < p; j++){
      K(i, j) = exp(-dot(gpbeta.t() % (params.row(i) - params.row(j)),(params.row(i) - params.row(j))));
      K(j, i) = K(i, j);
    }
  }
  double progress = 0;
  std::vector<std::tuple<mat,vec>> forwardBackwardRuns;
  for(int sp=0; sp<S; sp++){
    forwardBackwardRuns.push_back(
      std::make_tuple(ones(p2,T), ones(T))
    );
  }
  for(uword i=0; i<niter; i++) {
    // loop over each spatial location and run emulator
    progress = ((double) i) / ((double) niter);
    if(progress == .1){
      std::cout << "10% sampled" << std::endl;
    }else if(progress == .25){
      std::cout << "25% sampled" << std::endl;
    }else if(progress == .5){
      std::cout << "50% sampled" << std::endl;
    }else if(progress == .75){
      std::cout << "75% sampled" << std::endl;
    }else if(progress == .95){
      std::cout << "emulator fit!" << std::endl;
    }
    RcppThread::parallelFor(0, S, [&] (int sp){
      // forward filter, backward sample
      if(p2>1){
        for(int i=0;i<T;i++){
          // this can either include dynamic regression components through params
          // or just dynamic intercept through a vector of ones
          // FFH(sp).slice(i) = join_rows(FF(sp).slice(i), ones(params.n_rows,1));
          FFH(sp).slice(i) = join_rows(FF(sp).slice(i), params);
        }
      }else{
        FFH(sp) = FF(sp);
      }
      forwardBackwardRuns[sp] = forwardBackward(y.slice(sp),
                                                m0, C0, a0, b0, c0,
                                                FF(sp), eye(p2,p2),
                                                K, alpha0, beta0,
                                                delta_w, delta_v);
      theta.slice(sp) = std::get<0>(forwardBackwardRuns[sp]);
      v.col(sp) = std::get<1>(forwardBackwardRuns[sp]);
    });
    kernel(y, params, theta, gpbeta, FFH, K, tune, v, S, nugget);
    
    // // // save mcmc sweeps
    if(i >= burnin){
      mc_v.slice(i-burnin) = v;
      mc_beta.col(i-burnin) = gpbeta;
      mc_theta(i-burnin) = theta;
      mc_K.slice(i-burnin) = K;
    }
  }
  
  
  // calibration step
  int draw, draw2;
  rowvec calibrate = Rcpp::runif(params.n_cols,.5,.5);
  vec yeta_draw(T, fill::zeros);
  mat yeta_draw_mat(T, S, fill::zeros);
  
  // calibration step
  
  mat aw(S, T, fill::zeros), mw(p3, T, fill::zeros), 
  fw(S, T, fill::zeros), 
  errorw(S, T, fill::zeros);
  cube theta_w(p3, T, S, fill::zeros);
  mat vw(T, S, fill::ones);
  
  cube Rw(S, S, T, fill::ones), Qw(S, S, T, fill::ones), 
  Cw(p3, p3, T, fill::zeros);
  
  cube Ww(p3, p3, T, fill::ones), Vw(S, S, T, fill::ones);
  mat Vw_(S, S, fill::eye), Ww_(p3, p3, fill::eye);
  
  vec y_eta(T, fill::zeros);
  
  vec theta0(S, fill::zeros);
  cube FFw(1, p3, T, fill::zeros);
  
  vec m0w(p3, fill::zeros);
  mat C0w(p3, p3, fill::eye);
  mat w(T, S, fill::zeros);
  
  mat mc_calibrate(niter-burnin, params.n_cols, fill::zeros);
  mat mc_sigma(1, niter-burnin, fill::ones);
  cube mc_Vw(T, S,niter-burnin,fill::zeros),
  mc_Ww(p3, p3, niter-burnin,fill::zeros);
  mat mc_beta2(2, niter-burnin, fill::zeros);
  cube mc_yeta(T, niter-burnin, S, fill::zeros);
  
  // mat mc_calibrate(numSamples*numBlocks, params.n_cols, fill::zeros);
  // mat mc_sigma(1, numSamples*numBlocks, fill::ones);
  // cube mc_w(S, T+1, numSamples*numBlocks, fill::zeros);
  // cube mc_Vw(S,S,numSamples*numBlocks,fill::zeros),
  // mc_Ww(S,S,numSamples*numBlocks,fill::zeros);
  // mat mc_beta2(2, numSamples*numBlocks, fill::zeros);
  // cube mc_yeta(T, S, numSamples*numBlocks, fill::zeros);
  
  mat temp;
  mat H(p3, p3,fill::eye), h(p3, 1, fill::ones);
  // try non-blocked version because beta is fucked up along with discrepancy
  for(int j=0; j<niter; j++){
    progress = ((double) j) / ((double) niter);
    if(progress == .1){
      std::cout << "10%" << std::endl;
    }else if(progress == .25){
      std::cout << "25%" << std::endl;
    }else if(progress == .5){
      std::cout << "50%" << std::endl;
    }else if(progress == .75){
      std::cout << "75%" << std::endl;
    }else if(progress == .95){
      std::cout << "calibrator fit!" << std::endl;
    }
    draw = (int) Rcpp::runif(1, 0, niter-burnin)(0);
    calibrateParams(z, y,
                    mc_theta(draw),
                    params,
                    mc_v.slice(draw),
                    tune2, tune3,
                    mc_beta.col(draw),
                    mc_K.slice(draw),
                    FF,
                    calibrate, sigma_z, a0z, b0z, S);
    // after drawing calibration parameters, need to get an emulation run at those
    // values, so now emulate at each spatial location
    if(j>=burnin){
      for(int sp=0; sp<S; sp++){
        yeta_draw = emulate(y.slice(sp), calibrate,
                            params, mc_beta.col(draw),
                            mc_theta(draw).slice(sp),
                            mc_K.slice(draw),
                            FF(sp),
                            yinit(sp),
                            mc_v.slice(draw).col(sp));
        mc_yeta.slice(sp).col(j-burnin) = yeta_draw;
      }
      for(int sp=0; sp<S; sp++){
        yeta_draw_mat.col(sp) = mc_yeta.slice(sp).col(j-burnin);
      }
      sampleVar(z, yeta_draw_mat, sigma_z);
      
      mc_calibrate.row(j-burnin) = calibrate;
      mc_sigma.col(j-burnin) = sigma_z;
    }
  }
  
  // where y_draw is located in the following, fix it through looping and building appropriate
  // matrix containing 
  rowvec mc_cal_colMeans(params.n_cols, fill::zeros);
  for(int i=0; i<mc_calibrate.n_rows;i++){
    mc_cal_colMeans += mc_calibrate.row(i) / mc_calibrate.n_rows;
  }
  for(int j=0; j<niter; j++){
    if(discrep!=0){
      progress = ((double) j) / ((double) niter);
      if(progress == .1){
        std::cout << "10%" << std::endl;
      }else if(progress == .25){
        std::cout << "25%" << std::endl;
      }else if(progress == .5){
        std::cout << "50%" << std::endl;
      }else if(progress == .75){
        std::cout << "75%" << std::endl;
      }else if(progress == .95){
        std::cout << "discrepancy fit!" << std::endl;
      }
      draw = (int) Rcpp::runif(1, 0, niter-burnin)(0);
      for(int sp=0; sp<S; sp++){
        yeta_draw = emulate(y.slice(sp), mc_calibrate.row(draw), params, mc_beta.col(draw),
                            mc_theta(draw).slice(sp),
                            mc_K.slice(draw),
                            FF(sp),
                            yinit(sp),
                            mc_v.slice(draw).col(sp));
        mc_yeta.slice(sp).col(j) = yeta_draw;
      }
      Vw_ = sigma_z*eye(S,S);
      // consider d-inverse gamma model by looping over spatial locations
      
      // modify FFw to be dynamic regression on the calibration parameters
      // forward2(z.t()-yeta_draw_mat.t(), m0w, C0w, FFw, eye(S,S), Vw_, Ww_, aw, mw, fw,
      //          errorw, Cw, Qw, Rw);
      // //
      // // // // // // // backward sample
      // backward2(theta_w, h, H, aw, join_rows(m0w,mw), join_slices(C0w,Cw), Rw, eye(S,S));
      // //
      // // // // // //
      // emission(z.t()-yeta_draw_mat.t(), theta_w, FFw, eye(S,S), ones(S), diagmat(ones(S)), diagmat(ones(S)), Vw_, Ww_);
      //
      // // https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1541-0420.2008.01115.x
      // try a d-inverse gamma model form for W_ instead of a GP for debugging
      // kernel2(z.t()-yeta_draw_mat.t(), spatialLocs, theta_w, gpbeta2, Ww_, sptune,
      //         sigma_z, S, nugget);
      //
      for(int i=0;i<T;i++){
        // FFw.slice(i) = join_horiz(ones(1), mc_calibrate.row(draw));
        FFw.slice(i) = eye(1,1);
      }
      RcppThread::parallelFor(0, S, [&] (int sp){
        // forward filter, backward sample
        // change this to non-discount version?
        // or estimate V correlation through a GP? or estimate W correlation through a GP
        forwardBackwardRuns[sp] = forwardBackward(z.col(sp).t()-yeta_draw_mat.col(sp).t(),
                                                  m0w, C0w, a0, b0, c0,
                                                  FFw, eye(p3, p3),
                                                  eye(1,1), alpha0, beta0,
                                                  delta_vw, delta_ww);
        theta_w.slice(sp) = std::get<0>(forwardBackwardRuns[sp]);
        vw.col(sp) = std::get<1>(forwardBackwardRuns[sp]);
        if(j>=burnin){
          mc_w.slice(j-burnin).col(sp) = theta_w.slice(sp).t();// theta_w.slice(sp).t() * join_cols(ones(1),mc_calibrate.row(draw).t());
        }
      });
    }
    if(j>=burnin){
      mc_Vw.slice(j-burnin) = vw;
      mc_Ww.slice(j-burnin) = Ww_;
    }
  }
  
  
  Rcpp::List out(5);
  out["v"] = mc_v;
  out["mc_w"] = mc_w;
  out["mc_Ww"] = mc_Ww;
  out["mc_Vw"] = mc_Vw;
  out["mc_yeta"] = mc_yeta;
  out["yeta_draw"] = yeta_draw;
  out["mc_sigma2"] = mc_sigma;
  out["mc_theta"] = mc_theta;
  out["thetaw"] = theta_w;
  out["mc_w"] = mc_w;
  out["mc_beta"] = mc_beta;
  out["mc_beta2"] = mc_beta2;
  out["mc_K"] = mc_K;
  out["calibrate"] = mc_calibrate;
  return(out);
}