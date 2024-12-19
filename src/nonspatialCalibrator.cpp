#include <RcppDist.h>
#include <chrono>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace arma;
using namespace std::chrono;

void kalman(const mat& y, const vec& m0, const mat& C0,
            const cube& FF,
            const mat& Sigma, const mat& W,
            mat& a, cube& R, mat& m,
            cube& C, vec& alpha, vec& beta,
            const vec& alpha0, const vec& beta0, 
            const double delta_w, const double delta_v, const int S){
  uword T = y.n_cols;
  a.col(0) =  m0.col(0);
  R.slice(0) = C0/delta_w;
  // R.slice(0) = C0 + W;
  mat f = FF.slice(0) * a.col(0);
  mat V = kron(eye(S,S), Sigma);
  mat invV = kron(eye(S,S), inv(Sigma));
  mat Q = FF.slice(0) * R.slice(0) * FF.slice(0).t() + V;
  // mat Qinv = inv(Q);
  mat Qinv = invV - invV*FF.slice(0)*
    inv(inv(R.slice(0))+FF.slice(0).t()*invV*FF.slice(0))*
    FF.slice(0).t()*invV;
  m.col(0) = a.col(0) + R.slice(0) * FF.slice(0).t() * Qinv * (y.col(0) - f.col(0));
  C.slice(0) = R.slice(0) - R.slice(0) * FF.slice(0).t() * Qinv * FF.slice(0) * R.slice(0);
  // modify alpha and beta for discount factor
  alpha(0) = as_scalar(delta_v*alpha0) + 0.5*y.n_rows;
  beta(0) = as_scalar(delta_v*beta0) + 0.5*dot(y.col(0) - f.col(0), Qinv * (y.col(0) - f.col(0)));
  for(uword t=1; t<T; t++){
    // # // predict
    a.col(t) =  m.col(t-1);
    R.slice(t) =  C.slice(t-1)/delta_w;
    // R.slice(t) =  C.slice(t-1) + W;
    // # // marginalize
    f = FF.slice(t) * a.col(t);
    Q = FF.slice(t) * R.slice(t) * FF.slice(t).t() + V;
    // Qinv = inv(Q);
    Qinv = invV - invV*FF.slice(t)*
      inv(inv(R.slice(t))+FF.slice(t).t()*invV*FF.slice(t))*
      FF.slice(t).t()*invV;
    // # // discount factors
    alpha(t) = delta_v*alpha(t-1) + 0.5*y.n_rows;
    beta(t) = delta_v*beta(t-1) + 0.5*dot(y.col(t) - f, Qinv * (y.col(t) - f));
    // # // filter
    m.col(t) = a.col(t) + R.slice(t) * FF.slice(t).t() * Qinv * (y.col(t) - f);
    C.slice(t) = (R.slice(t) - R.slice(t) * FF.slice(t).t() * Qinv * FF.slice(t) * R.slice(t));
  }
}

void backwardVolatility(const vec& n, const vec& d, 
                        const double delta_v, vec& v){
  uword T = n.n_elem;
  double vinv, v_star;
  vinv = Rcpp::rgamma(1,n(T-1), 1.0/d(T-1))(0);
  v(T-1) = 1.0/vinv;
  for(int t=T-2; t>=0; t--){
    v_star = Rcpp::rgamma(1,(1-delta_v)*n(t+1), 1.0/d(t))(0);
    vinv = delta_v*vinv + v_star;
    v(t) = 1.0/vinv;
  }
}

void backwardState(mat& theta, 
                   const mat& a, 
                   const mat& m, 
                   const cube& CC, 
                   const cube& RR,
                   const vec v,
                   const double nugget){
  uword T = a.n_cols;
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
}

// eventually change this to suggestion by Higdon
// check Stan code
double priorBeta(vec beta){
  vec out(beta.n_elem, fill::zeros);
  for(uword i=0; i < out.n_elem; i++){
    out(i) = R::dgamma(beta(i),.1,.1,true);
  }
  return(sum(out));
}

// void kernel(const mat& y, const mat& params, const mat& theta, vec& beta,
//             const cube& FF, mat& K, const vec tune, const vec& v, 
//             const double nugget, const int S){
//   uword T = y.n_cols;
//   uword p = params.n_rows;
//   //compute proposed kernel
//   vec eta = log(beta);
//   vec beta_star = exp(eta+rmvnorm(1, zeros(beta.n_rows,1), diagmat(tune)).t());
//   mat K_star = eye(params.n_rows, params.n_rows);
//   for(uword i=0; i<params.n_rows; i++){
//     for (uword j = (i+1); j< params.n_rows; j++){
//       K_star(i, j) = exp(-dot(beta_star.t() % (params.row(i) - params.row(j)),(params.row(i) - params.row(j))));
//       K_star(j, i) = K_star(i, j);
//     }
//   }
//   mat Sigma_star = kron(eye(S, S), K_star);
//   mat Sigma = kron(eye(S, S), K);
//   
//   // compute likelihoods
//   vec L(T, fill::zeros), L_star(T, fill::zeros);
//   for(uword i=0; i < L.n_elem; i++){
//     L_star(i) = dmvnorm(y.col(i).t(), FF.slice(i)*theta.col(i), v(i)*(Sigma_star), true)(0,0);
//     L(i) = dmvnorm(y.col(i).t(), FF.slice(i)*theta.col(i), v(i)*Sigma, true)(0,0);
//   }
//   // // // compute MH ratio
//   double logr = sum(L_star) + priorBeta(beta_star) + sum(log(1/beta))
//     - sum(L) - priorBeta(beta) - sum(log(1/beta_star));
//   //std::cout << logr << std::endl;
//   if(log(R::runif(0,1)) < logr){
//     beta = beta_star;
//     K=K_star;
//   }
//   K = K + nugget * eye(K.n_cols, K.n_cols);
// }

void kernel(const cube& y, const mat& params, const mat& theta, vec& beta,
            const field<cube>& FF, mat& K, const vec tune, const vec& v, 
            const double nugget, const int S)
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
  mat Sigma_star = kron(eye(S, S), K_star);
  mat Sigma = kron(eye(S, S), K);
  // compute likelihoods
  mat L(S, T, fill::zeros), L_star(S, T, fill::zeros);
  for(uword i=0; i < T; i++){
    for(uword sp=0; sp < S; sp++){
      L_star(sp,i) = dmvnorm(y.slice(sp).col(i).t(), FF(sp).slice(i)*theta(sp, i), v(i)*K_star, true)(0,0);
      L(sp, i) = dmvnorm(y.slice(sp).col(i).t(), FF(sp).slice(i)*theta(sp, i), v(i)*K, true)(0,0);
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

double calPrior(rowvec calibrate){
  rowvec out(calibrate.n_elem, fill::zeros);
  for(uword i=0; i < out.n_elem; i++){
    out(i) = R::dbeta(calibrate(i),3,3,true);
  }
  return(sum(out));
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


void calibrateParams(const mat& z, const mat& y,const cube& yy,
                     const mat& theta, const mat& params,
                     const vec& v, const vec tune, const double tune2, 
                     const vec& beta, const mat& K,
                     const cube FF, const field<cube>& FFc, rowvec& calibrate, 
                     double& sigma2_z, double a0, double b0, const int Sp)
{
  // compute S_t and M_t where M_t=phi[n,]*z[n-1]+r'*L*(y[,n]-to_matrix(F[,,n-lags]) * phi[n,])
  // S_t = sig_z+1-r'*L*r
  int T = y.n_cols;
  vec r(params.n_rows, fill::zeros), r_star(params.n_rows, fill::zeros);
  mat M(Sp, T, fill::zeros), M_star(Sp, T, fill::zeros), 
  S(Sp, T, fill::ones), S_star(Sp, T, fill::ones);
  mat Kinv = inv(K);
  rowvec calibrate_star = calibrate;
  //logit transform
  rowvec calibrate_star_logit = log(calibrate)-log(1-calibrate)+rmvnorm(1, zeros(calibrate.n_elem), diagmat(tune));
  calibrate_star = exp(calibrate_star_logit)/(1+exp(calibrate_star_logit));
  
  // calibrate_star = calibrate + rmvnorm(1, zeros(calibrate.n_elem), diagmat(tune));
  
  for(int j=0; j < params.n_rows; j++){
    r(j) = exp(-dot(beta.t() % (calibrate - params.row(j)),(calibrate - params.row(j))));
    r_star(j) = exp(-dot(beta.t() % (calibrate_star - params.row(j)),(calibrate_star - params.row(j))));
  }
  // vec r = join_vert(rr,rr);
  // vec r_star = join_vert(rr_star,rr_star);
  /// METROPOLIS STEP FOR CALIBRATION PARAMETERS
  // for(int t=1; t<T; t++){
  //   for(int s=0; s<Sp; s++){
  //     M(s, t) = theta(s, t)*z(s, t-1) + dot(r, Kinv * (y.col(t) - FF.slice(t) * theta.col(t)));
  //     M_star(s, t) = theta(s, t)*z(s, t-1) + dot(r_star, Kinv * (y.col(t) - FF.slice(t) * theta.col(t)));
  //     S(s, t) = sigma2_z + v(t) - v(t)*dot(r, Kinv * r);
  //     S_star(s, t) = sigma2_z + v(t) - v(t)*dot(r_star, Kinv * r_star);
  //   }
  // }
  mat ll(1,1,fill::zeros), ll_star(1,1,fill::zeros);
  double sigma_t, mu_t;
  double ll_counter=0, ll_counter_star=0;
  
  for(int s=0; s<Sp; s++){
    for(int t=1; t<T; t++){
      sigma_t = sigma2_z + v(t) - v(t)*dot(r, Kinv * r);
      mu_t = theta(s,t)*z(s,t-1) + dot(r, Kinv * (yy.slice(s).col(t) - FFc(s).slice(t)*theta(s, t)));
      ll_counter += -0.5*log(sigma_t) - 0.5*pow((z(s,t) - mu_t)/sqrt(sigma_t),2);
      
      sigma_t = sigma2_z + v(t) - v(t)*dot(r_star, Kinv * r_star);
      mu_t = theta(s,t)*z(s,t-1) + dot(r_star, Kinv * (yy.slice(s).col(t) - FFc(s).slice(t)*theta(s, t)));
      ll_counter_star += -0.5*log(sigma_t) - 0.5*pow((z(s,t) - mu_t)/sqrt(sigma_t),2);
    }
  }
  ll_star = ll_counter_star, ll=ll_counter;
  
  // // compute MH ratio;
  mat logr = ll_star - ll
    + log(prod(calibrate))+log(prod(1-calibrate))
    -log(prod(calibrate_star))-log(prod(1-calibrate_star));
    if(log(R::runif(0,1)) < logr(0,0)){
      calibrate = calibrate_star;
      M = M_star;
    }
}

void sampleVar(mat z, mat emulated, double& sigma2){
  double SSR = accu(pow(z - emulated,2));
  double vinv = Rcpp::rgamma(1,z.n_elem/2.0,2.0/SSR)(0);
  sigma2 = 1.0 / vinv;
}


// [[Rcpp::export]]
Rcpp::List DLMGP(const mat& z, const mat& y,  const cube& yy, const vec yinit,
                              const uword& niter, const uword& burnin,
                              mat& m0, mat& C0, mat& W,
                              const cube& FF, const field<cube>& FFc,
                              const mat& params,
                              const rowvec pred_params,
                              const double delta_v, const double delta_w,
                              const vec& alpha0, const vec& beta0, 
                              vec& tune, const vec& tune2, double tune3, const vec& tune4,
                              double& sigma_z, const double nugget,
                              const double a0z, const double b0z, 
                              const int S, const int S2)
{
  int T = y.n_cols;
  int p = y.n_rows;
  int p2 = S;
  // initialize and allocate variables necessary for forward filtering and 
  // backward sampling
  // for code readability, order according to Kalman filter, then backward sampele, etc
  
  field<cube> mc_W(niter-burnin);
  for(int i=0; i<niter-burnin; i++){
    mc_W(i) = cube(p2,p2,T,fill::zeros);
  }
  cube mc_K(params.n_rows, params.n_rows, niter-burnin, fill::ones), 
  mc_theta(p2, T, niter-burnin,fill::zeros);
  mat V(p,p,fill::eye), 
  K(params.n_rows,params.n_rows,fill::eye);
  cube K_cube(params.n_rows, params.n_rows, S, fill::zeros),
  C(p2, p2, T, fill::zeros), R(p2, p2, T, fill::ones);
  ;
  for(uword i=0; i < S; i++){
    K_cube.slice(i) = K;
  }
  for(uword i=0; i < T; i++){
    R.slice(i) = eye(p2,p2);
    C.slice(i) = eye(p2,p2);
  }
  
  mat a(p2,T,fill::zeros), m(p2,T,fill::zeros),
  Q(p,p,fill::zeros), theta(p2,T,fill::zeros), h(p2, 1, fill::ones);
  
  cube mc_C(p2, p2, niter, fill::zeros);
  mat spparams = params;
  vec gpbeta(params.n_cols, fill::ones), spgpbeta(spparams.n_cols, fill::ones),
  alpha(T, fill::ones), beta(T, fill::ones), s(T, fill::ones), 
  f(p, fill::zeros), v(T, fill::zeros);
  
  gpbeta = pow(Rcpp::rnorm(params.n_cols, 5, 1), 2);
  spgpbeta = pow(Rcpp::rnorm(spparams.n_cols, 5, 1) ,2);
  
  mat mc_beta(params.n_cols, niter-burnin, fill::zeros);
  mat mc_v(T, niter-burnin, fill::zeros);
  double progress = 0;
  auto start=high_resolution_clock::now(), start_overall=high_resolution_clock::now();
  auto stop=high_resolution_clock::now(), stop_overall=high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  double kalman_time=0, backwardState_time=0, kernel_time=0;
  start_overall = high_resolution_clock::now();
  for(uword i=0; i<niter; i++) {
    // loop over each spatial location and run emulator
    progress = ((double) i) / ((double) niter);
    if(progress == .1){
      std::cout << "10%" << std::endl;
    }else if(progress == .25){
      std::cout << "25%" << std::endl;
    }else if(progress == .5){
      std::cout << "50%" << std::endl;
    }else if(progress == .75){
      std::cout << "75%" << std::endl;
    }else if(progress == .95){
      std::cout << "95%" << std::endl;
    }
    start = high_resolution_clock::now();
    kalman(y, m0, C0, FF, K, W, a, R, m, C,
           alpha, beta, alpha0, beta0, delta_w, delta_v, S2);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    kalman_time += duration.count() / ((double) niter);
    
    // // // backward sample volatility
    backwardVolatility(alpha, beta, delta_v, v);
    //
    // // // backward sample latent state
    start = high_resolution_clock::now();
    backwardState(theta, a, m, C, R, v, nugget);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    backwardState_time += duration.count() / ((double) niter);
    //
    // // // // // // sample kernel parameters
    start = high_resolution_clock::now();
    kernel(yy, params, theta, gpbeta, FFc, K, tune, v, nugget, S);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    kernel_time += duration.count() / ((double) niter);

    // // save mcmc sweeps
    if(i >= burnin){
      mc_v.col(i-burnin) = v;
      mc_beta.col(i-burnin) = gpbeta;
      mc_theta.slice(i-burnin) = theta;
      mc_K.slice(i-burnin) = K;
      mc_C.slice(i-burnin) = C.slice(3);
    }
  }
  std::cout << "AVG KALMAN TIME:" << kalman_time << std::endl;
  std::cout << "AVG BACKWARD TIME:" << backwardState_time << std::endl;
  std::cout << "AVG KERNEL TIME:" << kernel_time << std::endl;

  // calibration step
  mat mc_calibrate(niter-burnin, params.n_cols, fill::zeros);
  mat mc_sigma(1, niter-burnin, fill::ones);
  int draw;
  rowvec calibrate = Rcpp::runif(params.n_cols,.5,.5);
  
  cube yeta_draw(T,niter,S2,fill::zeros);
  mat yeta_draw_mat(T, S2, fill::zeros);
  for(int i=0; i<niter; i++) {
    draw = (int) Rcpp::runif(1, 0, niter-burnin)(0);
    for(int s=0;s<S2;s++){
      yeta_draw.slice(s).col(i) = emulate(yy.slice(s), pred_params,
                      params, mc_beta.col(draw),
                      mc_theta.slice(draw).row(0),
                      mc_K.slice(draw),
                      FFc(s),
                      yinit(s),
                      mc_v.col(draw));
    }
    calibrateParams(z, y, yy,
                    mc_theta.slice(draw),
                    params,
                    mc_v.col(draw),
                    tune2, tune3,
                    mc_beta.col(draw),
                    mc_K.slice(draw),
                    FF, FFc,
                    calibrate, sigma_z, a0z, b0z, S);
    
    for(int sp=0; sp<S2; sp++){
      yeta_draw_mat.col(sp) = yeta_draw.slice(sp).col(i);
    }
    sampleVar(z, yeta_draw_mat, sigma_z);
    
    if(i>=burnin){
      mc_sigma.col(i-burnin) = sigma_z;
      mc_calibrate.row(i-burnin) = calibrate;
    }
  }
  stop = high_resolution_clock::now();
  duration = duration_cast<seconds>(stop - start_overall);
  std::cout << "TOTAL TIME:" << duration.count() << std::endl;
  Rcpp::List out(5);
  out["yeta"] = yeta_draw;
  out["v"] = mc_v;
  out["mc_sigma2"] = mc_sigma;
  out["mc_theta"] = mc_theta;
  out["mc_beta"] = mc_beta;
  out["mc_K"] = mc_K;
  out["calibrate"] = mc_calibrate;
  out["mc_C"] = mc_C;
  return(out);
}



