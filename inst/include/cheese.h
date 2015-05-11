#ifndef __CHEESE
#define __CHEESE


#define LOG_PI 1.1447298858494001639 // log(pi)

#include <mb_base.h>
#include <except.h>
#include <utilfuncs.h>
#include <Eigen/Core>
#include <cppad_atomics.h>
#include <MVN_AD.h>
#include <wish_AD.h>


using Rcpp::List;
using Rcpp::as;
using Eigen::Matrix;
using Eigen::Array;
using Eigen::MatrixBase;
using Eigen::Dynamic;
using Eigen::Block;
using Eigen::Map;
using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::VectorXd;


template<typename T>
bool my_finite(const T& x) {
  return( (abs(x) <= __DBL_MAX__ ) && ( x == x ) );
}


class cheese  {

  typedef Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
  typedef Matrix<AScalar, Dynamic, 1> VectorXA;
  typedef Array<AScalar, Dynamic, Dynamic> ArrayXA;

 public:
  
  cheese(const List&);
  
  template<typename Tpars>
    AScalar eval_f(const MatrixBase<Tpars>&);
  
  AScalar eval_LL(); 
  
  template<typename Tpars>
    typename Tpars::Scalar eval_LL(const MatrixBase<Tpars>&);
  
  AScalar eval_LLi(); 


  template<typename TL>
    void eval_LLi(const MatrixBase<TL>&);


  template<typename TL>
    void eval_prior_i(const MatrixBase<TL>&);


  AScalar eval_prior();
  AScalar eval_hyperprior();
  
  template<typename Tpars>
    typename Tpars::Scalar eval_prior(const MatrixBase<Tpars>&);
  
 
  template<typename Tpars>
    typename Tpars::Scalar eval_hyperprior(const MatrixBase<Tpars>&);
  
  template<typename Tpars>
    void unwrap_params(const MatrixBase<Tpars>&);

  int get_nvars();
  int nvars; 
  
  // data and priors.

  MatrixXA volume;
  MatrixXA log_price;
  MatrixXA disp;
  MatrixXA log_volume;

  MatrixXA mu_prior_mean;
  MatrixXA mu_prior_chol_prec;
  AScalar nu;
  MatrixXA chol_A;
  AScalar alpha_prior;

  int N, k, maxT; // number of observations per unit, number of units, number of regressors
  
  VectorXi T;
  // declare parameters

  MatrixXA B;
  VectorXA mu;
  MatrixXA chol_G; // Cholesky of precision of B
  VectorXA log_alpha;
  VectorXA alpha;

};



cheese::cheese(const List& params)
{

  List & pars = static_cast<List &>(const_cast<List &>(params));

  List data = as<List>(pars["data"]);
  List priors = as<List>(pars["priors"]);

  // Map R to Rcpp
  NumericMatrix volume_((SEXP)data["volume"]);
  NumericMatrix log_price_((SEXP)data["log_price"]);
  NumericMatrix disp_((SEXP)data["disp"]);
  NumericVector mu_prior_mean_((SEXP)priors["mu_prior_mean"]);
  NumericMatrix mu_prior_chol_prec_((SEXP)priors["mu_prior_chol_prec"]);
  nu = as<double>((SEXP)priors["nu"]);
  NumericMatrix chol_A_((SEXP)priors["chol_A"]);
  alpha_prior = as<double>((SEXP)priors["alpha_prior"]);
  
  Rcpp::IntegerVector T_((SEXP)data["T"]);

  N = volume_.ncol();
  k = mu_prior_mean_.size();
  maxT = volume_.nrow();

  // Map Rcpp to Eigen Scalar
  Map<MatrixXd> log_price_d = MatrixXd::Map(log_price_.begin(),maxT,N);
  Map<MatrixXd> disp_d = MatrixXd::Map(disp_.begin(),maxT,N);
  Map<MatrixXd> volume_d = MatrixXd::Map(volume_.begin(),maxT,N);

  Map<MatrixXd> mu_prior_meand = MatrixXd::Map(mu_prior_mean_.begin(),k,1);
  Map<MatrixXd> mu_prior_chol_precd = MatrixXd::Map(mu_prior_chol_prec_.begin(),k,k);
  Map<MatrixXd> chol_Ad = MatrixXd::Map(chol_A_.begin(),k,k);
  T = VectorXi::Map(T_.begin(),N);

  // copy Eigen Scalar Map to Eigen AScalar objects
  volume = volume_d.cast<AScalar>();
  log_price = log_price_d.cast<AScalar>();
  disp = disp_d.cast<AScalar>();
  
  log_volume.resize(maxT,N);
  log_volume.array() = volume.array().log();

  mu_prior_mean = mu_prior_meand.cast<AScalar>();
  mu_prior_chol_prec = mu_prior_chol_precd.cast<AScalar>();
  chol_A = chol_Ad.cast<AScalar>();


  // reserve space for parameters
  chol_G.resize(k,k);
  B.resize(k,N);
  mu.resize(k);
  log_alpha.resize(N);

  alpha.resize(N);

}




template<typename Tpars>
void cheese::unwrap_params(const MatrixBase<Tpars>& P)
{

  // here we are copying the parameters.  Could we map instead? Doubtful because of the cast. 

  B = MatrixXA::Map(P.derived().data(), k, N);
  log_alpha = P.derived().segment(k*N,N);
  mu = P.derived().segment(k*N+N,k);

  chol_G.setZero();   
  int idx = N*k+N+k;
  for (int j=0; j<k; j++) {
    for (int i=j; i<k; i++) {
      chol_G(i,j) = P(idx); // fill lower triangle only
      idx++;
    }
    chol_G(j,j) = exp(chol_G(j,j));
  }

  alpha.array() = log_alpha.array().exp();

}




template<typename Tpars>
AScalar cheese::eval_f(const MatrixBase<Tpars>& P) {
  


  unwrap_params(P);

  AScalar LL = eval_LL();
  AScalar prior = eval_prior();
  AScalar hyperprior = eval_hyperprior();
  
  AScalar f = LL + prior + hyperprior;
  
  
  return(f);
}




template<typename TL>
void cheese::eval_LLi(const MatrixBase<TL>& LLi_)
{ 

  MatrixBase<TL>& LLi = const_cast<MatrixBase<TL>& >(LLi_);

  for (int i=0; i<N; i++) {
    int ti = T[i];
    VectorXA mi = VectorXA::Constant(ti,B(0,i));
 
    mi += B(1,i)*log_price.col(i).head(ti) + B(2,i)*disp.col(i).head(ti);
    VectorXA ri(ti);
    ri.array() = alpha(i)*mi.array().exp();
    VectorXA lgamma_ri(ti);
    for (int j=0; j<ti; j++) {
      lgamma_ri(j) = lgamma(ri(j));
    }
    VectorXA LLit(ti);
    LLit.array() = ri.array()*log_alpha(i)-lgamma_ri.array();
    LLit.array() +=  (ri.array()-AScalar(1.))*log_volume.col(i).head(ti).array() - alpha(i)*volume.col(i).head(ti).array();
    LLi(i) = LLit.sum();
     assert(my_finite(LLi(i)));
  }

}




AScalar cheese::eval_LL()
{ 

  VectorXA LLi(N); //vector to store individual-level likelihoods
  eval_LLi(LLi);
  AScalar LL = LLi.sum();
  return(LL);

}


template<typename TL>
void cheese::eval_prior_i(const MatrixBase<TL>& prior_i_) {
  
  MatrixBase<TL>& prior_i = const_cast<MatrixBase<TL>& >(prior_i_);

  // MVN prior on B

  prior_i.setConstant(chol_G.diagonal().array().log().sum() - k*M_LN_SQRT_2PI);

  ArrayXA Z = ((B.colwise() - mu).transpose() * chol_G).array();
  prior_i.array() -= AScalar(0.5)*((Z*Z).rowwise().sum());

  // half-Cauchy prior on alpha
  prior_i.array() += AScalar(M_LN2 - LOG_PI) - log(AScalar(1)+(alpha.array()/alpha_prior).pow(2));
}




AScalar cheese::eval_prior() {
  
  VectorXA log_prior_i(N);
  eval_prior_i(log_prior_i);

 // add log Jacobian for log_alpha -> alpha
  log_prior_i += log_alpha;

  AScalar log_prior = log_prior_i.sum();
  assert(my_finite(log_prior));
  return(log_prior);
}




AScalar cheese::eval_hyperprior() {
  
  MatrixXA mu_prior(1,1);
  MVN_logpdf(mu, mu_prior_mean, mu_prior_chol_prec, mu_prior, true);
  AScalar G_prior = Wishart_logpdf(chol_G, nu, chol_A);
  AScalar logJacG = k * M_LN2;
  for (int i=0; i<k; i++) {
    logJacG += (k-i+2)*log(chol_G(i,i));
  }
 
  AScalar res = mu_prior(0,0) + G_prior + logJacG;
  assert(my_finite(res));  

  return(res);

}



template<typename Tpars>
typename Tpars::Scalar cheese::eval_LL(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  typename Tpars::Scalar LL = eval_LL();
  
  return(LL);
 
}



template<typename Tpars>
typename Tpars::Scalar cheese::eval_prior(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  typename Tpars::Scalar prior = eval_prior();
  
  return(prior);
 
}



template<typename Tpars>
typename Tpars::Scalar cheese::eval_hyperprior(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  typename Tpars::Scalar hyperprior = eval_hyperprior();
  
  return(hyperprior);
 
}


AScalar cheese::eval_LLi() {
  return(AScalar(0));
}


#endif
