#ifndef __BINARY_EVAL
#define __BINARY_EVAL

#include <mb_base.h>
#include <Eigen/Core>
#include <cppad_atomics.h>
#include <MVN_AD.h>
#include <wish_AD.h>


using Eigen::Matrix;
using Eigen::MatrixBase;
using Eigen::Dynamic;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;


AScalar log_choose(const AScalar& n, const AScalar& k) {

  AScalar res = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
  return(res);

}

class binary  {

  typedef Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
  typedef Matrix<AScalar, Dynamic, 1> VectorXA;
  typedef Matrix<AScalar, 1, Dynamic> RowVectorXA;

 public:
  
  binary(const List&);
  
  template<typename Tpars>
    AScalar eval_f(const MatrixBase<Tpars>&);
  
  AScalar eval_LL(); 
  
  template<typename Tpars>
    AScalar eval_LL(const MatrixBase<Tpars>&);
  
  AScalar eval_LLi(); 

  template<typename TL>
    void eval_LLi(const MatrixBase<TL>&);

  template<typename Tpars, typename TL>
    void eval_LLi(const MatrixBase<Tpars>&, const MatrixBase<TL>&);

  AScalar eval_prior();
  AScalar eval_hyperprior();
  
  template<typename Tpars>
    AScalar eval_prior(const MatrixBase<Tpars>&);
  
 
  template<typename Tpars>
    AScalar eval_hyperprior(const MatrixBase<Tpars>&);
  
  template<typename Tpars>
    void unwrap_params(const MatrixBase<Tpars>&);

  
 protected:
  
  // data and priors.

  VectorXA Y;
  MatrixXA X;
  MatrixXA mu_prior_mean;
  MatrixXA mu_prior_chol_prec;
  double nu;
  MatrixXA chol_S;
  AScalar T;
  int N, k; // number of observations per unit, number of units, number of regressors
  
  // declare parameters

  MatrixXA B;
  VectorXA mu;
  MatrixXA chol_G;
  VectorXA logit_p; // N x 1
  VectorXA binom_i;

};


binary::binary(const List& params) {

  List & pars = static_cast<List &>(const_cast<List &>(params));

  List data = Rcpp::as<List>(pars["data"]);
  List priors = Rcpp::as<List>(pars["priors"]);

  // Map R to Rcpp
  NumericVector Y_((SEXP)data["Y"]);
  NumericMatrix X_((SEXP)data["X"]);
  T = Rcpp::as<double>((SEXP)data["T"]);
  NumericVector mu_prior_mean_((SEXP)priors["mu.prior.mean"]);
  NumericMatrix mu_prior_chol_prec_((SEXP)priors["mu.prior.chol.prec"]);
  nu = Rcpp::as<double>((SEXP)priors["nu"]);
  NumericMatrix chol_S_((SEXP)priors["chol.S"]);
  
  N = X_.ncol();
  k = X_.nrow();

  // Map Rcpp to Eigen Scalar
  Map<MatrixXd> Xd = MatrixXd::Map(X_.begin(),k,N);
  Map<VectorXd> Yd = VectorXd::Map(Y_.begin(),N);
  Map<MatrixXd> mu_prior_meand = MatrixXd::Map(mu_prior_mean_.begin(),k,1);
  Map<MatrixXd> mu_prior_chol_precd = MatrixXd::Map(mu_prior_chol_prec_.begin(),k,k);
  Map<MatrixXd> chol_Sd = MatrixXd::Map(chol_S_.begin(),k,k);

  // copy Eigen Scalar Map to Eigen AScalar objects
  X = Xd.cast<AScalar>();
  Y = Yd.cast<AScalar>();
  mu_prior_mean = mu_prior_meand.cast<AScalar>();
  mu_prior_chol_prec = mu_prior_chol_precd.cast<AScalar>();
  chol_S = chol_Sd.cast<AScalar>();


  // reserve space for parameters
  chol_G.resize(k,k);
  B.resize(k,N);
  mu.resize(k);
  logit_p.resize(N);

  binom_i.resize(N);
  for (int i=0; i<N; i++) {
    binom_i(i) = log_choose(T, Y(i));
  }

}



template<typename Tpars>
void binary::unwrap_params(const MatrixBase<Tpars>& P)
{

  // here we are copying the parameters.  Could we map instead? Doubtful because of the cast. 

  B = MatrixXA::Map(P.derived().data(), k, N);
  mu = P.derived().segment(k*N,k);

  chol_G.setZero();   
  int idx = N*k+k;
  for (int j=0; j<k; j++) {
    chol_G(j,j) = exp(P(idx));
    idx++;
    for (int i=j+1; i<k; i++) {
      chol_G(i,j) = P(idx); // fill lower triangle only
      idx++;
    }
  }
}



template<typename Tpars>
AScalar binary::eval_f(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);

  AScalar LL = eval_LL();
  AScalar prior = eval_prior();
  AScalar hyperprior = eval_hyperprior();
  
  AScalar f = LL + prior + hyperprior;

  return(f);
}



template<typename TL>
void binary::eval_LLi(const MatrixBase<TL>& LLi_)
{ 


  MatrixBase<TL>& LLi = const_cast<MatrixBase<TL>& >(LLi_);

  logit_p = (B.array() * X.array()).matrix().colwise().sum().transpose();
  LLi = binom_i;
  LLi.array() += Y.array() * logit_p.array();

  for (int i=0; i<N; i++) {
    LLi(i) +=  CppAD::CondExpGt(logit_p(i),AScalar(700),
				-T*logit_p(i),
				-T*log1p(exp(logit_p(i)))
				);
  }			   
}

template<typename Tpars, typename TL>
void binary::eval_LLi(const MatrixBase<Tpars>& P,
						const MatrixBase<TL>& LLi_) {
 
  MatrixBase<TL>& LLi = const_cast<MatrixBase<TL>& >(LLi_); 
  unwrap_params(P);
  eval_LLi(LLi);
}


AScalar binary::eval_LL() {

  VectorXA LLi(N);
  eval_LLi(LLi);
  AScalar LL = LLi.sum();
  return (LL);

}

AScalar binary::eval_prior() {
  
  VectorXA B_prior(N);
  MVN_logpdf(B, mu, chol_G, B_prior, true);
  return(B_prior.sum());
}

AScalar binary::eval_hyperprior() {
 
  VectorXA mu_prior(1);
  MVN_logpdf(mu, mu_prior_mean, mu_prior_chol_prec, mu_prior, true);
  
  AScalar G_prior = Wishart_logpdf(chol_G, nu, chol_S);

  AScalar logJacG = k * M_LN2;
  for (int i=1; i<=k; i++) {
    logJacG += AScalar(k-i+2)*log(chol_G(i-1,i-1));
  }
 
  AScalar res = mu_prior(0) + G_prior + logJacG;
  
 return(res);

}


template<typename Tpars>
AScalar binary::eval_LL(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  AScalar LL = eval_LL();
  
  return(LL);
 
}


template<typename Tpars>
AScalar binary::eval_prior(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  AScalar prior = eval_prior();
  
  return(prior);
 
}


template<typename Tpars>
AScalar binary::eval_hyperprior(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  AScalar hyperprior = eval_hyperprior();
  
  return(hyperprior);
 
}

#endif
