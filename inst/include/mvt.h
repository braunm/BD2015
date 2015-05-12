#ifndef MVT_EVAL
#define MVT_EVAL



#include <mb_base.h>
#include <except.h>
#include <utilfuncs.h>
#include <Eigen/Core>
#include <cppad_atomics.h>
#include <Eigen/Core>
#include <MVN_AD.h>
//#include "/Users/braunm/Documents/R_packages/CppADutils/inst/include/MVN_AD.h"

using Eigen::Matrix;
using Eigen::MatrixBase;
using Eigen::Dynamic;
using Eigen::Map;
using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::VectorXd;




class mvt {

  typedef int Index; // Tpars is the type of the parameter vector

  typedef Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
  typedef Matrix<AScalar, Dynamic, 1> VectorXA;
  typedef Matrix<AScalar, 1, Dynamic> RowVectorXA;


 public:
  
  mvt(const List&);
  
  template<typename Tpars>
    AScalar eval_f(const MatrixBase<Tpars>&);
  
  AScalar eval_LL(); 

  
  template<typename Tpars>
    AScalar eval_LL(const MatrixBase<Tpars>&);
  


  AScalar eval_prior();

  
  template<typename Tpars>
    AScalar eval_prior(const MatrixBase<Tpars>&);
  
 
  
  template<typename Tpars>
    void unwrap_params(const MatrixBase<Tpars>&);


  
 protected:
  
  // data and priors
  
  MatrixXA Y;
  MatrixXA X;
  VectorXA b0;
  MatrixXA chol_inv_V0;
  AScalar r0, s0;
  
  Index N, T, k; // number of observations per unit, number of units, number of regressors

  // parameters

  VectorXA B;
  AScalar logsig2;
  AScalar sig;
  MatrixXA sigMat;
    
};


mvt::mvt(const List& params) {

  using Rcpp::List;
  using Rcpp::NumericMatrix;
  using Rcpp::NumericVector;
  using Rcpp::as;

  List & pars = static_cast<List &>(const_cast<List &>(params));

  List data = as<List>(pars["data"]);
  List priors = as<List>(pars["priors"]);

  NumericMatrix Y_((SEXP)data["Y"]);
  NumericMatrix X_((SEXP)data["X"]);

  NumericVector b0_((SEXP)priors["b0"]);
  NumericMatrix chol_inv_V0_((SEXP)priors["chol.inv.V0"]);
  r0 = as<double>((SEXP)priors["r0"]);
  s0 = as<double>((SEXP)priors["s0"]);
  
  N = Y_.ncol();
  T = Y_.nrow();
  k = X_.nrow();
  T = 1;

  Map<MatrixXd> Xd = MatrixXd::Map(X_.begin(),k,N);
  Map<MatrixXd> Yd = MatrixXd::Map(Y_.begin(),T,N);
  Map<VectorXd> b0d = VectorXd::Map(b0_.begin(),k);
  Map<MatrixXd> chol_inv_V0d = MatrixXd::Map(chol_inv_V0_.begin(),k,k);
 
  X = Xd.cast<AScalar>();
  Y = Yd.cast<AScalar>();
  b0 = b0d.cast<AScalar>();
  chol_inv_V0 = chol_inv_V0d.cast<AScalar>();

  B.resize(k);
  sigMat = MatrixXA::Identity(T,T);

}


template<typename Tpars>
void mvt::unwrap_params(const MatrixBase<Tpars>& P) {

  B = VectorXA::Map(P.derived().data(), k);
  logsig2 = P(k);
  sig = exp(logsig2/2);
  sigMat.diagonal().setConstant(1/sig);
}



AScalar mvt::eval_LL() {
  

  RowVectorXA y_mean = B.transpose() * X;
  MatrixXA LLi(1,y_mean.size());

  
  MVN_logpdf(Y, y_mean, sigMat, LLi, true);
  return(LLi.sum());
}



AScalar mvt::eval_prior() {
  
  MatrixXA B_prior_i(X.rows(), 1);
  MVN_logpdf(B, b0, chol_inv_V0/sig, B_prior_i, true); // error?
  AScalar sig2 = exp(logsig2);
  AScalar sigprior = (r0/2)*(log(s0)-M_LN2) - lgamma(r0/2) - (1+r0/2)*logsig2 - s0/(2*sig2);
  AScalar sigJac = logsig2; // log Jacobian of transform
  AScalar prior = B_prior_i.sum() + sigprior + sigJac;


  return(prior);
}


template<typename Tpars>
AScalar mvt::eval_f(const MatrixBase<Tpars>& P) {

  unwrap_params(P);
  AScalar LL = eval_LL();
  AScalar prior = eval_prior();

  AScalar f = LL + prior;
  return(f);
}


template<typename Tpars>
AScalar mvt::eval_LL(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  AScalar LL = eval_LL();
  
  return(LL);
 
}


template<typename Tpars>
AScalar mvt::eval_prior(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  AScalar prior = eval_prior();
  
  return(prior);
 
}



#endif


 
