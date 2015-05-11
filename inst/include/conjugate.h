#ifndef __CONJUGATE
#define __CONJUGATE

#include <mb_base.h>
#include <except.h>
#include <utilfuncs.h>
#include <Eigen/Core>

using Eigen::Matrix;
using Eigen::MatrixBase;
using Eigen::Dynamic;
using Eigen::Block;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

template<typename T>
bool my_finite(const T& x) {
  return( (abs(x) <= __DBL_MAX__ ) && ( x == x ) );
}

class conjugate  {

  typedef int Index; // Tpars is the type of the parameter vector

  typedef Matrix<AScalar, Dynamic, Dynamic> MatrixXA;
  typedef Matrix<AScalar, Dynamic, 1> VectorXA;

 public:
  
  conjugate(const List&);
  
  template<typename Tpars>
    AScalar eval_f(const MatrixBase<Tpars>&);
  
  AScalar eval_LL(); 
  
  template<typename Tpars>
    typename Tpars::Scalar eval_LL(const MatrixBase<Tpars>&);
  
  AScalar eval_LLi(); 


  template<typename TL>
    void eval_LLi(const MatrixBase<TL>&);

  template<typename Tpars, typename TL>
    void eval_LLi(const MatrixBase<Tpars>&, const MatrixBase<TL>&);

  template<typename TL>
    void eval_prior_i(const MatrixBase<TL>&);

  template<typename Tpars, typename TL>
    void eval_prior_i(const MatrixBase<Tpars>&, const MatrixBase<TL>&);

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

  MatrixXA Y;

  Index N, T;
  
  VectorXA theta;
  AScalar mu;
  AScalar log_sig;
  AScalar log_tau;
  AScalar s2;
  AScalar t2;

};



conjugate::conjugate(const List& params)
{


  List & pars = static_cast<List &>(const_cast<List &>(params));

  List data = Rcpp::as<List>(pars["data"]);
  
  // Map R to Rcpp
  NumericMatrix Y_((SEXP)data["Y"]); // N x T

  T = Y_.ncol();
  N = Y_.nrow();

  Map<MatrixXd> Y_d = MatrixXd::Map(Y_.begin(),N,T);
  Y = Y_d.cast<AScalar>();

  // reserve space for parameters
  theta.resize(N);

}


template<typename Tpars>
void conjugate::unwrap_params(const MatrixBase<Tpars>& P)
{

  // here we are copying the parameters.  Could we map instead? Doubtful because of the cast. 

  theta = P.derived().head(N);
  mu = P(N);
  log_sig = P(N+1);
  log_tau = P(N+2);
  s2 = exp(2*log_sig);
  t2 = exp(2*log_tau);
}




template<typename Tpars>
AScalar conjugate::eval_f(const MatrixBase<Tpars>& P) {


  unwrap_params(P);

  AScalar LL = eval_LL();
  AScalar prior = eval_prior();
  AScalar hyperprior = eval_hyperprior();
  
  AScalar f = LL + prior + hyperprior;
  
  return(f);
}




template<typename TL>
void conjugate::eval_LLi(const MatrixBase<TL>& LLi_)
{ 

  MatrixBase<TL>& LLi = const_cast<MatrixBase<TL>& >(LLi_);

  LLi.setConstant(-T*(log_sig + M_LN_SQRT_2PI));

  MatrixXA Z = (Y.colwise() - theta).cwiseAbs2();
  Z.array() = 0.5 * Z.array() / s2;
  LLi -= Z.rowwise().sum();

}

AScalar conjugate::eval_LL()
{ 

  VectorXA LLi(N); //vector to store individual-level likelihoods
  eval_LLi(LLi);
  AScalar LL = LLi.sum();
  return(LL);

}


template<typename TL>
void conjugate::eval_prior_i(const MatrixBase<TL>& prior_i_) {
  
  MatrixBase<TL>& prior_i = const_cast<MatrixBase<TL>& >(prior_i_);

  prior_i.setConstant(-log_tau - M_LN_SQRT_2PI);
  prior_i.array() -= 0.5 * (theta.array()-mu).square()/t2;

}




AScalar conjugate::eval_prior() {
  
  VectorXA log_prior_i(N);
  eval_prior_i(log_prior_i);
  AScalar log_prior = log_prior_i.sum();
  assert(my_finite(log_prior));
  return(log_prior);
}




AScalar conjugate::eval_hyperprior() {
  
  // Jacobian, since we place a uniform prior on tau,
  //   not log_tau.
 
  AScalar res = log_tau;
  assert(my_finite(res));  

  return(res);

}



template<typename Tpars>
typename Tpars::Scalar conjugate::eval_LL(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  typename Tpars::Scalar LL = eval_LL();
  
  return(LL);
 
}



template<typename Tpars>
typename Tpars::Scalar conjugate::eval_prior(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  typename Tpars::Scalar prior = eval_prior();
  
  return(prior);
 
}



template<typename Tpars>
typename Tpars::Scalar conjugate::eval_hyperprior(const MatrixBase<Tpars>& P) {
  
  unwrap_params(P);
  typename Tpars::Scalar hyperprior = eval_hyperprior();
  
  return(hyperprior);
 
}


AScalar conjugate::eval_LLi() {
  

  return(AScalar(0));
}



template<typename Tpars, typename TL>
void conjugate::eval_LLi(const MatrixBase<Tpars>& P,
						const MatrixBase<TL>& LLi_) {
 
  MatrixBase<TL>& LLi = const_cast<MatrixBase<TL>& >(LLi_); 
  unwrap_params(P);
  eval_LLi(LLi);

}



template<typename Tpars, typename TL>
void conjugate::eval_prior_i(const MatrixBase<Tpars>& P,
						const MatrixBase<TL>& log_prior_i_) {
 
  MatrixBase<TL>& log_prior_i = const_cast<MatrixBase<TL>& >(log_prior_i_); 
  unwrap_params(P);
  eval_prior_i(log_prior_i);

}






#endif
