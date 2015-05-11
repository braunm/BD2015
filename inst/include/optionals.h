double get_f_direct(const NumericVector& P_) {

  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  
  VectorXA P = VectorXd::Map(P_.begin(), nvars()).template cast<AScalar>();
  AScalar fr = model -> eval_f(P); // inherited from base class
  double res = Value(fr);
  return(res);
}

double get_LL(const NumericVector& P_) {
  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  VectorXA P = VectorXd::Map(P_.begin(), nvars()).template cast<AScalar>();  
  AScalar fr = model -> eval_LL(P); // inherited from base class
  double res = Value(fr);
  return(res);
}

double get_hyperprior(const NumericVector& P_) {
  if (!tape_ready)
    throw MyException("tape not ready",__FILE__,__LINE__);

  VectorXA P = VectorXd::Map(P_.begin(), nvars()).template cast<AScalar>();  
  AScalar fr = model -> eval_hyperprior(P); // inherited from base class
  double res = Value(fr);
  return(res);
}


