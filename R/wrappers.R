#' @name sparseMVN_wrappers
#' @title sparseMVN wrappers
#' @description wrapper to sparseMVN package
#' @param n.draws Number of draws
#' @param params parameter list
#' @return random MVN samples
rmvn.sparse.wrap <- function(n.draws, params) {
  res <- rmvn.sparse(n.draws, params$mean, params$CH, prec=TRUE)
  return(res)
}


#' @param d Values at which to evaluate density
#' @return MVN densities
#' @rdname sparseMVN_wrappers
dmvn.sparse.wrap <- function(d, params) {
  res <- sparseMVN::dmvn.sparse(d, params$mean, params$CH, prec=TRUE)
  return(res)
}
