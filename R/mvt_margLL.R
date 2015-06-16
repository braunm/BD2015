#' @title True marginal likelihood for MVT example
#' @param data data
#' @param hyps hyperpriors
#' @export
get.marg.LL <- function(data, hyps) {

  b0 <- hyps$b0
  chol.inv.V0 <- hyps$chol.inv.V0
  inv.V0 <- tcrossprod(chol.inv.V0)
  V0 <- solve(inv.V0)
  r0 <- hyps$r0
  s0 <- hyps$s0

  y <- data$Y
  x <- data$X
  n <- length(y)
  xtx <- tcrossprod(x)
  k <- NROW(x)+1

  yxb <- y - b0 %*% x

  Q <- solve(diag(n) + t(x) %*% V0 %*% x)
  chol.D <- t(chol(diag(k) + V0 %*% xtx))
  log.D <- 2*sum(log(diag(chol.D)))

  log.const <- -n*log(pi*s0)/2 - log.D/2 +  lgamma((r0+n)/2) - lgamma(r0/2)

  W <- 1 + (yxb %*% Q %*% t(yxb))/s0
  log.const - (r0+n)*log(W)/2

}

