rm(list=ls())
gc()

library(Rcpp)
library(RcppEigen)
library(doParallel)
library(trustOptim)
library(numDeriv)
library(sparseMVN)
library(bayesGDS)
library(plyr)
library(Matrix)

run.par <- TRUE
if(run.par) registerDoParallel(cores=1) else registerDoParallel(cores=1)

chain <- 5

seed.id <- 1234*chain
set.seed(seed.id)
mod.name <- "conjugate"
file.path <- "~/Documents/gds-mksc/conjugate"



rmvn.sparse.wrap <- function(n.draws, params) {
  res <- rmvn.sparse(n.draws, params$mean, params$CH, prec=TRUE)
  return(res)
}

dmvn.sparse.wrap <- function(d, params) {
  res <- sparseMVN:::dmvn.sparse(d, params$mean, params$CH, prec=TRUE)
  return(res)
}

extract.gds <- function(B) {

  counts <- as.vector(laply(B, function(X) return(X$counts)))
  gt.1 <- as.vector(laply(B, function(X) return(X$gt.1)))
  log.thresholds <- as.vector(laply(B, function(X) return(X$log.thresholds)))
  log.post.dens <- as.vector(laply(B, function(X) return(X$log.post.dens)))
  log.prop.dens <- as.vector(laply(B, function(X) return(X$log.prop.dens)))
  log.phi <- as.vector(laply(B, function(X) return(X$log.phi)))
  draws <- foreach(i=1:length(B), .combine=rbind,
                   .multicombine=TRUE) %do% return(B[[i]]$draws)
  acc.rate <- 1/mean(counts,na.rm=TRUE)
  
  res <- list(draws=draws,counts=counts, gt.1=gt.1,
              log.thresholds=log.thresholds, log.post.dens=log.post.dens,
              log.prop.dens=log.prop.dens, log.phi=log.phi)
  return(res)
 
} 

mode.file <- paste0(file.path,"/results/",mod.name,"_mode.Rdata")
save.file <- paste0(file.path,"/results/",mod.name,"_timetest.Rdata")

M <- 1000
sc <- .98

load(mode.file)

var.names <- names(opt$solution)

lib <- dyn.load(paste(mod.name,".so",sep=""))
module <- Module(mod.name,PACKAGE=lib)

pm <- as.vector(opt$solution)
hs <- opt$hessian
nvars <- length(pm)
CH <- Cholesky(-sc*opt$hessian)
log.c1 <- opt$fval

cl <- eval(parse(text=paste0("new(module$",mod.name,", L)")))
cl$record.tape(pm)

get.f <- function(P, ...) return(cl$get.f(P))

prop.params <- list(mean = pm,
                    CH = CH
                    )

log.c2 <- dmvn.sparse.wrap(pm, prop.params)
log.const <- log.c1 - log.c2

cat("Collecting GDS Proposal Draws\n")

cat("Sampling MVN proposals\n")
st1 <- system.time(draws.m <- as(rmvn.sparse.wrap(M,prop.params),"matrix"))

cat("Evaluating log posterior\n")
st2 <- system.time(log.post.m <- aaply(draws.m, 1,get.f))

cat("Evaluating MVN density\n")
st3 <- system.time(log.prop.m <- dmvn.sparse.wrap(draws.m,params=prop.params))


save(st1, st2, st3, M, file=save.file)
