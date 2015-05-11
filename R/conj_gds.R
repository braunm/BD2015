
rm(list=ls())
gc()

library(Rcpp)
library(RcppEigen)
library(Matrix)
library(bayesm)
library(reshape2)
library(plyr)
library(bayesGDS)
library(foreach)
library(doParallel)
library(sparseMVN)

rmvn.sparse.wrap <- function(n.draws, params) {
  res <- rmvn.sparse(n.draws, params$mean, params$CH, prec=TRUE)
  return(res)
}

dmvn.sparse.wrap <- function(d, params) {
  res <- sparseMVN::dmvn.sparse(d, params$mean, params$CH, prec=TRUE)
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


run.par <- TRUE
if(run.par) registerDoParallel(cores=20) else registerDoParallel(cores=1)

n.draws <- 240
batch.size <- 12
est.acc.rate <- 0.01
max.AR.tries <- 5000000


set.seed(12345)


sc <- .98
M <- 70000

file.path <- "~/Documents/gds-mksc/conjugate/"

mod.name <- "conjugate"
lib <- dyn.load(paste(mod.name,".so",sep=""))
module <- Module(mod.name, PACKAGE=lib)

mode.file <- paste0(file.path,"results/",mod.name,"_mode.Rdata")
save.file <- paste0(file.path,"results/",mod.name,"_gds2.Rdata")

load(mode.file)
Y <- L$data$Y

post.mode <- unlist(sol)

N <- NROW(Y)
T <- NCOL(Y)

nvars <- N + 3


cl <- eval(parse(text=paste0("new(module$",mod.name,", L)")))
cl$record.tape(post.mode)

get.f <- function(P) return(cl$get.f(P))
log.c1 <- opt$fval

chol.hess <- Cholesky(-sc*opt$hessian)
prop.params <- list(mean = post.mode,
                    CH = chol.hess
                    )

log.c2 <- dmvn.sparse.wrap(matrix(post.mode,nrow=1), prop.params)
log.const <- log.c1 - log.c2

cat("Collecting GDS Proposal Draws\n")

st1 <- Sys.time()

draws.m <- as(rmvn.sparse.wrap(M,prop.params),"matrix")
log.post.m <- aaply(draws.m, 1,get.f, .parallel=TRUE)
log.prop.m <- dmvn.sparse.wrap(draws.m,params=prop.params)
log.phi <- log.post.m - log.prop.m +log.c2 - log.c1

cat("Are any log.phi > 0?  ",any(log.phi>0),"\n")
if (all(log.phi <= 0)) {
    
    n.batch <- floor(n.draws / batch.size)
    
    cat("Generating DS draws - accept-reject phase\n")
    rs <- .Random.seed
    draws <- foreach(i=1:n.batch, .inorder=FALSE) %dopar% sample.GDS(
                                      n.draws=batch.size,
                                      log.phi=log.phi,
                                      post.mode=matrix(post.mode,nrow=1),
                                      fn.dens.post = get.f,
                                      fn.dens.prop = dmvn.sparse.wrap,
                                      fn.draw.prop = rmvn.sparse.wrap,
                                      prop.params = prop.params,
                                      max.tries=max.AR.tries,
                                      report.freq=100,
                                      announce=TRUE,
                                      thread.id=i,
                                      seed=i)
    
    st2 <- Sys.time()
    st.diff <- difftime(st2,st1,units="mins")
    dd <- extract.gds(draws)
    if (any(is.na(dd$counts))) {
        LML <- NA
    } else {
        LML <- get.LML(counts=dd$counts,log.phi=log.phi,
                       post.mode=post.mode,
                       fn.dens.post=get.f,
                       fn.dens.prop=dmvn.sparse.wrap,
                       prop.params=prop.params)
        
    }
    ## Section H:  Compute log marginal likelihood
    
    acc.rate <- 1/mean(dd$counts)

    dd$LML <- LML
    dd$acc.rate <- acc.rate   
}

save(dd, file=save.file)







