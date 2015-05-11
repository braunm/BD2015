
rm(list=ls())
gc()

library(Rcpp)
library(RcppEigen)
library(Matrix)
library(bayesm)
library(reshape2)
library(plyr)
library(trustOptim)

set.seed(123)

file.path <- "~/Documents/gds-mksc/conjugate/"

mod.name <- "conjugate"
lib <- dyn.load(paste(mod.name,".so",sep=""))
module <- Module(mod.name, PACKAGE=lib)

save.file <- paste0(file.path,"results/",mod.name,"_mode.Rdata")
data.file <- paste0(file.path,"data/",mod.name,".Rdata")
load(data.file)

N <- NROW(Y)
T <- NCOL(Y)

data <- list(Y=Y)

priors <- NULL

L <- list(data=data, priors=priors)

method <- "Sparse"
control.list <- list(start.trust.radius=5,
                     stop.trust.radius = 1e-5,
                     prec=1e-7,
                     report.freq=1L,
                     report.level=4L,
                     report.precision=1L,
                     maxit=200L,
                     function.scale.factor = as.numeric(-1),                           
                     preconditioner=1L
                     ) 

nvars <- N + 3

start.list <- as.relistable(list(theta=rep(0,N),
                                 mu=2,
                                 log_sig=0,
                                 log_tau=0)
                            )

start <- unlist(start.list)

cl <- eval(parse(text=paste0("new(module$",mod.name,", L)")))
cl$record.tape(start)
cat("Initializing Hessian\n")
cl$hessian.init.nopattern(start)


get.f <- function(P) return(cl$get.f(P))
get.df <- function(P) return(cl$get.fdf(P)$grad)
get.fdf <- function(P) return(cl$get.fdf(P))
get.hessian <- function(P) return(cl$get.hessian.sparse(P))

opt <- trust.optim(start, fn=get.f,
                   gr = get.df,
                   hs = get.hessian,  ## used only for Sparse method
                   method = "Sparse",
                   control = control.list
                   )

sol <- relist(opt$solution,skeleton=start.list)
save(sol, opt, L, file=save.file)
