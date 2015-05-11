
rm(list=ls())
gc()

library(Rcpp)
library(RcppEigen)
library(Matrix)
library(bayesm)
library(reshape2)
library(plyr)
library(MCMCpack)


set.seed(1234)

file.path <- "~/Documents/gds-mksc/conjugate/"

mod.name <- "conjugate"
lib <- dyn.load(paste(mod.name,".so",sep=""))
module <- Module(mod.name, PACKAGE=lib)

mode.file <- paste0(file.path,"results/",mod.name,"_mode.Rdata")
save.file <- paste0(file.path,"results/",mod.name,"_gibbs.Rdata")

load(mode.file)
Y <- L$data$Y
Ymean <- rowMeans(Y)

N <- NROW(Y)
T <- NCOL(Y)

nvars <- N + 3

ndraws <- 2000
report <- 500
warmup <- 1000

draws <- matrix(0,ndraws-warmup,nvars)
LL <- vector("numeric",length=ndraws-warmup)


theta.start <- rep(0,N)
mu.start <- 1
sig2.start <- .1
tau2.start <- .1

theta <- theta.start
mu <- mu.start
sig2 <- sig2.start
tau2 <- tau2.start

##start <- c(theta, mu, log(sig2)/2, log(tau2)/2)
start <- opt$solution

cl <- eval(parse(text=paste0("new(module$",mod.name,", L)")))
cl$record.tape(start)


sig.df <- T*N

for (d in 1:ndraws) {

    if (d %% report ==0) {
        cat("Iteration: ",d,"\n")
    }
    
    theta.prec <- 1/tau2 + T/sig2
    theta.sd <- 1/sqrt(theta.prec)
    theta.mean <- (mu/tau2 + T*Ymean/sig2)/theta.prec
    theta <- rnorm(N,theta.mean, theta.sd)

    mu.mean <- mean(theta)
    mu.sd <- sqrt(tau2/N)
    mu <- rnorm(1, mu.mean, mu.sd)

    sig2.scale <- mean((Y-theta)^2) ## theta should repeat by column
    tmp <- rchisq(1,sig.df)
    sig2 <- sig.df*sig2.scale/tmp

    tau2.scale <- sum((theta-mu)^2)/(N-1)
    tmp <- rchisq(1, N-1)
    tau2 <- (N-1)*tau2.scale/tmp

    if (d > warmup) {
        res <- c(theta, mu, log(sig2)/2, log(tau2)/2)
        draws[d-warmup,] <- res
        LL[d-warmup] <- cl$get.f(res)
    }
}

save(draws, LL, file=save.file)






