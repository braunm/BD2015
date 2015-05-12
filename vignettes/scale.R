
library(plyr)
library(Matrix)
library(mvtnorm)
library(numDeriv)
library(trustOptim)
library(dplyr)
library(ggplot2)


nreps <- 50
N.vec <- c(500,1000,2000,5000,10000,15000,20000,30000,40000,50000)
nx <- length(N.vec)
k <- 3

kk12 <- k*(k+1)/2

T <- 52
B.mean.true <- seq(-10,10,length=k)
B.cov.true <- 0.1*diag(k)
x.mean <- rep(0,k-1)
x.cov <- diag(k-1)

R <- 200

res <- array(dim=c(nx,10,nreps))
dimnames(res) <- list(N=N.vec,
                      time=c("logpost","gradient","hessian","mode",
                             "rnorm","chol","solve","add","mult","tot_mvn"),
                      rep=1:nreps)

for (i in 1:nx) {

    N <- N.vec[i]
    cat("N = ",N,"\t")
    Nk <- N*k
    nvars <- Nk + k + kk12

    sx <- Matrix(0,nrow=nvars,ncol=nvars)
    sx[1:Nk,1:Nk] <- bdiag(llply(1:N,function(i) return(tril(Matrix(TRUE,nrow=k,ncol=k)))))
    sx[(Nk+1):nvars,1:Nk] <- TRUE
    sx[1:Nk,(Nk+1):nvars] <- TRUE
    sx[(Nk+1):nvars,(Nk+1):nvars] <- tril(Matrix(TRUE,nrow=k+kk12, ncol=k+kk12))
    sx <- as(sx, "TsparseMatrix")

    hess.struct <- list(rows=sx@i+1, cols=sx@j+1)

    set.seed(round(123*k*N/100))
    B <- t(rmvnorm(N, mean=B.mean.true, sigma=B.cov.true)) ## k x N
    X <- rbind(rep(1,N),t(rmvnorm(N, mean=x.mean, sigma=x.cov))) ## k x N
    XB <- colSums(X * B)
    log.p <- XB - log1p(exp(XB))
    Y <- laply(log.p, function(q) return(rbinom(1,T,exp(q))))

    cat("Building data and priors list\n")
    data <- list(Y=Y, X=X, T=T)
    priors <- list(mu.prior.mean=rep(0,k),
                   mu.prior.chol.prec=t(chol(.0001*diag(k))),
                   chol.S=t(chol(diag(k))),
                   nu=k+6
                   )

    DL <- list(data=data, priors=priors)


    chol.G.start <- diag(k)
    vech.chol.G.start <- t(t(chol.G.start[!upper.tri(chol.G.start)]))

    start.list <- as.relistable(list(beta=matrix(0,k,N),
                                     mu=rep(0,k),
                                     chol.G=vech.chol.G.start
                                     )
                                )

    start <- unlist(start.list)

    cl <- new("binary", DL)
    cl$record.tape(start)
    cl$init.hessian(start)

    get.f <- function(P) return(cl$get.f(P))
    get.df <- function(P) return(cl$get.df(P))
    get.hessian <- function(P) return(cl$get.hessian.sparse(P))


    for (rep in 1:nreps) {

        t4 <- system.time(
            opt <- trust.optim(start, fn=get.f,
                               gr = get.df,
                               hs = get.hessian,
                               method = "Sparse",
                               control =  list(start.trust.radius=5,
                                               stop.trust.radius = 1e-5,
                                               prec=1e-7,
                                               report.freq=100L,
                                               report.level=0L,
                                               report.precision=1L,
                                               maxit=1000L,
                                               function.scale.factor = -1,
                                               preconditioner=0L,
                                               trust.iter=20000
                                               )
                               )
            )

        t1 <- system.time(f <- get.f(opt$solution))
        t2 <- system.time(df <- get.df(opt$solution))
        t3 <- system.time(hess <- get.hessian(opt$solution))

        res[i,"logpost",rep] <- t1["elapsed"]
        res[i,"gradient",rep] <- t2["elapsed"]
        res[i,"hessian",rep] <- t3["elapsed"]
        res[i,"mode",rep] <- t4["elapsed"]

        mu <- opt$solution
        t5 <- system.time(z <- rnorm(R*nvars))
        dim(z) <- c(nvars,R)
        t6 <- system.time(C <- Cholesky(-hess))
        ex <- expand(C)
        L <- ex$L
        P <- ex$P
        t7 <- system.time(y <- solve(t(L),z))
        t8 <- system.time(y <- y+mu)
        t9 <- system.time(L %*% z)

        res[i,"rnorm", rep] <- t5["elapsed"]
        res[i,"chol", rep] <- t6["elapsed"]
        res[i,"solve", rep] <- t7["elapsed"]
        res[i,"add", rep] <- t8["elapsed"]
        res[i,"mult", rep] <- t9["elapsed"]
        res[i,"tot_mvn", rep] <- sum(res[i,5:9, rep])
    }
    rm(cl)
    gc()
}


D <- reshape2::melt(res) %>%
  mutate(plotID=as.factor(ifelse(time %in% c("logpost","gradient"),"logpost/gradient",
                                 ifelse(time == "hessian","hessian",
                                        ifelse(time == "tot_mvn","total_mvn","mvn_steps"))
                                 )
                          )
         )


save(D,file="vignettes/results/scale.Rdata")

theme_set(theme_bw())
P1 <- D %>%
  group_by(N, time, plotID) %>%
  summarize(seconds=mean(value)) %>%
  ggplot(aes(x=N, y=seconds, color=time)) %>%
  + geom_line() %>%
  + geom_point() %>%
  + facet_wrap(~plotID, scales="free")

pdf(file="vignettes/results/scale.pdf")
P1
dev.off()


