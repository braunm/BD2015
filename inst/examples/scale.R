
library(plyr)
library(Matrix)
library(mvtnorm)
library(numDeriv)
library(trustOptim)
library(dplyr)
library(ggplot2)
library(tidyr)
library(doParallel)
library(BD2015)

run.par <- TRUE
if(run.par) registerDoParallel(cores=10) else registerDoParallel(cores=1)

nreps <- 20
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

run_times <- function(N) {
    output <- NULL
    cat("N = ",N,"\n")
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
                                               report.freq=1L,
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

        mu <- opt$solution
        t5 <- system.time(z <- rnorm(R*nvars))
        dim(z) <- c(nvars,R)
        t6 <- system.time(C <- Cholesky(-hess))
        ex <- Matrix::expand(C)
        L <- ex$L
        P <- ex$P
        t7 <- system.time(y <- solve(t(L),z))
        t8 <- system.time(y <- y+mu)
        t9 <- system.time(L %*% z)

        output <- rbind(output,
                        data_frame(rep=rep,N=N,
                                   logpost=t1["elapsed"],
                                   gradient=t2["elapsed"],
                                   hessian=t3["elapsed"],
                                   mode=t4["elapsed"],
                                   rnorm=t5["elapsed"],
                                   chol=t6["elapsed"],
                                   solve=t7["elapsed"],
                                   add=t8["elapsed"],
                                   mult=t9["elapsed"],
                                   tot_mvn=rnorm+chol+solve+add+mult
                                   )
                        )
    }
    return(output)
}


res <- ldply(N.vec, run_times, .parallel=run.par)
D <- gather(res, stat, time, logpost:tot_mvn) %>%
  mutate(plotID=as.factor(ifelse(stat %in% c("logpost","gradient"),
                                 "logpost/gradient",
                                 ifelse(stat == "hessian",
                                        "hessian",
                                        ifelse(stat == "tot_mvn",
                                               "total_mvn",
                                               "mvn_steps"))
                                 )
                          )
         )  %>%
  group_by(N, stat, plotID) %>%
  summarize(seconds=mean(time))


theme_set(theme_bw())
P1 <-  ggplot(filter(D, plotID=="logpost/gradient"),
              aes(x=N, y=seconds, color=stat)) %>%
  + geom_line() %>%
  + geom_point()

P2 <-  ggplot(filter(D, plotID=="hessian"),
              aes(x=N, y=seconds, color=stat)) %>%
  + geom_line() %>%
  + geom_point()

P3 <-  ggplot(filter(D, plotID=="total_mvn"),
              aes(x=N, y=seconds, color=stat)) %>%
  + geom_line() %>%
  + geom_point()

P4 <-  ggplot(filter(D, plotID=="mvn_steps"),
              aes(x=N, y=seconds, color=stat)) %>%
  + geom_line() %>%
  + geom_point()



