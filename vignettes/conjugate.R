library(trustOptim)
library(sparseMVN)
library(bayesGDS)
library(doParallel)
library(ggplot2)
library(gridExtra)

set.seed(1234)

N <- 1500
T <- 10

tau <- 3
sigma <- 2
mu <- -1
theta <- rnorm(N,mu,tau)

Y <- laply(theta, function (a) return(rnorm(T,a,sigma)))

data <- list(Y=Y)
priors <- NULL

L <- list(data=data, priors=priors)

nvars <- N + 3
start.list <- as.relistable(list(theta=rnorm(N, 0, 3),
                                 mu=0,
                                 log_sig=0,
                                 log_tau=1.1)
                            )

start <- unlist(start.list)
cl <- new("conjugate", L)
cl$record.tape(start)
cl$init.hessian(start)

get.f <- function(P) return(cl$get.f(P))
get.df <- function(P) return(cl$get.df(P))
get.fdf <- function(P) return(cl$get.fdf(P))
get.hessian <- function(P) return(cl$get.hessian.sparse(P))

## Posterior mode
opt <- trust.optim(start,
                   fn=get.f,
                   gr = get.df,
                   hs = get.hessian,
                   method = "Sparse",
                   control = list(
                       report.freq=10L,
                       report.level=3L,
                       report.precision=1L,
                       maxit=500L,
                       function.scale.factor = as.numeric(-1),
                       preconditioner=1L
                       )
                   )

post.mode <- opt$solution
sol <- relist(post.mode, skeleton=start.list)

## Rejection Sampling

run.par <- TRUE
if(run.par) registerDoParallel(cores=12) else registerDoParallel(cores=1)

n.draws <- 3
batch.size <- 1

sc <- .94
M <- 20000

log.c1 <- get.f(post.mode)

chol.hess <- Cholesky(-sc*opt$hessian)
prop.params <- list(mean = post.mode,
                    CH = chol.hess
                    )

log.c2 <- dmvn.sparse.wrap(matrix(post.mode,nrow=1), prop.params)
log.const <- log.c1 - log.c2

cat("Collecting GDS Proposal Draws\n")

draws.m <- as(rmvn.sparse.wrap(M,prop.params),"matrix")
log.post.m <- plyr::aaply(draws.m, 1,get.f, .parallel=TRUE)
log.prop.m <- dmvn.sparse.wrap(draws.m,params=prop.params)
log.phi <- log.post.m - log.prop.m +log.c2 - log.c1

stopifnot(all(log.phi <= 0))

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
    max.tries=5000000,
    report.freq=1000,
    announce=TRUE,
    thread.id=i,
    seed=i)

dd <- Reduce(function(x,y) Map(rbind,x,y), draws)
dd$LML <- get.LML(counts=dd$counts,log.phi=log.phi,
                  post.mode=post.mode,
                  fn.dens.post=get.f,
                  fn.dens.prop=dmvn.sparse.wrap,
                  prop.params=prop.params)
dd$acc.rate <- 1/mean(dd$counts)

## Gibbs sampler

ngibbs <- 2000
report <- 500
warmup <- 1000

gibbs <- matrix(0,ngibbs-warmup,nvars)
LL <- vector("numeric",length=ngibbs-warmup)

theta <- rep(0,N)
mu <- 1
sig2 <- .1
tau2 <- .1

sig.df <- T*N

for (d in 1:ngibbs) {

    if (d %% report ==0) {
        cat("Iteration: ",d,"\n")
    }

    theta.prec <- 1/tau2 + T/sig2
    theta.sd <- 1/sqrt(theta.prec)
    theta.mean <- (mu/tau2 + T*rowMeans(Y)/sig2)/theta.prec
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
        gibbs[d-warmup,] <- res
        LL[d-warmup] <- cl$get.f(res)
    }
}

R <- NROW(gibbs)
var.names <- c(paste0("theta.",1:N),"mu","log.sig","log.tau")
nvars <- length(var.names)

colnames(gibbs) <- var.names
draws.gibbs <- as.data.frame(gibbs) %>%
  mutate(iterations=1:R,
         tau2=exp(2*log.tau),
         sig2=exp(2*log.sig)) %>%
  select(iterations, mu, tau2, sig2)  %>%
  melt(id.vars="iterations", variable.name="parameters") %>%
  mutate(method="Gibbs")

logpost.gibbs <- data.frame(iterations=1:R,parameters="log.post",
                            value=LL, method="Gibbs")

MM <- NROW(dd$draws)
colnames(dd$draws) <- var.names
draws.gds <- as.data.frame(dd$draws) %>%
  mutate(iterations=1:MM,
         tau2=exp(2*log.tau),
         sig2=exp(2*log.sig)) %>%
  select(iterations, mu, tau2, sig2) %>%
  melt(id.vars="iterations", variable.name="parameters") %>%
  mutate(method="Ours")

names(dd$log.post.dens) <- 1:MM
logpost.gds <- data.frame(iterations=1:MM,
                          parameters="log.post",
                          value=dd$log.post.dens,
                          method="Ours")

d <- rbind(draws.gibbs, draws.gds, logpost.gds, logpost.gibbs)
d$parameters <- reorder(d$parameters,new.order=c("mu","tau2","sig2","log.post"))

B <- ggplot(data=d,aes(x=method,y=value)) + geom_boxplot(outlier.size=.6,size=.2)
B <- B + theme(axis.text=element_text(size=8),
               axis.title=element_blank(),
               strip.background=element_rect(fill="white"))
B <- B + facet_wrap(~parameters, scales="free",nrow=1,ncol=4)

save(d, file="conj_gibbs.Rdata")
print(B)
