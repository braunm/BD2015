library(trustOptim)
library(sparseMVN)
library(bayesGDS)
library(doParallel)
library(bayesm)
library(plyr)
library(reshape2)
library(dplyr)

set.seed(123)

registerDoParallel(cores=10)



data(cheese)
x <- cheese
colnames(x) <- c("Retailer","Volume","Disp","Price")
x <- transform(x,RetID=as.integer(Retailer))
x2 <- ddply(x,"RetID",function(j) subset(transform(j,week=1:length(j$Disp)),
                                         select=c("RetID","week","Volume","Disp","Price"))
            )

T <- as.integer(daply(x2,"RetID",function(j) length(j$RetID)))
N <- length(levels(as.factor(x2$RetID)))

x3 <- melt(x2,id.var=c("RetID","week"))
y <- acast(x3,RetID~week~variable)
y[is.na(y)] <- 0

volume <- y[,,"Volume"]/1000
disp <- y[,,"Disp"]
log_price <- log(y[,,"Price"])
log_price[!is.finite(log_price)] <- 0

k <- 2+1
maxT <- max(T)
nu <- k+3
mu_prior_mean <- rep(0,k)
mu_prior_cov <- 100*diag(k)
mu_prior_chol_prec <- t(chol(solve(mu_prior_cov)))
alpha_prior <- 5
A <- diag(k)
chol_A <- t(chol(A))

data <- list(volume=t(volume), log_price=t(log_price),
             disp=t(disp), T=T)

priors <- list(mu_prior_mean=mu_prior_mean,
               mu_prior_chol_prec=mu_prior_chol_prec,
               nu=nu,
               alpha_prior=alpha_prior,
               chol_A=chol_A)

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

nvars <- N*k + k + N + k*(k+1)/2

start.list <- as.relistable(list(beta=matrix(0,k,N),
                                 log_alpha= rep(0,N),
                                 mu= rep(0,k),
                                 chol.G=diag(k)
                                 )
                            )

start <- unlist(start.list)
cl <- new("conjugate", L)
cl$record.tape(start)
cl$init.hessian(start)

get.f <- function(P) return(cl$get.f(P))
get.df <- function(P) return(cl$get.fdf(P)$grad)
get.fdf <- function(P) return(cl$get.fdf(P))
get.hessian <- function(P) return(cl$get.hessian.sparse(P))


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
