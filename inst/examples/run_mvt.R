library(plyr)
library(foreach)
library(doParallel)
library(trustOptim)
library(bayesGDS)
library(BD2015)

mod.name <- "mvt"


run.par <- TRUE
ncores <- 12
if(run.par) registerDoParallel(cores=ncores) else registerDoParallel(cores=1)

set.seed(12345)

n.iter <- 20
obs.list <- c(200,2000)
nx.list <- c(5,25,100)
M.list <- c(1000,10000)
scale.list <- c(0.8, 0.7, 0.6, 0.5)


res.names <- c("iter","n.x", "n.obs","k","scale","M","LML.true","LML.ds",
               "acc.rate")

res <- matrix(nrow=n.iter*length(obs.list)*length(nx.list)*length(M.list)*length(scale.list), ncol=length(res.names))
colnames(res) <- res.names

n.draws <- 24
batch.size <- ceiling(n.draws/ncores)


control.list <- list(start.trust.radius=.01,
                     report.freq=20,
                     report.level=5,
                     report.precision=5,
                     maxit=20000,
                     function.scale.factor = as.numeric(-1)
                     )


idx <- 0

for (iter in 1:n.iter) {
    for (n.obs in obs.list) {
        for (n.x in nx.list) {
            X <- matrix(rnorm(n.obs*n.x, 0,1),n.x,n.obs)
            Y <- 5 +
              matrix(seq(-5,5,length=n.x),1,n.x) %*% X[1:n.x,] + matrix(rnorm(n.obs),1,n.obs)

            cat("iter: ",iter,"\tn.obs: ",n.obs,"\tnx: ",n.x,"\t")
            x.mod <- rbind(1,X[1:n.x,])
            k <- n.x + 1
            data <- list(X=x.mod,Y=Y)

            ## hyperparameters

            b0 <- rep(0,k)
            chol.inv.V0 <- 0.2*diag(k)
            r0 <- 2
            s0 <- 1

            hyps <- list(b0=b0, chol.inv.V0=chol.inv.V0, r0=r0, s0=s0)
            params.list <- list(data=data,priors=hyps)

            startX <- rep(0,k+1)

            cat("setting up function\n")
            cl <- new("mvt", params.list)
            cl$record.tape(startX)
            cl$init.hessian(startX)
            f <- cl$get.f(startX)
            print(f)

            get.f <- function(P) return(cl$get.f(P))
            get.df <- function(P) return(cl$get.fdf(P)$grad)
            get.fdf <- function(P) return(cl$get.fdf(P))
            get.hessian <- function(P) return(cl$get.hessian(P))

            opt <- trust.optim(startX, fn=get.f, gr=get.df, ##hs = get.hessian,
                               method="BFGS",
                               control=control.list
                               )

            post.mode <- opt$solution
            gr <- opt$gradient
            B <- as(-get.hessian(post.mode),"sparseMatrix")
            LML.true <- get.marg.LL(data, hyps)

            for (M in M.list) {
                for (prop.scale in scale.list) {
                    cat("scale: ",prop.scale,"\tM: ",M,"\n")
                    idx <- idx+1

                    ## Direct Sampling

                    prop.params <- list(mean = post.mode,
                                        CH = Matrix::Cholesky(prop.scale*B))


                    log.c1 <- opt$fval
                    log.c2 <- dmvn.sparse.wrap(post.mode, prop.params)
                    log.const <- log.c1 - log.c2

                    cat("Collecting GDS Proposal Draws\n")

                    draws.m <- as(rmvn.sparse.wrap(M,prop.params),"matrix")
                    log.post.m <- aaply(draws.m, 1,get.f, .parallel=TRUE)
                    log.prop.m <- dmvn.sparse.wrap(draws.m,params=prop.params)
                    log.phi <- log.post.m - log.prop.m +log.c2 - log.c1

                    cat("Are any log.phi > 0?  ",any(log.phi>0),"\n")
                    if (all(log.phi<=0)) {

                        n.batch <- floor(n.draws / batch.size)
                        rm(draws.m)
                        gc()
                        cat("Generating DS draws - accept-reject phase\n")
                        rs <- .Random.seed
                        draws.list <- foreach(i=1:n.batch, .inorder=FALSE) %dopar% sample.GDS(
                            n.draws = batch.size,
                            log.phi=log.phi,
                            post.mode=post.mode,
                            fn.dens.post = get.f,
                            fn.dens.prop = dmvn.sparse.wrap,
                            fn.draw.prop = rmvn.sparse.wrap,
                            prop.params = prop.params,
                            max.tries=5000000,
                            announce=TRUE,
                            report.freq=1000,
                            thread.id=i,
                            seed=rs+iter+i)

                        draws.sum <- Reduce(function(x,y) Map(rbind,x,y), draws.list)

                        counts <- draws.sum$counts
                        gt.1 <- draws.sum$gt.1
                        draws <- t(draws.sum$draws)

                        log.prop.dens <- draws.sum$log.prop.dens


                        ## get MVT likelihood

                        LML.ds <- get.LML(counts=counts,
                                          log.phi=log.phi,
                                          post.mode=post.mode,
                                          fn.dens.post=get.f,
                                          fn.dens.prop=dmvn.sparse.wrap,
                                          prop.params=prop.params)
                        acc.rate <- 1/mean(counts)


                        res[idx,]<- c(iter, n.x, n.obs, k, prop.scale, M, LML.true, LML.ds,
                                      acc.rate)
                    } else {
                        res[idx,1:7] <- c(iter, n.x, n.obs, k, prop.scale, M, LML.true)
                    }
                }
            }
        }
    }
}




