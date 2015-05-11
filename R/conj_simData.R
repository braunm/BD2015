rm(list=ls())
library(plyr)

N <- 2000
T <- 10

mod.name <- "conjugate"
file.path <- "~/Documents/gds-mksc/conjugate/"
save.file <- paste0(file.path,"data/",mod.name,".Rdata")

## true values for simulated data

tau <- 3
sigma <- 2
mu <- -1
theta <- rnorm(N,mu,tau)
true.values <- list(theta=theta, mu=mu, sigma=sigma, tau=tau)

Y <- laply(theta, function (a) return(rnorm(T,a,sigma)))

save(Y, true.values, file=save.file)
