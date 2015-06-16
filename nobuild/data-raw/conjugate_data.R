set.seed(1234)
N <- 1500
T <- 10

## true values for simulated data

tau <- 3
sigma <- 2
mu <- -1
theta <- rnorm(N,mu,tau)
conj_true <- list(theta=theta, mu=mu, sigma=sigma, tau=tau)
conj_Y <- plyr::laply(theta, function (a) return(rnorm(T,a,sigma)))

devtools::use_data(conj_Y, conj_true, overwrite=TRUE)

