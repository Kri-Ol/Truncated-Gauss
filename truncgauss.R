require(rootSolve)
require(msm)

phi <- function(z) {
    dnorm(z)
}

Phi <- function(z) {
    pnorm(z)
}

Mean <- function(mu, sigma, a, b) {
    alfa <-  (a - mu) / sigma
    beta <-  (b - mu) / sigma

    Z <-  Phi(beta) - Phi(alfa)

    mu + sigma*(phi(alfa) - phi(beta))/Z
}

f <- function(mu, mean, sigma, a, b) {
    mean - Mean(mu, sigma, a, b)
}

a <-  50000.0
b <-  250000.0
mean  <- 70000.0
sigma <- 24000.0

q <- uniroot(f, c(a, b), mean, sigma, a, b)
mu <- q$root

print(sprintf("Found mu = %f for the desired mean %f and sigma %f", mu, mean, sigma))

set.seed(32345)
N = 100000
r <-rtnorm(N, mean=mu, sd=sigma, lower=a, upper=b)

print(sprintf("Sampled %d truncated gaussians and got observed mean = %f", N, mean(r)))
