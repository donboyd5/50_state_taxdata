
source(here::here("include", "libraries.r"))
# remotes::install_github("tidyverse/dplyr") if needed
devtools::session_info()

library(mvtnorm)
library(maxLik)


## estimate mean and variance of normal random vector
set.seed(123)

x <- rnorm(50, 1, 2 )
x <- rnorm(50e3, 1, 2 )

## log likelihood function.
## Note: 'param' is a vector
llf <- function( param ) {
  mu <- param[ 1 ]
  sigma <- param[ 2 ]
  llValue <- dnorm(x, mean=mu, sd=sigma, log=TRUE)
  return(sum(llValue))
}

## Estimate it. Take standard normal as start values
ml <- maxLik( llf, start = c(mu=0, sigma=1) )
print(summary(ml))

ml2 <- maxNR( llf, start = c(mu=0, sigma=1) )
print(summary(ml2))

## Estimates close to c(1,2) :-)
## Example how to use maxLik in your own function and allow users
## to override the default parameters
##
## 'estimate': user contructed estimation routine
## Note: it accepts both 'control' and '...'
estimate <- function(control=NULL, ...) {
  return(maxLik(llf, start=c(1,1),
                control=c(list(iterlim=100), control),
                # user-supplied 'control' overrides default
                # 'iterlim=100'
                ...))
}
m <- estimate(control=list(iterlim=1), fixed=2)
# user can override default 'iterlim' and
# supply additional parameters ('fixed')
show(maxControl(m))
# iterlim should be 1
print(coef(m))
# sigma should be 1.00


# Maximum Likelihood estimation of Poissonian distribution
n <- rpois(100, 3)
n <- rpois(50e3, 3)
n

loglik <- function(l) n*log(l) - l - lfactorial(n)

# we use numeric gradient
summary(maxBFGS(loglik, start=1))

# you would probably prefer mean(n) instead of that ;-)

# Note also that maxLik is better suited for Maximum Likelihood
###
### Now an example of constrained optimization
###
f <- function(theta) {
  x <- theta[1]
  y <- theta[2]
  exp(-(x^2 + y^2))
  ## you may want to use exp(- theta %*% theta) instead
}

## use constraints: x + y >= 1
A <- matrix(c(1, 1), 1, 2)
A
B <- -1
res <- maxNM(f, start=c(1,1), constraints=list(ineqA=A, ineqB=B),
             control=list(printLevel=1))
print(summary(res))



# a simple two-dimensional exponential hat
f <- function(a) exp(-a[1]^2 - a[2]^2)
#
# maximize wrt. both parameters
free <- maxNR(f, start=1:2)
summary(free) # results should be close to (0,0)
activePar(free)

# keep the first parameter constant
cons <- maxNR(f, start=1:2, fixed=c(TRUE,FALSE))
summary(cons) # result should be around (1,0)
activePar(cons)

hatf <- function(theta) {
  x <- theta[1]
  y <- theta[2]
  exp(-(x^2 + y^2))
  ## Note: you may prefer exp(- theta %*% theta) instead
}

## use constraints: x + y = 1
A <- matrix(c(1, 1), 1, 2)
A
B <- -1
res <- sumt(hatf, start=c(0,0), maxRoutine=maxNR,
            constraints=list(eqA=A, eqB=B))
print(summary(res))



### Now an example of constrained optimization
###
f <- function(theta) {
  x <- theta[1]
  y <- theta[2]
  exp(-(x^2 + y^2))
  ## you may want to use exp(- theta %*% theta) instead
}

## use constraints: x + y >= 1
A <- matrix(c(1, 1), 1, 2)
A
B <- -1
res <- maxNM(f, start=c(1,1), constraints=list(ineqA=A, ineqB=B),
             control=list(printLevel=1))
print(summary(res))


res1 <- maxNR(f, start=c(1,1), constraints=list(eqA=A, eqB=B),
             control=list(printLevel=1))
print(summary(res1))


# Maximum Likelihood estimation of Poissonian distribution
n <- rpois(100, 3)
n <- rpois(50e3, 3)
n

loglik <- function(l) n*log(l) - l - lfactorial(n)

# we use numeric gradient
summary(maxBFGS(loglik, start=1))



# MLE poisson ----
h <- 10 # number households
k <- 3 # number characteristics
s <- 5 # number states

# x is an h x k matrix of characteristics
sigma <- matrix(c(1, .8, .6,
                  .8, 1, .4,
                  .6, .4, 1), ncol=3, byrow = TRUE)
sigma
set.seed(1234)
x <- rmvnorm(n=h, mean=c(10, 20, 30), sigma)
cor(x)

# delta is an h-length vector of individual constants in the formula for state weights
delta <- 1:h / 100

# beta is an s x k matrix of coefficients
beta <- matrix(c(1, 2, 3,
                 4, 5, 6,
                 7, 8, 9,
                 10, 11, 12,
                 13, 14, 15), ncol=3, byrow = TRUE) / 100
beta

getw <- function(ih, is, beta, delta, x){
  # get weight for household ih, state is
  exp(t(beta[is, ]) %*% x[ih, ] + delta[ih]) %>% as.numeric
}
getw(ih=2, is=3, beta, delta, x)

get_weights <- function(beta, delta, x){
  # get all weights
  beta_x <- beta %*% t(x)
  # add delta to every row of beta_x and transpose
  beta_xd <- apply(beta_x, 1 , function(m) m + delta) 
  exp(beta_xd)
}
get_weights(beta, delta, x)

get_delta <- function(w, beta, x){
  beta_x <- exp(beta %*% t(x))
  log(w / colSums(beta_x))
}

ih <- 10; is <- 5
getw(ih, is, beta, delta, x)
get_weights(beta, delta, x)[ih, is]


# now choose weights one of 2 ways -- either from a poisson model, or randomly
sweights <- get_weights(beta, delta, x) # poisson model
weights <- rowSums(sweights) # weights for each individual
states <- colSums(sweights)

targets <- t(sweights) %*% x
targets

beta
delta
beta %*% t(x)
theta <- c(as.vector(beta), delta)
theta
.theta <- theta

fw <- function(.theta, s, k, x){
  # get weights, given theta, s, k, x
  # theta is beta (5 x 3 ie s x k), delta (1 x h)
  .beta <- matrix(.theta[1:(s * k)], ncol=k, byrow = FALSE)
  .delta <- .theta[(s * k + 1) : length(.theta)]
  beta_x <- .beta %*% t(x)
  # add delta to every row of beta_x and transpose
  beta_xd <- apply(beta_x, 1 , function(m) m + .delta) 
  exp(beta_xd)
}
fw(theta, s, k, x)


hatf <- function(theta) {
  x <- theta[1]
  y <- theta[2]
  exp(-(x^2 + y^2))
  ## Note: you may prefer exp(- theta %*% theta) instead
}

## use constraints: x + y = 1
A <- matrix(c(1, 1), 1, 2)
A
B <- -1
res <- sumt(hatf, start=c(0,0), maxRoutine=maxNR,
            constraints=list(eqA=A, eqB=B))
print(summary(res))




ewhs <- get_weights(ebeta, edelta, x)


dtheta <- function(.theta, .s, .k, .x, .targets){
  # get weights, given theta, s, k, x
  # theta is beta (5 x 3 ie s x k), delta (1 x h)
  beta <- matrix(.theta[1:(.s * .k)], ncol=.k, byrow = FALSE)
  delta <- .theta[(.s * .k + 1) : length(.theta)]
  esw <- get_weights(beta, delta, .x)
  etargets <- t(esw) %*% .x
  dist <- .targets - etargets
  sse <- sum(dist^2) # a measure of error
  -sse
}


t(sweights) %*% x

.theta <- theta
theta <- c(as.vector(ibeta), idelta)
dtheta(theta, s, k, x, targets)

res <- maxNR(fn=dtheta, start=theta, .s=s, .k=k, .x=x, .targets=targets)
print(summary(res))

