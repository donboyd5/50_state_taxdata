
# This program is based largely on the methodology in:
#   Khitatrakun, Surachai, Gordon B T Mermin, and Norton Francis. “Incorporating State Analysis into the 
#   Tax Policy Center’s Microsimulation Model: Documentation and Methodology.” Working Paper, March 2016.
#   https://www.taxpolicycenter.org/sites/default/files/alfresco/publication-pdfs/2000697-Incorporating-State-Analysis-into-the-TPCs-Microsimulation-Model.pdf.


source(here::here("include", "libraries.r"))
library(mvtnorm)

devtools::session_info()
(.packages()) %>% sort


# define needed functions ----
get_delta <- function(w, beta, x){
  beta_x <- exp(beta %*% t(x))
  log(w / colSums(beta_x))
}

get_weights <- function(beta, delta, x){
  # get all weights
  beta_x <- beta %*% t(x)
  # add delta to every row of beta_x and transpose
  beta_xd <- apply(beta_x, 1 , function(m) m + delta) 
  exp(beta_xd)
}


# create the poisson problem ----
# xscale <- 10 and step_scale <- 10 works well with this
h <- 10 # number households
k <- 3 # number characteristics
s <- 5 # number states

# for the poisson version we need beta and delta
# beta is an s x k matrix of coefficients
beta <- matrix(c(1, 2, 3,
                 4, 5, 6,
                 7, 8, 9,
                 10, 11, 12,
                 13, 14, 15), ncol=3, byrow = TRUE) / 100
beta

# delta is an h-length vector of individual constants in the formula for state weights
delta <- 1:h / 100


# create multivariate normal x for the 3-characteristics case ----
# x is an h x k matrix of characteristics
sigma <- matrix(c(1, .8, .6,
                  .8, 1, .4,
                  .6, .4, 1), ncol=3, byrow = TRUE)
sigma

set.seed(1234)
x <- rmvnorm(n=h, mean=c(10, 20, 30), sigma)

# get household-state weights if poisson model ----
hsweights <- get_weights(beta, delta, x) # poisson model
# end poisson ----

# Alternatively, create a large problem ----
# xscale <- 100 and step_scale <- nrow(x) works well with this
h <- 40e3 # number households
k <- 30 # number of characteristics
s <- 50 # number of states


#.. get household-state weights randomly IF not poisson ----
set.seed(1234)
hsweights <- matrix(runif(n=h * s, min=4, max=8000), nrow=h, ncol=s)
#.. end random ----


# alternative way to get x for large or small problem ----
xbar <- seq(10, 100, length.out=k)
xsd <- xbar / 10

set.seed(1234)
x <- matrix(rnorm(n=h * k, mean=xbar, sd=xsd), nrow=h, ncol=k)
# end alternative way to get x ----


# look at x and its correlations ----
x
cor(x)

# get total household weights and state weights ----
hsweights
hweights <- rowSums(hsweights) # sum of weights for each individual
hweights
sweights <- colSums(hsweights) # sum of weights for each state
sweights

# done with data setup ----


# get problem ready to solve ----
xsave <- x # save initial x values in case we scale them


# scale and go ----
# need to get a rule of thumb for this
# xscale <- 1
(xscale <- sum(xsave) / 1000)

x <- xsave / xscale

targets <- t(hsweights) %*% x
targets

# now search for the beta coefficients and deltas that are consistent with these data and that meet the targets
# we know the data and we know the aggregate weights, but not the state weights -- e prefix will mean estimated

# define xprime-x and its inverse before entering the loop as it will not change within the loop
xpx <- t(x) %*% x
# inverse of matrix that has been multiplied by a non-zero scalar equals inverse of the scalar multiplied by inverse of the matrix,
# so solve xpx once at the start
(invxpx <- solve(xpx)) # this needs to be invertible
# to get the inverse of the jacobian we will multiply by the inverse of the scalar


# define initial values before entering loop ----
ibeta <- matrix(0, nrow=s, ncol=k) # tpc uses 0 as beta starting point
idelta <- get_delta(hweights, ibeta, x) # tpc uses initial delta based on initial beta 

ihsweights <- get_weights(ibeta, idelta, x)
ihsweights

# look at distance from targets using initial values
idist <- targets - t(ihsweights) %*% x
idist

isse <- sum(idist^2)
isse
 
# before start, set values to use when entering loop for first time ----
ebeta <- ibeta
edelta <- idelta

# note that step_scale has a MAJOR impact on iterations
(step_scale <- nrow(x))

a <- proc.time()
for(i in 1:500){
  ehsweights <- get_weights(ebeta, edelta, x) # estimated weights for each household for each state
  ehweights <- rowSums(ehsweights) # estimated total weights for each individual
  esweights <- colSums(ehsweights) # estimated total weights for each state
  
  etargets <- t(ehsweights) %*% x
  dist <- targets - etargets
  sse <- sum(dist^2) # sum of squared errors
  sseu <- sse * (xscale ^2) # what the unscaled sse would be
  
  if(sseu < 1e-6 | (sse < 1e-10 & sseu < 1e-2)) {
    print(sprintf("DONE at iteration %i: scaled sse: %.5e, unscaled sse: %.5e", i, sse, sseu))
    break
  }
  
  if(i<=10 | i %% 5 ==0) {
    print(sprintf("iteration %i: scaled sse: %.5e, unscaled sse: %.5e", i, sse, sseu))
  }
  
  # get the step
  # step <- matrix(nrow=s, ncol=k)
  # for(i in 1:s) step[i, ] <- t((1 /esweights[i]) * invxpx %*% dist[i, ]) * step_scale
  step <- (1 / esweights) * dist %*% invxpx * step_scale
  
  ebeta <- ebeta + step
  edelta <- get_delta(ehweights, ebeta, x)
}
b <- proc.time()
b - a


sse; sseu
dist
quantile(dist)

# calculate final weights 
fbeta <- ebeta
fdelta <- edelta
fhsweights <- get_weights(fbeta, fdelta, x)
(fhweights <- rowSums(fhsweights))
(fsweights <- colSums(fhsweights))


fbeta %>% round(3)
beta %>% round(3)
ibeta %>% round(3)

fdelta %>% round(3)
delta %>% round(3)
idelta %>% round(3)

fhsweights %>% round(3)
hsweights %>% round(3)
ihsweights %>% round(3)
(fhsweights - hsweights) %>% round(6)
fhsweights / hsweights
fhsweights / ihsweights

quantile((fhsweights - hsweights) / hsweights * 100 - 100)

fhweights
hweights
(fhweights / hweights * 100 - 100) %>% round(3)

fsweights
sweights
(fsweights / sweights * 100 - 100) %>% round(3)


t1 <- targets * xscale
t2 <- (t(fhsweights) %*% x) * xscale
t1 %>% round(3)
t2 %>% round(3)
(t2 - t1) %>% round(3)
(t2 / t1 * 100 - 100) %>% round(4)


# regression of initial household weights hweights on x variables
mod <- glm(round(hweights) ~ x, family="poisson")
summary(mod)
summary(glm(round(hsweights[, 5]) ~ x, family="poisson")) # 4, 5 are poisson

# poisson regressions with final estimated weights
mod2 <- glm(round(fhweights) ~ x, family="poisson")
summary(mod2)

# now pick a state and estimate that model
fhweights_s <- fhsweights[, 3]
mod3 <- glm(round(fhweights_s) ~ x, family="poisson")
summary(mod3)



