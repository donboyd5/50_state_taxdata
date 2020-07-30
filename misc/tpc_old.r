
source(here::here("include", "libraries.r"))
library(mvtnorm)

h <- 10 # number households
k <- 3 # number characteristics
s <- 5 # number states

# try a large problem
# step_scale <- 1e4 and xscale <- 1e4 works well with this
h <- 40e3 # number households
k <- 30 # number characteristics
s <- 50 # number states

# x is an h x k matrix of characteristics
sigma <- matrix(c(1, .8, .6,
                  .8, 1, .4,
                  .6, .4, 1), ncol=3, byrow = TRUE)
sigma
set.seed(1234)
x <- rmvnorm(n=h, mean=c(10, 20, 30), sigma)

# alternative way to get x
xbar <- seq(10, 100, length.out=k)
xsd <- xbar / 10
set.seed(1234)
x <- matrix(rnorm(n=h * k, mean=xbar, sd=xsd), nrow=h, ncol=k)
x

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

ih <- 10; is <- 5
getw(ih, is, beta, delta, x)
get_weights(beta, delta, x)[ih, is]


get_delta <- function(w, beta, x){
  beta_x <- exp(beta %*% t(x))
  log(w / colSums(beta_x))
}


# now choose weights one of 2 ways -- either from a poisson model, or randomly
sweights <- get_weights(beta, delta, x) # poisson model

# randomly
set.seed(1234)
sweights <- matrix(runif(n=h * s, min=4, max=8000), nrow=h, ncol=s)
sweights

# proceed with whatever set of weights we have
weights <- rowSums(sweights) # weights for each individual
weights
states <- colSums(sweights)
states

# create data
colnames(x) <- paste0("x", 1:ncol(x))
colnames(sweights) <- paste0("s", 1:ncol(sweights))
data <- as_tibble(cbind(pid=1:nrow(x), x, sweights))
data
  
# calculate targets - sum of x1-x3 in each state
targets_df <- data %>%
  pivot_longer(cols=starts_with("s"), names_to="state", values_to="weight") %>%
  group_by(state) %>%
  summarise(across(x1:x3, ~ sum(. * weight)))
targets_df

xsave <- x


xscale <- 1e4
x <- xsave / xscale

targets <- t(sweights) %*% x
targets

# now search for the beta coefficients and deltas that are consistent with these data and that meet the targets
# we know the data and we know the aggregate weights, but not the state weights -- e prefix will mean estimated

xpx <- t(x) %*% x
# inverse of matrix that has been multiplied by a non-zero scalar equals inverse of the scalar multiplied by inverse of the matrix,
# so solve xpx once at the start
(ixpx <- solve(xpx)) # this needs to be invertible
# to get the inverse of the jacobian we will multiply by the inverse of the scalar

iebeta <- matrix(0, nrow=s, ncol=k) # tpc uses 0 as a starting point
iedelta <- get_delta(weights, iebeta, x)

iweights <- get_weights(iebeta, iedelta, x)
idist <- targets - t(iweights) %*% x
isse <- sum(idist^2)
iweights
isse
 
# before start, set initial values
ebeta <- iebeta
edelta <- iedelta
step_scale <- 1e4 # djb note that this has a MAJOR impact on iterations

a <- proc.time()
for(i in 1:100){
  ewhs <- get_weights(ebeta, edelta, x) # estimated weights for each household for each state
  ew <- rowSums(ewhs) # estimated total weights for each individual
  es <- colSums(ewhs) # estimated total weights for each state
  
  etargets <- t(ewhs) %*% x
  dist <- targets - etargets
  sse <- sum(dist^2) # a measure of error
  
  if(sse < 1e-8) {
    print(sprintf("DONE at iteration %i: sse is %.5e", i, sse))
    break
  }
  
  if(i<=10 | i %% 5 ==0) {
    print(sprintf("iteration %i: sse is %.5e", i, sse))
  }
  
  # get the step
  step <- matrix(nrow=s, ncol=k)
  for(i in 1:s) step[i, ] <- t((1 /es[i]) * ixpx %*% dist[i, ]) * step_scale
  
  ebeta <- ebeta + step
  edelta <- get_delta(ew, ebeta, x)
}
b <- proc.time()
b - a

sum(dist^2)
dist
quantile(dist)

fbeta <- ebeta
fdelta <- edelta
fweights <- get_weights(fbeta, fdelta, x)

fbeta %>% round(3)
beta %>% round(3)
iebeta %>% round(3)

fdelta %>% round(3)
delta %>% round(3)
iedelta %>% round(3)

fweights %>% round(3)
sweights %>% round(3)
iweights %>% round(3)
(fweights - sweights) %>% round(6)
fweights / sweights
fweights / iweights

rowSums(fweights)
weights
(rowSums(fweights) / weights * 100 - 100) %>% round(3)

colSums(fweights)
states
(colSums(fweights) / states * 100 - 100) %>% round(3)


t1 <- targets * xscale
t2 <- (t(fweights) %*% x) * xscale
t1 %>% round(3)
t2 %>% round(3)
(t2 - t1) %>% round(3)
(t2 / t1 * 100 - 100) %>% round(4)
