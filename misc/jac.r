 
# jacobian
# t is a set of targets
h <- 4 # number households
k <- 2 # number characteristics
s <- 3 # number states

msize <- k * s
# vlen <- k + s


# functions ----
vtm <- function(vec, s){
  matrix(vec, nrow=s, byrow=FALSE)
}

get_delta <- function(wh, beta, x){
  beta_x <- exp(beta %*% t(x))
  log(wh / colSums(beta_x))
}

get_weights <- function(beta, delta, x){
  # get all weights
  beta_x <- beta %*% t(x)
  # add delta to every row of beta_x and transpose
  beta_xd <- apply(beta_x, 1 , function(m) m + delta) 
  # beta_xd <- t(beta_x) # uncomment to see the impact of not updating delta
  exp(beta_xd)
}

# analysis ----
xmat <- matrix(c(5, 6,
              7, 4,
              9, 2,
              3, 1), 
            nrow=h, byrow=TRUE)
xxp <- xmat %*% t(xmat)
xxp

betavec <- seq(.1, .9, length.out=msize)
beta <- vtm(betavec, s)
beta

delta <- 1:h

whs <- get_weights(beta, delta, xmat)

(wh=rowSums(whs))
(ws=colSums(whs))

(ctargets <- t(whs) %*% xmat) # s x k targets in the data

set.seed(1234)
tadj <- runif(msize, .9, 1.1) # calculated targets
targets <- ctargets * tadj # targets I create

d <- targets - ctargets
d


f <- function(betavec, wh, xmat, targets){
  # return the distances from targets as a vector
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  as.vector(d)
}

f(betavec, wh, xmat, targets)
d
diff <- 1e-8
v <- c(diff, 0, 0, 0, 0, 0)
v <- c(0, diff, 0, 0, 0, 0)
v <- c(0, 0, 0, 0, 0, diff)
(f(betavec + v, wh, xmat, targets) - f(betavec, wh, xmat, targets)) / diff


# finite-differences jacobian
jfd <- jacobian(f, x=as.vector(beta), method="Richardson", wh=wh, xmat=xmat, targets=targets)
jfd

# calculated jacobian
x <- xmat
sum(whs[, 1] * x[, 1] * x[, 1]) # b11
sum(whs[, 2] * x[, 1] * x[, 1]) # b21
sum(whs[, 3] * x[, 1] * x[, 1]) # b31

sum(whs[, 1] * x[, 2] * x[, 2]) # b11
sum(whs[, 2] * x[, 2] * x[, 2]) # b21
sum(whs[, 3] * x[, 2] * x[, 2]) # b31

# let's figure row 1 of the matrix: change in d[j] wrt b11

# for each person h we have, for state 1
# s1k1
w_hs1 <- expression(exp(b11*xh1 + b12*xh2 + ch) * xh1)
D(w_hs1, "b11") # w_hs1 * xh1 * xh1
# for element 1 of the row (s1k1), sum these over all h
sum(whs[, 1] * x[, 1] * x[, 1]) # b11

# element 2 of the row
# s2k1
w_hs2 <- expression(exp(b21*xh1 + b22*xh2 + ch) * xh1)
D(w_hs2, "b11") # w_hs1 * xh1 * xh1



D(w_h1s2, "b22") # w_h1s2 * x12
exp(b21 * x11 + b22 * x12 + c1) * x12


