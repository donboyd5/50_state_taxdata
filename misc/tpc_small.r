
# This program is based largely on the methodology in:
#   Khitatrakun, Surachai, Gordon B T Mermin, and Norton Francis. “Incorporating State Analysis into the 
#   Tax Policy Center’s Microsimulation Model: Documentation and Methodology.” Working Paper, March 2016.
#   https://www.taxpolicycenter.org/sites/default/files/alfresco/publication-pdfs/2000697-Incorporating-State-Analysis-into-the-TPCs-Microsimulation-Model.pdf.


source(here::here("include", "libraries.r"))
library(mvtnorm)

library(numDeriv)
grad(sin, pi)
grad(get_weights, beta, delta, x)

fvec <- function(betavec, wh, xmat, targets){
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  sse <- sum(d^2)
  sse
}

fvec(as.vector(beta2), wh, xmat, targets)

betavec <- as.vector(beta2)
xmat=x
fvec(as.vector(beta2), wh, xmat=xmat, targets=targs)
sum(d2^2)
grad(fvec, x=as.vector(beta2), wh=wh, xmat=xmat, targets=targs)

f1(wh, beta2, x, targs)

devtools::session_info()
(.packages()) %>% sort


# function to get step
stepc <- function(whs, d, beta, xmat){
  xmat <- x
  wh <- rowSums(whs)
  betavec <- as.vector(beta)
  
  g <- grad(fvec, x=betavec, wh=wh, xmat=xmat, targets=targs)
  h <- hessian(fvec, x=betavec, wh=wh, xmat=xmat, targets=targs)
  st <- solve(h * fvec(betavec, wh=wh, xmat=xmat, targets=targs)) %*% g
  #st <- solve(h * sum(d^2)) %*% g
  stm <- matrix(st, nrow=nrow(targets), byrow=FALSE)
  -stm
}
stepc(whs2, d2, beta2, x)

stepd <- function(whs, d, beta, xmat){
  xmat <- x
  wh <- rowSums(whs)
  betavec <- as.vector(beta)
  
  g <- grad(fvec, x=betavec, wh=wh, xmat=xmat, targets=targs)
  h <- hessian(fvec, x=betavec, wh=wh, xmat=xmat, targets=targs)
  st <- solve(h) %*% g
  #st <- solve(h * sum(d^2)) %*% g
  stm <- matrix(st, nrow=nrow(targets), byrow=FALSE)
  -stm
}
stepd(whs2, d2, beta2, x)


# define needed functions ----
get_delta <- function(wh, beta, x){
  beta_x <- exp(beta %*% t(x))
  log(wh / colSums(beta_x))
}

get_weights <- function(beta, delta, x){
  # get all weights
  beta_x <- beta %*% t(x)
  # add delta to every row of beta_x and transpose
  beta_xd <- apply(beta_x, 1 , function(m) m + delta) 
  exp(beta_xd)
}

z <- c(5, 7, 9, 3)
z * t(z)
z %*% t(z)

# create the poisson problem 1 ----
# xscale <- 10 and step_scale <- 10 works well with this
h <- 4 # number households
k <- 1 # number characteristics
s <- 3 # number states

x <- matrix(c(5,
              7,
              9,
              3), nrow=h, byrow=TRUE)

whs <- matrix(c(4, 3, 4,
                7, 12, 1,
                20, 6, 4,
                6, 4, 5), nrow=4, byrow=TRUE)

# end prob 1 ----

# alternative create poisson problem 2 ----
# xscale <- 10 and step_scale <- 10 works well with this
h <- 4 # number households
k <- 2 # number characteristics
s <- 3 # number states

x <- matrix(c(5, 6,
              7, 4,
              9, 2,
              3, 1), nrow=h, byrow=TRUE)

whs <- matrix(c(4, 3, 4,
                7, 12, 1,
                20, 6, 4,
                6, 4, 5), nrow=4, byrow=TRUE)

wh <- rowSums(whs)

# end prob 2 -----

stepa <- function(whs, d, beta, x){
  dmat <- rowSums(t(whs) * as.numeric(xpx))
  dmat
  step <- 1 / dmat * d
  step
}

stepb <- function(whs, d, beta, x){
  ws <- colSums(whs)
  step <- 1 / ws * d %*% solve(xpx)
  step
}

stepc(whs2, d2, beta2, x)


# functions for my derivative calculation
f <- function(wh, beta, x, targets){
  # this is what we want to minimize -- the sum of squared errors
  delta <- get_delta(wh, beta, x)
  whs <- get_weights(beta, delta, x)
  etargets <- t(whs) %*% x
  d <- targets - etargets
  sum(d^2)
}
f(wh, beta1, x, targs)
f(wh, beta5, x, targs)
f(wh, beta9, x, targs)
f(wh, ebeta, x, targs)

f1 <- function(wh, beta, x, targets){
  # first derivative of sum of squared errors wrt beta
  # will return 1 value per beta (will have s x k values)
  # 2x whs d simplified
  delta <- get_delta(wh, beta, x) # h values
  whs <- get_weights(beta, delta, x) # h x s values
  etargets <- t(whs) %*% x  # s x k values
  d <- targets - etargets # s x k values
  f1vec <- 2 * as.vector(etargets) * as.vector(d)
  matrix(f1vec, nrow=nrow(beta), byrow=TRUE)
}
beta <- beta1
targets <- targs
f1(wh, beta1, x, targs)
f1(wh, beta9, x, targs)
f1(wh, ebeta, x, targs)


f1 <- function(wh, beta, x, targets){
  # first derivative of sum of squared errors wrt beta
  # will return 1 value per beta (will have s x k values)
  # f1 -2x^2(e^(xbeta+d))*(t - xe^(xbeta+d))
  delta <- get_delta(wh, beta, x) # h values
  whs <- get_weights(beta, delta, x) # h x s values
  etargets <- t(whs) %*% x  # s x k values
  d <- targets - etargets # s x k values
  f1.1 <- -2*(x[, 1]^2) %*% (whs[, 1]) %*% (targs[1, 1] - d[1, 1])# beta[1] = state1, chars[1]
  f1.1 <- -2*(x[, 1] %*% x[, 1]) %*% colSums(whs)[ 1] %*% (targs[1, 1] - d[1, 1])# beta[1] = state1, chars[1]
  f1.1 <- -2*(x[, 1] * x[, 1]) %*% t(whs[, 1]) %*% (targs[1, 1] - d[1, 1])# beta[1] = state1, chars[1]
  f1.1
  
  f1vec <- 2 * as.vector(etargets) * as.vector(d)
  matrix(f1vec, nrow=nrow(beta), byrow=TRUE)
}

beta <- beta2

f(wh, beta2, x, targs)

g <- grad(fvec, x=as.vector(beta2), wh=wh, xmat=xmat, targets=targs)
h <- hessian(fvec, x=as.vector(beta2), wh=wh, xmat=xmat, targets=targs)
st <- solve(h * f(wh, beta2, x, targs)) %*% g
st


f1(wh, beta2, x, targs)

f2 <- function(wh, beta, x, targets){
  # 2nd derivative of sum of squared errors wrt beta
  # will return 1 value per beta (will have s x k values)
  # 2 * x^2 * whs * (etargs + d) simplified
  delta <- get_delta(wh, beta, x) # h values
  whs <- get_weights(beta, delta, x) # h x s values
  etargets <- t(whs) %*% x  # s x k values
  d <- targets - etargets # s x k values
  xxp <- x %*% t(x)
  xpx <- t(x) %*% x
  f2vec <- 
  f2vec <- 1
  
}

# start ----
h; k; s
(xxp <- x %*% t(x))
(xpx <- t(x) %*% x)
x; xxp; xpx
whs

(wh=rowSums(whs))
(ws=colSums(whs))

(targs <- t(whs) %*% x) # s x k

(t(whs) %*% x) %*% t(x)
t(whs) %*% (x %*% t(x))


stepfn <- stepa
stepfn <- stepb
stepfn <- stepc


# iteration #1 ----
beta1 <- matrix(0, nrow=s, ncol=k) # tpc uses 0 as beta starting point
delta1 <- get_delta(wh, beta1, x) # tpc uses initial delta based on initial beta 
beta1; delta1

whs1 <- get_weights(beta1, delta1, x)
whs1

(ws1=colSums(whs1))
(wh1=rowSums(whs1))

(targs1 <- t(whs1) %*% x)
d1 <- targs - targs1
d1

step1 <- stepfn(whs1, d1, beta1, x)
step1
stepb(whs1, d1, beta1, x)


# step1 <- 1 / ws1 * d1 %*% solve(xpx) # ???


# iter2 ----
beta2 <- beta1 + step1
delta2 <- get_delta(wh1, beta2, x) 
whs2 <- get_weights(beta2, delta2, x)
whs2

(ws2=colSums(whs2))
(wh2=rowSums(whs2))

d2 <- targs - t(whs2) %*% x
d2

step2 <- stepfn(whs2, d2, beta2, x)
# step2 <- 1 / ws2 * d2 %*% solve(xpx) # ???

# iter3
beta3 <- beta2 + step2
delta3 <- get_delta(wh2, beta3, x) 
whs3 <- get_weights(beta3, delta3, x)
whs3

(ws3=colSums(whs3))
(wh3=rowSums(whs3))

d3 <- targs - t(whs3) %*% x
d3

step3 <- stepfn(whs3, d3, beta3, x)
step3

# iter4
beta4 <- beta3 + step3
delta4 <- get_delta(wh3, beta4, x) 
whs4 <- get_weights(beta4, delta4, x)
whs4

(ws4=colSums(whs4))
(wh4=rowSums(whs4))

d4 <- targs - t(whs4) %*% x
d4

step4 <- stepfn(whs4, d4, beta4, x)
# step4 <- 1 / ws4 * d4 %*% solve(xpx) # ???
step4

# iter5
beta5 <- beta4 + step4
delta5 <- get_delta(wh4, beta5, x) 
whs5 <- get_weights(beta5, delta5, x)
whs5

(ws5=colSums(whs5))
(wh5=rowSums(whs5))

d5 <- targs - t(whs5) %*% x
d5

step5 <- stepfn(whs5, d5, beta5, x)
step5

# iter6
beta6 <- beta5 + step5
delta6 <- get_delta(wh5, beta6, x) 
whs6 <- get_weights(beta6, delta6, x)
whs6

(ws6=colSums(whs6))
(wh6=rowSums(whs6))

d6 <- targs - t(whs6) %*% x
d6

step6 <- stepfn(whs6, d6, beta6, x)
step6

# iter7 ----
beta7 <- beta6 + step6
delta7 <- get_delta(wh6, beta7, x) 
whs7 <- get_weights(beta7, delta7, x)
whs7

(ws7=colSums(whs7))
(wh7=rowSums(whs7))

d7 <- targs - t(whs7) %*% x
d7

step7 <- stepfn(whs7, d7, beta7, x)
step7

# iter8 ----
beta8 <- beta7 + step7
delta8 <- get_delta(wh7, beta8, x) 
whs8 <- get_weights(beta8, delta8, x)
whs8

(ws8=colSums(whs8))
(wh8=rowSums(whs8))

d8 <- targs - t(whs8) %*% x
d8

step8 <- stepfn(whs8, d8, beta8, x)
step8


# iter9 ----
beta9 <- beta8 + step8
delta9 <- get_delta(wh8, beta9, x) 
whs9 <- get_weights(beta9, delta9, x)
whs9

(ws9=colSums(whs9))
(wh9=rowSums(whs9))

d9 <- targs - t(whs9) %*% x
d9

step9 <- stepfn(whs9, d9, beta9, x)
step9


# end iters
d1; d2; d3; d4; d5; d6; d7; d8; d9
sum(d1^2); sum(d2^2); sum(d3^2); sum(d4^2); sum(d5^2); sum(d6^2); sum(d7^2); sum(d8^2); sum(d9^2)
step1; step2; step3; step4; step5; step6; step7; step8; step9
whs; whs1; whs2; whs3; whs4; whs5; whs6; whs7; whs8; whs9
targs; t(whs9) %*% x



# loop ----
stepfn <- stepa
stepfn <- stepb # better
stepfn <- stepc
stepfn <- stepd

#.. setup ----
# xscale <- 10 and step_scale <- 10 works well with this
h <- 4 # number households
k <- 2 # number characteristics
s <- 3 # number states

x <- matrix(c(5, 6,
              7, 4,
              9, 2,
              3, 1), nrow=h, byrow=TRUE)

whs <- matrix(c(4, 3, 4,
                7, 12, 1,
                20, 6, 4,
                6, 4, 5), nrow=4, byrow=TRUE)

# (xxp <- x %*% t(x))
# (xpx <- t(x) %*% x)
(wh=rowSums(whs))
(ws=colSums(whs))

(targs <- t(whs) %*% x) # s x k

#.. end setup ----
betavec <- as.vector(ebeta)
xmat <- x
targets <- targs
f <- function(betavec, wh, xmat, targets){
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  as.vector(d)
}


fvec <- function(betavec, wh, xmat, targets){
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  sse <- sum(d^2)
  sse
}
g <- grad(f, x=betavec, wh=wh, xmat=xmat, targets=targs)
h <- hessian(f, x=betavec, wh=wh, xmat=xmat, targets=targs)

ebeta <- matrix(0, nrow=s, ncol=k) # tpc uses 0 as beta starting point
edelta <- get_delta(wh, ebeta, x) # tpc uses initial delta based on initial beta 
f(as.vector(ebeta), wh, xmat=x, targs)
gval <- grad(f, x=as.vector(ebeta), wh=wh, xmat=xmat, targets=targs)
xxp <- x %*% t(x)
xpx <- t(x) %*% x

fpod <- matrix(gval, nrow=nrow(targets), byrow=FALSE)

1 / fpod %*% t(d)

# i <- 1

for(i in 1:20){
  print(i)
  ewhs <- get_weights(ebeta, edelta, x)
  ews <- colSums(ewhs)
  ewh <- rowSums(ewhs)
  
  etargs <- t(ewhs) %*% x
  d <- targs - etargs
  d
  sum(d^2)
  # exit if good
  
  gval <- grad(f, x=as.vector(ebeta), wh=wh, xmat=xmat, targets=targs)
  fpod <- matrix(gval, nrow=nrow(targets), byrow=FALSE)
  step <- sweep(1 / fpod, 1, t(d), "*")
  
  # step <- stepfn(ewhs, d, ebeta, x)
  step
  ebeta <- ebeta + step
  edelta <- get_delta(ewh, ebeta, x)
}
d

whs
ewhs %>% round(1)
ws; ews %>% round(1)
wh; ewh



# nlm ----
f <- function(x) sum((x-1:length(x))^2)
nlm(f, c(10,10))
nlm(f, c(10,10), print.level = 2)
utils::str(nlm(f, c(5), hessian = TRUE))

f <- function(x, a) sum((x-a)^2)
nlm(f, c(10,10), a = c(3,5))
f <- function(x, a)
{
  res <- sum((x-a)^2)
  attr(res, "gradient") <- 2*(x-a)
  res
}
nlm(f, c(10,10), a = c(3,5))



fvec <- function(betavec, wh, xmat, targets){
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  sse <- sum(d^2)
  sse
}

beta0 <- matrix(0, nrow=s, ncol=k) # tpc uses 0 as beta starting point
delta0 <- get_delta(ebeta0, x) # tpc uses initial delta based on initial beta 

betavec <- as.vector(beta0)
fvec(betavec, wh, xmat=x, targets=targs)

res <- nlm(fvec, betavec, wh=wh, xmat=x, targets=targs)
res

bres <- matrix(res$estimate, nrow=nrow(targs), byrow=FALSE)
dres <- get_delta(wh, bres, x)
wres <- get_weights(bres, dres, x)
wres
colSums(wres)
rowSums(wres)
tres <- t(wres) %*% x
d <- targs - tres
d
sum(d^2)

