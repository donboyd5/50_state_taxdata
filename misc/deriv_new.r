

source(here::here("include", "libraries.r"))
library(numDeriv) # grad, jacobian, and hessian
library(ipoptr)
library(nloptr)
library(minpack.lm) # nls.lm 



f <- function(betavec, delta, wh, xmat, targets){
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  as.vector(d)
}


fd <- function(betavec, delta, wh, xmat, targets){
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  as.vector(d)
}


# start #1 ----
k <- 2
s <- 3
h <- 4
ks <- k * s

# h x k
xhk <- matrix(c(1, 2, 
               3, 4,
               5, 6,
               7, 8), ncol=k, byrow=TRUE)
xhk

wh <- c(6, 15, 21, 30)

beta <- matrix(c(.1, .2,
                 .3, .4,
                 .5, .6), ncol=k, byrow=TRUE)

delta <- get_delta(wh, beta, xhk)
delta

whs <- get_weights(beta, delta, xhk)
whs # h x s

# 1 x h
rowSums(whs)
# 1 x s
(ws <- colSums(whs))

etargs <- t(whs) %*% xhk
etargs

targs <- matrix(c(4, 8,
                  25, 40,
                  320, 400), ncol=k, byrow=TRUE)
targs


dm <- targs - etargs
dm

dv <- as.vector(dm)
dv

bv <- as.vector(beta)

jval <- jacobian(f, x=as.vector(beta), delta=delta, wh=wh, xmat=xhk, targets=targs)
jval
jvald <- jacobian(fd, x=as.vector(beta), delta=delta, wh=wh, xmat=xhk, targets=targs)
jvald
solve(jval)
solve(jvald)

# xpx <- t(xhs) %*% xhs
# xxp <- xhs %*% t(xhs)

#.. this is it #1 ----
(x2 <- xhk * xhk)
(x12 <- xhk[, 1] * xhk[, 2])
x2; x12

(p1 <- (t(whs) %*% x2) %>% as.vector)
(p2 <- t(whs) %*% x12)

jval

# do it efficiently
(p1 <- (t(whs) %*% x2))
(p2 <- t(whs) %*% x12)
t(whs) %>% cbind(p1, p2)

t(whs[, 1]) %*% cbind(x2, x12)
t(whs[, 2]) %*% cbind(x2, x12, x2)

diag(as.vector(t(whs) %*% x2))
p1 <- t(whs) %*% x2
p2 <- t(whs) %*% x12
diag(p1[, 1])
diag(p1[, 2])
diag(p2[, 1])



# start #2 add 1 to each k s h ----
k <- 3
s <- 4
h <- 5
ks <- k * s

# h x k
xhk <- matrix(1:(h*k), ncol=k, byrow=TRUE)
xhk

wh <- seq(10, 20, length.out=h)

beta <- matrix(seq(.1, ks*.1, length.out=ks), ncol=k, byrow=TRUE)
beta

delta <- get_delta(wh, beta, xhk)
delta

whs <- get_weights(beta, delta, xhk)
whs # h x s

# 1 x h
rowSums(whs)
# 1 x s
(ws <- colSums(whs))

etargs <- t(whs) %*% xhk
etargs

set.seed(1234)
targs <- etargs + runif(ks, -2, 2)
targs


dm <- targs - etargs
dm

dv <- as.vector(dm)
dv

bv <- as.vector(beta)

jval <- jacobian(f, x=as.vector(beta), delta=delta, wh=wh, xmatrix=xhk, targets=targs)
jval

# xpx <- t(xhs) %*% xhs
# xxp <- xhs %*% t(xhs)

# this is it #2 ----
s; k; h
(x2 <- xhk * xhk)
(x12 <- xhk[, 1] * xhk[, 2])
x2; x12

(p1 <- (t(whs) %*% x2) %>% as.vector)
(p2 <- t(whs) %*% x12)

jval

# do it efficiently
(p1 <- (t(whs) %*% x2))
(p2 <- t(whs) %*% x12)
t(whs) %>% cbind(p1, p2)

t(whs[, 1]) %*% cbind(x2, x12)
t(whs[, 2]) %*% cbind(x2, x12, x2)

diag(as.vector(t(whs) %*% x2))
p1 <- t(whs) %*% x2
p2 <- t(whs) %*% x12
diag(p1[, 1])
diag(p1[, 2])
diag(p2[, 1])

jac <- bind_rows()

# fns ----

f2 <- function(betavec, wh, xmat, targets, delta=NULL){
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  if(is.null(delta)) delta <- get_delta(wh, beta, xmat) 
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  sum(d^2)
}

gf2 <- function(betavec, wh, xmat, targets, delta=NULL){
  numDeriv::grad(func=f2, x=betavec, wh=wh, xmat=xmat, targets=targets, delta=delta)
}

hf2 <- function(betavec, wh, xmat, targets, delta=NULL){
  numDeriv::hessian(func=f2, x=betavec, wh=wh, xmat=xmat, targets=targets, delta=delta)
}


# loop ----
xhk.s <- sweep(xhk, 2, colSums(xhk), "/") * 10
xhk.s; colSums(xhk.s)
xhk.s / xhk
targs.s <- sweep(targs, 2, colSums(xhk), "/") * 10
targs.s / targs

xhk.s <- xhk; targs.s <- targs

xpx <- t(xhk.s) %*% xhk.s
invxpx <- solve(xpx)

n.beta <- vtom(rep(0, ks), nrow=s)
n.beta
wh
round(whs, 2)
for(i in 1:100){
  n.delta <- get_delta(wh, n.beta, xhk.s) # wh and xhk are const
  n.delta
  # delta_prior <- delta
  n.whs <- get_weights(n.beta, n.delta, xhk.s)
  n.whs
  n.ws <- colSums(n.whs)
  n.wh <- rowSums(n.whs) # not needed I think
  n.wh
  n.etargs <- t(n.whs) %*% xhk.s
  d <- targs.s - n.etargs
  sse <- sum(d^2)
  sse
  d
  #if(i <=10 | i %% 10==0) {print(i); print(sse)}
  #if(sse < 2) break
  # jac <- jacobian(f, x=as.vector(n.beta), delta=n.delta, wh=wh, xmat=xhk, targets=targs)
  # jac
  # step <- solve(jac) %*% as.vector(d)
  # # step <- as.vector(d) %*% solve(jac)
  # step <- vtom(step, nrows = s)
  # step
  betavec <- as.vector(n.beta)
  # g2 <- gf2(betavec, wh, xhk.s, targs.s, delta=delta_prior)
  # h2 <- hf2(betavec, wh, xhk.s, targs.s, delta=delta_prior)
  # ih2 <- solve(h2)
  # step <- ih2 %*% g2
  # step <- vtom(step, nrows = s)
  # step
  step <- -(1 / n.ws) * d %*% invxpx * 1.5
  # -(1 / step_inputs$ews) * step_inputs$d %*% step_inputs$invxpx * step_inputs$step_scale
  # step
  n.beta <- n.beta - step
  n.beta
  # delta_prior <- delta
}
i
sse
targs.s
n.etargs

out <- nls.lm(par=rep(0, ks), fn = fz, wh=wh, xmat=xhk.s, targets=targs.s)
out
out$rsstrace

n.beta
f.beta <- vtom(out$par, s)
f.delta <- get_delta(wh, f.beta, xhk.s)
f.whs <- get_weights(f.beta, f.delta, xhk.s)
f.etargs <- t(f.whs) %*% xhk.s
targs.s
f.etargs %>% round(1)


f2(betavec, delta, wh, xmat, targets)
# n.beta <- vtom(rep(0, ks), nrow=s)
betavec <- rep(0, ks)
n.beta <- vtom(betavec, nrow=s)
wh
n.delta <- get_delta(wh, n.beta, xhk.s) # wh and xhk are const
f2(betavec, n.delta, wh, xhk.s, targs.s)
  
out <- nls.lm(par=betavec, fn = fz, wh=wh, xmat=xhk.s, targets=targs.s)
str(out)
out$par
out$rsstrace

fbeta <- vtom(out$par, nrows = s)
fdelta <- get_delta(wh, fbeta, xhk.s)
fwhs <- get_weights(fbeta, fdelta, xhk.s)
fetargs <- t(fwhs) %*% xhk.s
fetargs %>% round(1)
targs.s
d <- targs.s - fetargs
d


prob <- pran

prob <- pacs
prob <- scale_problem(pacs, scale_goal=100)

betavec <- rep(0, length(prob$targets))
d0 <- fz(betavec, wh=prob$wh, xmat=prob$x, targets=prob$targets)
sum(d0^2)

betavec <- as.vector(result$best_ebeta)

a <- proc.time()
# out <- nls.lm(par=betavec, fn = fz, 
#               control=nls.lm.control(maxiter=50, nprint=1),
#               wh=prob$wh, xmat=prob$x, targets=prob$targets)
out3 <- nls.lm(par=out2$par, fn = fz, 
              control=nls.lm.control(maxiter=50, nprint=1),
              wh=prob$wh, xmat=prob$x, targets=prob$targets)
b <- proc.time()
b - a
# 5.1 hours 40k records 50 states, 13 targets per state

out2$rsstrace

str(out)
out$niter
out$deviance
out$rsstrace

d0 <- fz(betavec, wh=prob$wh, xmat=prob$x, targets=prob$targets)
sum(d0^2)
df <- fz(out$par, wh=prob$wh, xmat=prob$x, targets=prob$targets)
# df; out$fvec
sum(df^2)

fbeta <- vtom(out$par, nrows = nrow(prob$targets))
fdelta <- get_delta(prob$wh, fbeta, prob$x)
fwhs <- get_weights(fbeta, fdelta, prob$x)
fetargs <- t(fwhs) %*% prob$x
fetargs %>% round(1)
prob$targets %>% round(1)
d <- prob$targets - fetargs
d %>% round(2)
(d / prob$targets * 100) %>% round(2) 


fz <- function(betavec, wh, xmat, targets){
  # return the distances from targets as a vector using the OLD delta
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  as.vector(d)
}

