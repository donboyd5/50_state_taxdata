
# libraries ----
source(here::here("include", "libraries.r"))
library(numDeriv) # grad, jacobian, and hessian
library(ipoptr)
library(nloptr)
library(minpack.lm) # nls.lm 


# functions ----
source(here::here("include", "functions_prep_dev.r")) # soon we will replace functions_prep.r with the dev version
source(here::here("include", "functions_prep_data.r"))
source(here::here("include", "functions_poisson_model.r"))
source(here::here("include", "functions_utilities.r"))


# get problem ----
pran <- make_problem(h=6, k=3, s=4)
pran <- make_problem(h=100, k=5, s=10) # ss 100
pran <- make_problem(h=1000, k=10, s=20)
pran <- make_problem(h=5e3, k=20, s=30)
pran <- make_problem(h=10e3, k=20, s=50)


# problem notes ----
# 6, 3, 4 sse 122.7315 0.1651873 2.941963e-07 8.263819e-18 1.123338e-27 INCREASE 3.382636e-27 4.897248e-27


# loop ----
prob <- pran

prob <- pacs
prob <- scale_problem(pacs, scale_goal = 1000)
prob <- scale_problem_mdn(pacs, scale_goal = 1)

(h <- prob$h)
(k <- prob$k)
(s <- prob$s)

wh <- prob$wh
whs <- prob$whs # what we are seeking
xmat <- prob$xmat

xpx <- t(xmat) %*% xmat # can be done before loop
ixpx <- solve(xpx) # also can be done before loop

targets <- prob$targets

ebeta0 <- vtom(rep(0, k * s), nrow=s)
ebeta0

step_scale <- 1
step_scale <- 30e3 # 0.5
ebeta <- ebeta0
#  12k 2.06703e+04"
# 13k 1.54422e+04

# ebeta <- fbeta
step_scale <- .00001

for(i in 1:1000){
  edelta <- get_delta(wh, ebeta, xmat) # wh and xmat are constant
  edelta
  
  ewhs <- get_weights(ebeta, edelta, xmat)
  ewhs
  ews <- colSums(ewhs) # needed for analytic jacobian calculation but not for finite differences
  
  etargets <- t(ewhs) %*% xmat
  etargets
  
  d <- targets - etargets
  sse <- sum(d^2)
  d
  d / targets * 100
  sse
  print(sprintf("i: %3i  sse: %.5e", i, sse))
  if(sse < 1e-6) break

  # compute the step
  
  # TPC jacobian approach, but with numeric jacobian
  # diff_vec takes arguments betavec, wh, xmat, targets
  
  # betavec <- as.vector(ebeta)
  # jac <- jacobian(diff_vec, x=betavec, wh=wh, xmat=xmat, targets=targets)
  # jac
  # ijac <- solve(jac)
  # ijac
  # step.j <- ijac %*% as.vector(d)
  # step.j <- vtom(step.j, nrows=s)

  # ad hoc step
  # xpx <- t(xmat) %*% xmat # can be done before loop
  # ixpx <- solve(xpx) # also can be done before loop
  step.gr1 <- numDeriv::grad(sse_fn, as.vector(ebeta), wh=wh, xmat=xmat, targets=targets) # gr(as.vector(ebeta), wh, xmat, targets)
  step.gr2 <- vtom(step.gr1, s)
  step.gr <- step.gr2 # * d / 100000
  
  step.ah <- -(1 / ews) * d %*% ixpx
  
  # step <- step.ah
  step <- step.gr
  step
  
  ebeta <- ebeta - step * step_scale
  ebeta
}

temp <- ebeta

ebeta <- temp

i
sse
ebeta
targets %>% round(2)
etargets %>% round(2)

(targets - etargets) %>% round(2)

targets.u <- sweep(prob$targets, 2, prob$scale_factor, "/")
etargets.u <- sweep(etargets, 2, prob$scale_factor, "/")
targets.u %>% round(2)
etargets.u %>% round(0)

(rowSums(ewhs) - wh) %>% round(1)


step.j; step.ah
step.j / step.ah


ebeta
edelta

gr <- function(par, wh, xmat, targets){
  
}

prob <- pran
res <- optim(rep(0, prob$s * prob$k), sse_fn, gr = gr, method = "BFGS", 
             control=list(trace=2, fnscale=10), wh=prob$wh, xmat=prob$xmat, targets=prob$targets)

pran <- make_problem(h=5e3, k=20, s=30)
pran <- make_problem(h=10e3, k=20, s=50)

prob <- scale_problem(pran, scale_goal = 1000)

prob <- pran

a <- proc.time()
res2 <- optim(rep(0, prob$s * prob$k), sse_fn, method = "BFGS", 
             control=list(trace=2, fnscale=10), wh=prob$wh, xmat=prob$xmat, targets=prob$targets)
b <- proc.time()
b - a

# scaling notes
# https://www.nag.com/numeric/fn/manual/pdf/c09/c09int_fn04.pdf
# Normally, users should restrict themselves to linear transformations of variables, although occasionally
# nonlinear transformations are possible. The most common such transformation (and often the most
# appropriate) is of the form
# 
#    xnew = Dxold,
# 
# where D is a diagonal matrix with constant coefficients. Our experience suggests that more use should
# be made of the transformation 
# 
#   xnew = Dxold + v,
# 
# where v is a constant vector.

# nls.lm ----
# get a starting point
sp <- nls.lm(par=rep(0, prob$k * prob$s), fn = diff_vec, 
             control=nls.lm.control(maxiter=100, nprint=2),
             wh=prob$wh, xmat=prob$xmat, targets=prob$targets)

x0 <- rep(0, prob$k * prob$s)
x0 <- as.vector(ebeta)

a <- proc.time()
sp <- nls.lm(par=x0, fn = diff_vec, 
             control=nls.lm.control(factor=20, maxiter=100, nprint=2),
             wh=prob$wh, xmat=prob$xmat, targets=prob$targets)
b <- proc.time()
b - a

sp <- nls.lm(par=rep(0, prob$k * prob$s), fn = diff_vec, 
             control=nls.lm.control(factor=18, maxiter=100, nprint=2),
             wh=prob$wh, xmat=prob$xmat, targets=prob$targets)

sp2 <- nls.lm(par=sp$par, fn = diff_vec, 
             control=nls.lm.control(maxiter=100, nprint=2),
             wh=prob$wh, xmat=prob$xmat, targets=prob$targets)

diff_vec(as.vector(ebeta), prob$wh, prob$xmat, prob$targets)

str(sp)
sp$rsstrace
sp$par
sse_fn(sp$par, prob$wh, prob$xmat, prob$targets)

prob$s; prob$k; prob$h

fbeta <- vtom(sp$par, prob$s)
fbeta

fdelta <- get_delta(prob$wh, fbeta, prob$xmat) # wh and xmat are constant
fdelta

ewhs <- get_weights(fbeta, fdelta, prob$xmat)
ewhs

targets.u <- sweep(prob$targets, 2, prob$scale_factor, "/")
etargets <- t(ewhs) %*% prob$xmat
etargets.u <- sweep(etargets, 2, prob$scale_factor, "/")
targets.u %>% round(2)
etargets.u %>% round(0)
d <- targets.u - etargets.u
sse <- sum(d^2)
d
(d / targets.u * 100) %>% round(2)
sse


ssed <- deriv(~ sum(diff_vec(betavec, wh, xmat, targets)^2), 
              c("betavec", "xmat", "wh", "xmat", "targets"), 
              function.arg = TRUE)


sse2 <- deriv(~ -log(lambda * dnorm((x-mu1)/s1)/s1 + (1-lambda)*dnorm((x-mu2)/s2)/s2), 
                 c("mu1", "mu2", "s1", "s2", "lambda"), 
                 function.arg = TRUE)

prob <- pran
system.time(res <- optim(rep(0, prob$s * prob$k), sse_fn, gr = NULL, method = "BFGS", 
             control=list(trace=2, fnscale=1000), wh=prob$wh, xmat=prob$xmat, targets=prob$targets))
str(res)
