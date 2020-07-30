

library(nlsr)


pran <- make_problem(h=6, k=3, s=4)
names(pran)

nlfb(start, resfn, jacfn=NULL, trace=FALSE, lower=-Inf, upper=Inf,
     maskidx=NULL, weights=NULL, data=NULL, control, ...)

pran <- make_problem(h=100, k=5, s=10)
ans <- nlfb(start=rep(0, pran$k * pran$s), trace=TRUE, sse_fn, wh=pran$wh, xmat=pran$xmat, targets=pran$targets)
str(ans)

out <- nls.lm(par=rep(0, pran$k * pran$s), fn = diff_vec, 
               control=nls.lm.control(maxiter=50, nprint=1),
               wh=pran$wh, xmat=pran$x, targets=pran$targets)
str(out)
cbind(out$par, ans$coefficients) %>% round(3)

sse_fn(out$par, pran$wh, pran$xmat, pran$targets)
sse_fn(ans$coefficients, pran$wh, pran$xmat, pran$targets)


shobbs.res <- function(x){ # scaled Hobbs weeds problem -- residual
  # This variant uses looping
  if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
  y <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
         38.558, 50.156, 62.948, 75.995, 91.972)
  tt <- 1:12
  res <- 100.0*x[1]/(1+x[2]*10.*exp(-0.1*x[3]*tt)) - y
}

shobbs.jac <- function(x) { # scaled Hobbs weeds problem -- Jacobian
  jj <- matrix(0.0, 12, 3)
  tt <- 1:12
  yy <- exp(-0.1*x[3]*tt)
  zz <- 100.0/(1+10.*x[2]*yy)
  jj[tt,1] <- zz
  jj[tt,2] <- -0.1*x[1]*zz*zz*yy
  jj[tt,3] <- 0.01*x[1]*zz*zz*yy*x[2]*tt
  return(jj)
}

cat("try nlfb\n")
st <- c(b1=1, b2=1, b3=1)
low <- -Inf
up <- Inf


ans1 <- nlfb(st, shobbs.res, shobbs.jac, trace=TRUE)
ans1

cat("No jacobian function -- use internal approximation\n")
ans1n <- nlfb(st, shobbs.res, trace=TRUE, control=list(watch=TRUE)) # NO jacfn
ans1n

str(ans1n)
a <- shobbs.res(ans1n$coefficients)
str(a)

nlfb(start, resfn, jacfn=NULL, trace=FALSE, lower=-Inf, upper=Inf,
     maskidx=NULL, weights=NULL, data=NULL, control, ...)




library(BB)
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
                        seed=1234))
rosbkext <- function(x){
  # Extended Rosenbrock function
  n <- length(x)
  j <- 2 * (1:(n/2))
  jm1 <- j - 1
  sum(100 * (x[j] - x[jm1]^2)^2 + (1 - x[jm1])^2)
}
p0 <- rnorm(50)
spg(par=p0, fn=rosbkext)
BBoptim(par=p0, fn=rosbkext)


# examine BB ----
pran <- make_problem(h=6, k=3, s=4)
pran <- make_problem(h=100, k=10, s=15)
pran <- make_problem(h=1000, k=20, s=25)
pran <- make_problem(h=5000, k=25, s=35)
names(pran)

t1 <- proc.time()
opt.bb <- BBoptim(par=rep(0, pran$k * pran$s), fn=sse_fn, 
                  method=3,
                  wh=pran$wh, xmat=pran$x, targets=pran$targets)
t2 <- proc.time()
t2 - t1

gr <- function(x, wh, xmat, targets){
  numDeriv::grad(sse_fn, x, wh=wh, xmat=xmat, targets=targets)
}
gr(rep(0, pran$k * pran$s), wh=pran$wh, xmat=pran$x, targets=pran$targets)

t1 <- proc.time()
opt.bb <- BBoptim(par=rep(0, pran$k * pran$s), fn=sse_fn, # gr=gr,
                  method=3, control=list(M=10),
                  wh=pran$wh, xmat=pran$x, targets=pran$targets)
t2 <- proc.time()
t2 - t1

pran <- p
t1 <- proc.time()
opt.bb <- BBoptim(par=rep(0, pran$k * pran$s), fn=sse_fn, gr=grad_sse,
                  method=3, control=list(M=10),
                  wh=pran$wh, xmat=pran$x, targets=pran$targets)
t2 <- proc.time()
t2 - t1

t3 <- proc.time()
opt.lm <- nls.lm(par=rep(0, pran$k * pran$s), fn = diff_vec, 
              control=nls.lm.control(factor=20, maxiter=50, nprint=1),
              wh=pran$wh, xmat=pran$x, targets=pran$targets)
t4 <- proc.time()
t4 - t3

cbind(opt.bb$par, opt.lm$par)


sse_fn(opt.bb$par, pran$wh, pran$xmat, pran$targets)
sse_fn(opt.lm$par, pran$wh, pran$xmat, pran$targets)
