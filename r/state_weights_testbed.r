
# code folding ----
# alt-o, shift-alt-o
# alt-l, shift-alt-l
# alt-r

# notes ----


# libraries ----
source(here::here("include", "libraries.r"))
library(numDeriv) # grad, jacobian, and hessian
library(ipoptr)
library(nloptr)
library(minpack.lm) # nls.lm 
library(nleqslv) # nleqslv

devtools::session_info()
(.packages()) %>% sort


# globals ----
dbox <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/"
(fns <- paste0(c("acs_10krecs_5states", "acs_100krecs_20states", "acs_200krecs_50states", "acs_400krecs_50states"), ".rds"))


# functions ----
source(here::here("include", "functions_prep_dev.r")) # soon we will replace functions_prep.r with the dev version
source(here::here("include", "functions_prep_data.r"))
source(here::here("include", "functions_poisson_model.r"))
source(here::here("include", "functions_utilities.r"))


# prepare an acs problem ---
#.. get and structure the desired data ----
fns
acslist <- prep_ACS_sample(fns[1])
possible_target_vars <- acslist$possible_target_vars

#.. define target vars as a subset of possible vars -- for example, one of the following ----
possible_target_vars
target_vars <- possible_target_vars[c(1, 3, 6, 7)]
target_vars <- possible_target_vars[1:10] # best 6; can't do 1:7
ivars <- c(1:3, 5:8)
target_vars <- possible_target_vars[ivars]

target_vars <- setdiff(possible_target_vars, c("pap_sum", "ssip_sum", "intp_sum", "otherincp_sumneg"))
target_vars <- setdiff(possible_target_vars, c("pap_sum", "ssip_sum", "intp_sum"))
target_vars <- setdiff(possible_target_vars, c("pap_sum", "ssip_sum"))
target_vars <- setdiff(possible_target_vars, c("pap_sum"))
target_vars <- possible_target_vars

#.. prepare the problem ----
pacs <- make_acs_problem(acslist, target_incgroup=2, target_vars)


# make a desired random problem (not ACS) ----
pran <- make_problem(h=8, k=2, s=3)
pran <- make_problem(h=100, k=5, s=10)


# SOLVE the chosen problem ----
res <- solve_poisson(pran)
res <- solve_poisson(pran, scale=TRUE, scale_goal = 100)
res <- solve_poisson(pran, step_method="finite_diff")
res <- solve_poisson(pran, step_method="finite_diff", start=res$best_ebeta)
res <- solve_poisson(pran, scale=TRUE, scale_goal = 100, step_method="finite_diff")
res <- solve_poisson(pran, step_scale = 100)
res <- solve_poisson(pran, scale=TRUE, scale_goal = 100, step_scale=7)

res <- solve_poisson(pran, maxiter=20)
res <- solve_poisson(pran, maxiter=20, step_method="finite_diff")



#.. acs 1 ----
resah <- solve_poisson(pacs, scale=TRUE, scale_goal = 100, step_scale=2e3) # acs 1, 4 targets
resah <- solve_poisson(pacs, scale=TRUE, scale_goal = 100, step_scale=2e3, maxiter=30) # acs 1, 4 targets
resfd <- solve_poisson(pacs, scale=TRUE, scale_goal = 100, step_method="finite_diff")
resfd <- solve_poisson(pacs, scale=TRUE, scale_goal = 100, step_method="finite_diff", start=resah$best_ebeta)
resfd <- solve_poisson(pacs, step_method="finite_diff")

resah <- solve_poisson(pacs, step_scale=.75e3, scale=TRUE, scale_goal=100) # acs 1, 13 targets
resah <- solve_poisson(pacs, step_scale=1000, scale=TRUE, scale_goal=100, maxiter = 200) # acs 1, 13 targets

resah <- solve_poisson(pacs, step_scale=.75e3, scale=TRUE, scale_goal=100, maxiter=20)
resfd <- solve_poisson(pacs, step_scale=.75e3, scale=TRUE, scale_goal=100, step_method = "finite_diff")
resfd <- solve_poisson(pacs, step_scale=.75e3, scale=TRUE, scale_goal=100, step_method = "finite_diff", start=resah$best_ebeta)

#.. acs 2 ----
res <- solve_poisson(pacs, step_scale=6e3) # acs 2, 4 targets
system.time(res <- solve_poisson(pacs, step_scale=8000, maxiter=20000)) # acs 2, 13 targets
system.time(res <- solve_poisson(pacs, step_scale=8000, maxiter=20000, scale=TRUE, scale_goal=100)) # acs 2, 13 targets
system.time(res <- solve_poisson(pacs, step_scale=8e3, scale=TRUE, scale_goal=100, maxiter=20000, tol=1e-3)) # acs 3, 13 targets

#.. acs 3 ----
res <- solve_poisson(pacs, step_scale=11e3, tol=1e-3) # acs 3, 4 targets
res <- solve_poisson(pacs, step_scale=11e3, scale=TRUE, scale_goal=100, maxiter=200, tol=1e-3) # acs 3, 13 targets

#.. acs 4 ----
# 44 iterations 25 seconds
system.time(res <- solve_poisson(pacs, step_scale=15e3, scale=TRUE, scale_goal=100, maxiter=500, tol=1e-4)) # acs 3, 13 targets


# examine results ----
result <- res
names(result)
# str(result)
result$total_seconds
result$problem_unscaled$h; result$problem_unscaled$k; result$problem_unscaled$s
result$iter
result$max_rel_err
result$sse
result$sse_vec
result$d %>% round(4)
max(abs(result$d))

# check weights
(result$ewh - result$problem_unscaled$wh) %>% round(2) # total household weights

# check targets
result$problem_unscaled$targets %>% round(0)
result$etargets %>% round(0)
pdiff <- (result$etargets / result$problem_unscaled$targets * 100 - 100) %>% round(2)
pdiff
max(ifelse(is.infinite(abs(pdiff)), 0, abs(pdiff)))
sort(-abs(pdiff))

res$etargets

res$best_ebeta
res$ewhs

# solve a problem using nls.lm ----
prob <- pran
prob <- pacs; betavec <- rep(0, length(prob$targets))

prob <- scale_problem(pacs, scale_goal=1e5)
prob <- scale_problem_mdn(pacs, scale_goal=1000)
prob <- p

betavec <- as.vector(ebeta)
betavec <- rep(0, length(prob$targets))
d0 <- diff_vec(betavec, wh=prob$wh, xmat=prob$x, targets=prob$targets)
sum(d0^2)

betavec <- as.vector(result$best_ebeta)

reldiff <- function(par, prob){
  ((diff_vec(par, prob$wh, prob$xmat, prob$targets) / prob$targets) * 100) %>% round(3)
}

prob$targets
dw <- prob$targets / 100
dw <- ifelse(prob$targets!=0, 100 / prob$targets, 1)
dw
prob$targets * dw

a <- proc.time()
out3 <- nls.lm(par=tmp, fn=diff_vec, 
               control=nls.lm.control(maxiter=50, nprint=1, factor=100),
               wh=prob$wh, xmat=prob$x, targets=prob$targets, dweights=as.vector(dw))
b <- proc.time()
b - a # 5398.83  secs 1.5 hours
out3$rsstrace

a2 <- proc.time()
out32 <- nls.lm(par=out3en$x, fn=diff_vec, 
               control=nls.lm.control(maxiter=50, nprint=1, factor=100),
               wh=prob$wh, xmat=prob$x, targets=prob$targets, dweights=as.vector(dw))
b2 <- proc.time()
b2 - a2 # 5398.83  secs 1.5 hours

reldiff(out3$par, pacs) %>% round(2)
sse_fn(out3$par, wh=prob$wh, xmat=prob$x, targets=prob$targets) 


# djb this may be best when Jacobian is too expensive to calculate on every iteration ----
a5 <- proc.time()
out3e <- nleqslv(x=out3e$x, fn=diff_vec, jac=NULL,
                 wh=prob$wh, xmat=prob$x, targets=prob$targets, dweights=as.vector(dw),
                 method = c("Broyden"),
                 global = c("dbldog"),
                 xscalm = c("auto"),
                 jacobian=FALSE,
                 control=list(maxit=2000, trace=1, allowSingular=TRUE))
b5 <- proc.time()
b5 - a5 # 3 mins 3 iterations
names(out3e)
out3e$fvec %>% round(2)
prob$targets %>% round(1)
out3e$termcd
out3e$message

tmp <- out3e$x

sse_fn(out3e$x, wh=prob$wh, xmat=prob$x, targets=prob$targets) 
reldiff(out3e$x, pacs) %>% round(2)


# djb this appears to be the best approach so far for problems where it is reasonable to use Newton jacobian ----
a5n <- proc.time()
out3en <- nleqslv(x=out3en$x, fn=diff_vec, jac=NULL,
                 wh=prob$wh, xmat=prob$x, targets=prob$targets, dweights=as.vector(dw),
                 method = c("Newton"),
                 global = c("dbldog"),
                 xscalm = c("auto"),
                 jacobian=FALSE,
                 control=list(maxit=50, trace=1, allowSingular=TRUE))
b5n <- proc.time()
b5n - a5n 

sse_fn(out3en$x, wh=prob$wh, xmat=prob$x, targets=prob$targets) 
reldiff(out3en$x, pacs) %>% round(2)
# 2.970924e+04


a5nb <- proc.time()
out3enb2 <- nleqslv(x=out3e$x, fn=diff_vec, jac=NULL,
                  wh=prob$wh, xmat=prob$x, targets=prob$targets,
                  method = c("Newton"),
                  global = c("dbldog"),
                  xscalm = c("auto"),
                  jacobian=FALSE,
                  control=list(maxit=50, trace=1, allowSingular=TRUE))
b5nb <- proc.time()
b5nb - a5nb
names(out3enb)
out3enb2$fvec
# crossprod(out3enb2$fvec)/2
# Fnorm is sum(out3enb2$fvec^2)/2
reldiff(out3enb2$x, pacs)

anl1 <- proc.time()
nl1 <- nleqslv(x=out3e$x, fn=diff_vec, jac=NULL,
                  wh=prob$wh, xmat=prob$x, targets=prob$targets, dweights=as.vector(dw),
                  method = c("Newton"),
                  global = c("none"),
                  xscalm = c("auto"),
                  jacobian=FALSE,
                  control=list(maxit=20, trace=1, allowSingular=TRUE))
bnl1 <- proc.time()
bnl1 - anl1


az <- proc.time()
out3z <- nls.lm(par=rep(1, length(prob$targets)), fn = diff_vec, 
               control=nls.lm.control(maxiter=30, nprint=1, factor=100),
               wh=prob$wh, xmat=prob$x, targets=prob$targets)
bz <- proc.time()
bz - az # 5398.83  secs 1.5 hours
out3z$rsstrace
reldiff(out3z$par, pacs)


a1 <- proc.time()
out3a <- nls.lm(par=betavec, fn = diff_vec, 
               control=nls.lm.control(maxiter=50, nprint=1, factor=20),
               wh=prob$wh, xmat=prob$x, targets=prob$targets)
b1 <- proc.time()
b1 - a1
out3a$rsstrace


library(pracma)
a2 <- proc.time()
out3b <- lsqnonlin(fun=diff_vec, x0=betavec, options = list(tau=1e-4, maxeval=100), wh=prob$wh, xmat=prob$x, targets=prob$targets)
b2 <- proc.time()
b2 - a2 # 
# names(out3b)
out3b$errmess
reldiff(out3b$x, pacs)


a3 <- proc.time()
out3c <- newtonsys(Ffun= diff_vec, x0=betavec, maxiter=3, wh=prob$wh, xmat=prob$x, targets=prob$targets)
b3 <- proc.time()
b3 - a3 # 18 mins 3 iterations
names(out3c)
out3c$zero
out3c$fnorm
out3c$niter
reldiff(out3c$zero, pacs)


a4 <- proc.time()
out3d <- broyden(Ffun=diff_vec, x0=betavec, maxiter=20, wh=prob$wh, xmat=prob$x, targets=prob$targets)
b4 <- proc.time()
b4 - a4 # 6 mins 2 iterations
names(out3d)
out3d$zero
out3d$fnorm
out3d$niter
diff_vec(out3d$zero, wh=prob$wh, xmat=prob$x, targets=prob$targets)


library(nleqslv)


a6 <- proc.time()
out3f <- nleqslv(x=rep(0, length(prob$targets)), fn=diff_vec, jac=NULL,
                 wh=prob$wh, xmat=prob$x, targets=prob$targets,
                 method = c("Broyden"),
                 global = c("dbldog"),
                 xscalm = c("fixed"),
                 jacobian=FALSE,
                 control=list(maxit=50, trace=1, allowSingular=TRUE))
b6 <- proc.time()
b6 - a6

out3f$fvec
out3f$termcd
out3f$message
dv <- diff_vec(out3f$fvec, wh=prob$wh, xmat=prob$x, targets=prob$targets)
dv %>% round(3)
(dv / as.vector(prob$targets) * 100) %>% round(2)
sse_fn(out3f$fvec, wh=prob$wh, xmat=prob$x, targets=prob$targets) 


# nleqslv(x, fn, jac=NULL, ...,
#         method = c("Broyden", "Newton"),
#         global = c("dbldog", "pwldog",
#                    "cline", "qline", "gline", "hook", "none"),
#         xscalm = c("fixed","auto"),
#         jacobian=FALSE,
#         control = list())





# newtonsys(Ffun, x0, Jfun = NULL, ..., maxiter = 100, tol = .Machine$double.eps^(1/2))

# 5.1 hours 40k records 50 states, 13 targets per state
# secs 645.89 iter 50 rss 4208 using tpc staring point

# jacobian(diff_vec, x=as.vector(ebeta), wh=wh, xmat=xmat, targets=targets)
jfn <- function(betavec, wh, xmat, targets){
  jacobian(diff_vec, x=betavec, wh=wh, xmat=xmat, targets=targets)
}
a1 <- proc.time()
out3a <- nls.lm(par=betavec, fn = diff_vec, jac=jfn,
               control=nls.lm.control(maxiter=50, nprint=1, factor=20),
               wh=prob$wh, xmat=prob$x, targets=prob$targets)
b1 <- proc.time()
b1 - a1


# OTHER APPROACHES ----

#.. maxNR minimize sse ----
f_nlmxg <- function(betavec, wh, xmat, targets){
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  sse <- sum(d^2)
  -sse
}

pran
prob <- pran
prob <- pacs
names(prob)
opt <- maxNR(f_nlmxg, start=rep(0, prob$s * prob$k), print.level=2, wh=prob$wh, xmat=prob$x, targets=prob$targets)

opt$estimate
obeta <- matrix(opt$estimate, nrow=nrow(prob$targets), byrow=FALSE)
odelta <- get_delta(prob$wh, obeta, prob$x)
owhs <- get_weights(obeta, odelta, prob$x)
otargets <- t(owhs) %*% prob$x
otargets; prob$targets

#.. ipopt minimize sse -----

fwrap <- function(x, inputs){
  f_nlmxg <- function(betavec, wh, xmat, targets){
    
    get_delta1 <- function(wh, beta, xmat){
      beta_x <- exp(beta %*% t(xmat))
      log(wh / colSums(beta_x))
    }
    
    get_weights1 <- function(beta, delta, xmat){
      # get all weights
      beta_x <- beta %*% t(xmat)
      # add delta to every row of beta_x and transpose
      beta_xd <- apply(beta_x, 1 , function(m) m + delta) 
      exp(beta_xd)
    }
    
    beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
    delta <- get_delta1(wh, beta, xmat)
    whs <- get_weights1(beta, delta, xmat)
    etargets <- t(whs) %*% xmat
    d <- targets - etargets
    sse <- sum(d^2)
    sse
  }
  
  f_nlmxg(x, wh=inputs$wh, xmat=inputs$xmat, targets=inputs$targets)
}

gwrap <- function(x, inputs){
  grad(fwrap, x=x, method="simple", inputs=inputs) # Richardson, simple complex (danger)
}

hwrap <- function(x, obj_factor, hessian_lambda, inputs){
  # hessian(sse_fn, bvec,  wh=step_inputs$wh, xmat=step_inputs$xmat, targets=step_inputs$targets)
  obj_factor * numDeriv::hessian(fwrap, x=x, inputs=inputs) # Richardson, complex (danger)
}


pacs
prob <- scale_problem(pacs, 100)
inputs <- list()
inputs$wh <- prob$wh
inputs$xmat <- prob$x
inputs$targets <- prob$targets

x0 <- rep(0, length(inputs$targets))
x0 <- as.vector(resah$best_ebeta)

opts <- list("print_level" = 0,
             "file_print_level" = 5, # integer
             "max_iter"= 10e3,
             "linear_solver" = "ma57", # mumps pardiso ma27 ma57 ma77 ma86 ma97
             "output_file" = here::here("out", "v8.out"))

#eval_h_structure <- list()
# tmp <- hwrap(x0, obj_factor=1, hessian_lambda=1, inputs)
# str(tmp)
# eval_h_structure <- llply(1:65, .fun=function(x) 1:65)

a <- proc.time()
v1 <- ipoptr(x0=x0,
             #lb=rep(-500, length(x0)),
             #ub=rep(500, length(x0)),
             eval_f=fwrap,
             eval_grad_f=gwrap,
             #eval_h=hwrap,
             # eval_h_structure = eval_h_structure,
             opts=opts,
             inputs=inputs)
b <- proc.time()
b - a

str(v1)

bvals <- v1$solution
fbeta <- matrix(bvals, nrow=nrow(inputs$targets), byrow=FALSE)
fdelta <- get_delta(inputs$wh, fbeta, inputs$xmat)
fwhs <- get_weights(fbeta, fdelta, inputs$xmat)
round(inputs$wh - rowSums(fwhs), 2)
ftargets <- t(fwhs) %*% inputs$xmat
ftargets
inputs$targets
(ftargets / inputs$targets * 100 - 100) %>% round(2)


# nlm approach ----
f_nlmxg <- function(betavec, wh, xmat, targets){
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  sse <- sum(d^2)
  sse
}

prob <- scale_problem(pacs, 100e3)
prob$h; prob$k; prob$s
prob$targets
system.time(nlmsol2 <- nlm(f_nlmxg, rep(0, length(prob$targets)),
                         wh=prob$wh, xmat=prob$x, targets=prob$targets, iterlim = 10,
                         print.level=2))
nlmsol2 # 43 secs all vars 5k records

bvals <- nlmsol2$estimate
fbeta <- matrix(bvals, nrow=nrow(prob$targets), byrow=FALSE)
fdelta <- get_delta(prob$wh, fbeta, prob$x)
fwhs <- get_weights(fbeta, fdelta, prob$x)
round(prob$wh - rowSums(fwhs), 2)
ftargets <- t(fwhs) %*% prob$x
ftargets
prob$targets
(ftargets / prob$targets * 100 - 100) %>% round(2)


f_nlm <- function(betavec, wh, xmat, targets){
  sse_fn <- function(betavec, wh, xmat, targets){
    beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
    delta <- get_delta(wh, beta, xmat)
    whs <- get_weights(beta, delta, xmat)
    etargets <- t(whs) %*% xmat
    d <- targets - etargets
    sse <- sum(d^2)
    sse
  }
  
  
  sse <- sse_fn(betavec, wh, xmat, targets)
  g <- numDeriv::grad(sse_fn, x=betavec, wh=wh, xmat=xmat, targets=targets)
  h <- numDeriv::grad(sse_fn, x=betavec, wh=wh, xmat=xmat, targets=targets)
  
  attr(sse, "gradient") <- g
  attr(sse, "hessian") <- h
  sse
}

prob <- scale_problem(pacs, 100)
prob$h; prob$k; prob$s
prob$targets
system.time(nlmsolgh <- nlm(f_nlm, rep(0, length(prob$targets)),
                           wh=prob$wh, xmat=prob$x, targets=prob$targets, iterlim = 500,
                           print.level=2))
nlmsolgh

bvals <- nlmsol$estimate
fbeta <- matrix(bvals, nrow=nrow(prob$targets), byrow=FALSE)
fdelta <- get_delta(prob$wh, fbeta, prob$x)
fwhs <- get_weights(fbeta, fdelta, prob$x)
round(prob$wh - rowSums(fwhs), 2)
ftargets <- t(fwhs) %*% prob$x
ftargets
prob$targets
(ftargets / prob$targets * 100 - 100) %>% round(2)





system.time(res3 <- nlm(f_nlm, as.vector(beta0), wh=prob$wh, xmat=sxmat, targets=stargets, iterlim = 500))
# 283 iter = 2 mins
res3 # 462 iter, 7 mins

system.time(res4 <- mma(as.vector(beta0), fn=f_nlm, wh=hweights, xmat=sxmat, targets=stargets))

system.time(bb1 <- bobyqa(x0=as.vector(beta0), fn=f_nlm, lower = NULL, upper = NULL, nl.info = FALSE,
                          control = list(), wh=hweights, xmat=sxmat, targets=stargets))
str(bb1)

system.time(cl1 <- cobyla(x0=as.vector(beta0), fn=f_nlm, lower = NULL, upper = NULL, hin = NULL,
                          nl.info = FALSE, control = list(), wh=hweights, xmat=sxmat, targets=stargets))

system.time(slsqp(x0=as.vector(beta0), fn=f_nlm, gr = NULL, lower = NULL, upper = NULL, hin = NULL,
                  hinjac = NULL, heq = NULL, heqjac = NULL, nl.info = FALSE,
                  control = list(),  wh=hweights, xmat=sxmat, targets=stargets)) # pretty good

system.time(tnewton(x0=as.vector(beta0), fn=f_nlm,
                    control = list(),  wh=hweights, xmat=sxmat, targets=stargets))

system.time(nl <- neldermead(x0=as.vector(beta0), fn=f_nlm,
                             control = list(),  wh=hweights, xmat=sxmat, targets=stargets))

system.time(sb <- sbplx(x0=as.vector(beta0), fn=f_nlm,
                        control = list(),  wh=hweights, xmat=sxmat, targets=stargets))

library(alabama)
meths <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")
system.time(op <- optim(as.vector(beta0), fn=f_nlm, method=meths[4], wh=hweights, xmat=sxmat, targets=stargets))


# bobyqa(x0, fn, lower = NULL, upper = NULL, nl.info = FALSE,
#        control = list(), ...)

# mma(x0, fn, gr = NULL, lower = NULL, upper = NULL, hin = NULL, hinjac = NULL, nl.info = FALSE, control = list(), ...)

res <- res3x
fbeta <- matrix(res$estimate, nrow=nrow(targets), byrow=FALSE)
fdelta <- get_delta(hweights, fbeta, sxmat)
fwhs <- get_weights(fbeta, fdelta, sxmat)
round(hweights - rowSums(fwhs), 2)
ftargets <- t(fwhs) %*% xmat
ftargets
targets
(ftargets / targets * 100 - 100) %>% round(2)


D(w_h1s2, "b22") # w_h1s2 * x12
exp(b21 * x11 + b22 * x12 + c1) * x12


