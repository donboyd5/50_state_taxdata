
library(magrittr)
library(plyr) # needed for ldply; must be loaded BEFORE dplyr
library(tidyverse)
options(tibble.print_max = 65, tibble.print_min = 65) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats
library(scales)
library(hms) # hms, for times
library(lubridate) # lubridate, for date/times
library(vctrs)

library(btools) # has a few useful functions -- devtools::install_github("donboyd5/btools")

library(microweight)

library(optimx)

library(mize)

p <- make_problem(h=10, s=3, k=2)
p <- make_problem(h=100, s=10, k=5)
p <- make_problem(h=5000, s=20, k=8)


dw <- get_dweights(p$targets)
betavec <- rep(0, length(p$targets))
x <- rep(1, length(p$targets))

library(numDeriv)

system.time(hmat <- hessian(sse, betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw))
hmat
hmat2

(hmat - hmat2) %>% round(3)



# gradient of the sse function
system.time(a <- numDeriv::grad(sse, x, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw))
system.time(b <- pracma::grad(sse, x, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)) # much faster
a
b
sum(a)
sum(b)
sum((a-b)^2)

library(mize)
ff <- function(x, wh, xmat, targets, dweights=NULL){
  sse(betavec=x, wh=wh, xmat=xmat, targets=targets, dweights = dweights)
}

gf <- function(x, wh, xmat, targets, dweights=NULL){
  pracma::grad(sse, x0=x, wh=wh, xmat=xmat, targets=targets, dweights=dweights)
}

ff(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)
gf(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)

fgm <- list(fn = ff,
           gr = gf)
mize(par=betavec, fg=fgm, method = "L-BFGS", wh = p$wh, xmat = p$xmat, targets = p$targets,
     dweights = dw)

res1 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
  dweights = dw)

res2 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
  dweights = dw, method = 'Newton')

res3 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
  dweights = dw, method = 'LM')


# try something new
library(nloptr)
# allmeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", "nlminb",
#              "lbfgsb3", "Rcgmin", "Rtnmin", "Rvmmin", "spg", "ucminf",
#              "newuoa", "bobyqa", "nmkb", "hjkb", "hjn", "lbfgs", "subplex")
library(lbfgsb3c)
library(subplex)
library(ucminf)
# newuoa slow; spg slow; bobyqa slow; nmkb no progress; ucminf slow; nlm slow; Rcgmin [needs gradient]

betavec <- rep(0, length(p$targets))
a <- proc.time()
res4 <- optimx::optimr(par = betavec,
                       fn = sse,
                       method = "L-BFGS-B", # L-BFGS-B
                       control = list(trace=1,
                                      maxit=2000,
                                      kkt=FALSE,
                                      starttests=FALSE),
                       wh = p$wh, xmat = p$xmat, targets = p$targets,
                       dweights = dw)
b <- proc.time()
b - a
str(res4)


gf <- function(betavec, wh, xmat, targets, dweights=NULL){
  pracma::grad(sse, x0=betavec, wh=wh, xmat=xmat, targets=targets, dweights=dweights)
}
gf(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)

a <- proc.time()
res5 <- optimx::optimr(par = betavec,
                       fn = sse,
                       method = "L-BFGS-B", # L-BFGS-B
                       gr=gf, # "grcentral", #gf,
                       control = list(trace=1,
                                      maxit=2000,
                                      kkt=FALSE,
                                      starttests=FALSE),
                       wh = p$wh, xmat = p$xmat, targets = p$targets,
                       dweights = dw)
b <- proc.time()
b - a
str(res5)

betavec <- rep(0, length(p$targets))
a <- proc.time()
res5 <- subplex(par = betavec,
                       fn = sse,
                       control = list(maxit=200),
                       wh = p$wh, xmat = p$xmat, targets = p$targets,
                       dweights = dw)
b <- proc.time()
b - a

res4$par
as.vector(res3$beta_opt_mat)




library(dfoptim)
p
system.time(resrv <- hjk(rep(0, length(p$targets)), sse, control = list(info=TRUE), wh=p$wh, xmat=p$xmat, targets=p$targets))
system.time(resrv2 <- nmk(rep(0, length(p$targets)), sse, control = list(), wh=p$wh, xmat=p$xmat, targets=p$targets))
resrv$par
resrv2$par

library(rootSolve)
A <- matrix(nrow = 500, ncol = 500, runif(500*500))
B <- runif(500)
jfun <- function (x) A
fun <- function(x) A %*%x - B
fun <- function(x) (A %*%x - B) %>% as.vector
fun2 <- function(x, z) (A %*%x - B) %>% as.vector
system.time(X <- multiroot(start = 1:500, f = fun,
                           jactype = "fullusr", jacfunc = jfun))
system.time(X <- multiroot(start = 1:500, f = fun))

system.time(X <- multiroot(start = 1:500, f = fun, parms=list(z=10)))


assert_that(is.character(x))
assert_that(is.matrix(A))

f <- function(x, z){
  c(x[1]+10, x[2]-5)^2 * z
}
f(c(2, 3), 10)
multiroot(start = c(2, 3), f, z=10)

system.time(X <- multiroot(start = 1:500, f = fun))

sum(fun(X$root)^2)
fun(1:500)

p2 <- make_problem(h=4000, s=30, k=15)

system.time(v1 <- multiroot(start=rep(0, length(p2$targets)), f=diff_vec,
                            atol=1e-4,
                            wh=p2$wh, xmat=p2$xmat, targets=p2$targets))
max(abs(v1$f.root))

v2 <- geoweight(wh = p2$wh, xmat = p2$xmat, targets = p2$targets,
          method = 'Newton')
max(abs(v2$output$fvec))
v2$etime

p2$targets
targets_mat(as.vector(v2$beta_opt_mat), p2$wh, p2$xmat, p2$s)
targets_mat(as.vector(v1$root), p2$wh, p2$xmat, p2$s)



v2 <- geoweight(wh = p2$wh, xmat = p2$xmat, targets = p2$targets,
                dweights = rep(1, length(p2$targets)), method = 'Newton')

v3 <- multiroot(diff_vec, rep(0, length(p$targets)), maxiter = 100,
                wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=get_dweights(p$targets))

v4 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
                dweights = rep(1, length(p$targets)), method = 'Newton')

v4 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets, method = 'Newton')


p2$targets
targets_mat(as.vector(v2$beta_opt_mat), p2$wh, p2$xmat, p2$s)
targets_mat(as.vector(v1$root), p2$wh, p2$xmat, p2$s)


multiroot(diff_vec, rep(0, length(p$targets)), maxiter = 100,
          rtol = 1e-6, atol = 1e-8, ctol = 1e-8,
          useFortran = TRUE, positive = FALSE,
          jacfunc = NULL, jactype = "fullint",
          verbose = FALSE, bandup = 1, banddown = 1,
          parms = NULL, wh=p$wh, xmat=p$xmat, targets=p$targets)

rosenbrock <- function(x){
  n <- length(x)
  sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}
par0 <- rep(0, 10)
hjk(par0, rosenbrock)

hjkb(c(0, 0, 0), rosenbrock, upper = 0.5)

library(magrittr)
library(plyr) # needed for ldply; must be loaded BEFORE dplyr
library(tidyverse)
options(tibble.print_max = 65, tibble.print_min = 65) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats
library(scales)
library(hms) # hms, for times
library(lubridate) # lubridate, for date/times
library(vctrs)

library(btools) # has a few useful functions -- devtools::install_github("donboyd5/btools")

# library(readxl)
# library(grDevices)
# library(knitr)
# library(kableExtra)
# library(janitor)
# library(ipoptr)


glimpse(acs)
count(acs, stabbr)
count(acs, incgroup)

acsdf <- acs
acsdf <- acsbig
glimpse(acsdf)

acs_prep <- acsdf %>%
  filter(incgroup == 2) %>%
  mutate(pop = 1,
         young = agep < 21,
         workage = agep %in% 22:64,
         married = mar==1,
         female = sex==1,
         ssip_nnz = ssip != 0) %>%
  # convert to indicator variables
  mutate(across(c(young, workage, married, female, ssip_nnz), as.integer)) %>%
  select(serialno, stabbr, pwgtp,
         pop, pincp, wagp, ssip, ssip_nnz, young, workage, married, female)
glimpse(acs_prep)

xmat <- acs_prep %>% select(-serialno, -stabbr, -pwgtp) %>%
  as.matrix()

wh <- acs_prep$pwgtp

targs_df <- acs_prep %>%
  group_by(stabbr) %>%
  summarise(across(pop:female, ~ sum(. * pwgtp)), .groups = "drop")
targets <- as.matrix(targs_df[ , -1])
rownames(targets) <- targs_df$stabbr

identical(colnames(xmat), colnames(targets))

colnames(xmat) <- NULL
colnames(targets) <- NULL

p <- list()
p$h <- nrow(xmat)
p$s <- nrow(targets)
p$k <- ncol(targets)
p$wh <- acs_prep$pwgtp
p$xmat <- xmat

p$h; p$s; p$k

set.seed(1234)
noise <- rnorm(n=length(targets), mean=0, sd=.01)
# noise <- 0
p$targets <- targets * (1 + noise)

res1 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'Broyden')
res2 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'Newton')
res3 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'LM')

rbind(res1$solver_message, res2$solver_message, res3$solver_message)
cbind(res1$sse_unweighted, res2$sse_unweighted, res3$sse_unweighted)
cbind(res1$sse_weighted, res2$sse_weighted, res3$sse_weighted)

res1$targets_diff %>% round(1)
res2$targets_diff %>% round(1)
res3$targets_diff %>% round(1)

res1$targets_pctdiff %>% round(1)
res2$targets_pctdiff %>% round(1)
res3$targets_pctdiff %>% round(1)
c(max(abs(res1$targets_pctdiff)),
  max(abs(res2$targets_pctdiff)),
  max(abs(res3$targets_pctdiff))) %>% round(3)

cbind(res1$etime, res2$etime, res3$etime)

quantile(res3$beta_opt_mat)


spoint <- get_starting_point(p)
spoint$result$solver_message
spoint$result$etime

res1s <- geoweight(wh=p2$wh, xmat=p2$xmat, targets=p2$targets, method = 'Broyden')
res1s$sse_weighted
res1s$sse_unweighted
res1s$etime
res1$beta_opt_mat
res1s$beta_opt_mat

res2s <- geoweight(wh=p2$wh, xmat=p2$xmat, targets=p2$targets, method = 'Newton')
res2s$sse_weighted
res2s$sse_unweighted
res2s$etime
res2$beta_opt_mat
res2s$beta_opt_mat

res3s <- geoweight(wh=p2$wh, xmat=p2$xmat, targets=p2$targets, method = 'LM')
res3s$sse_weighted
res3s$sse_unweighted
res3s$etime
res3$beta_opt_mat
res3s$beta_opt_mat

res2_2step <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'Newton', betavec=spoint$spoint)
res2_2step$etime

res3_2step <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'LM', betavec=spoint$spoint)
res3_2step$etime

rbind(res1$solver_message, res2$solver_message, res3$solver_message,
      res2_2step$solver_message, res3_2step$solver_message)

cbind(res1$etime, res2$etime, res3$etime, res2_2step$etime, res3_2step$etime)

cbind(res1$sse_unweighted, res2$sse_unweighted, res3$sse_unweighted,
      res2_2step$sse_unweighted, res3_2step$sse_unweighted)

cbind(res1$sse_weighted, res2$sse_weighted, res3$sse_weighted,
      res2_2step$sse_weighted, res3_2step$sse_weighted)


res1$targets_diff %>% round(1)
res2$targets_diff %>% round(1)
res3$targets_diff %>% round(1)
res2_2step$targets_diff %>% round(1)
res3_2step$targets_diff %>% round(1)

res1$targets_pctdiff %>% round(1)
res2$targets_pctdiff %>% round(1)
res3$targets_pctdiff %>% round(1)
res3_2step$targets_pctdiff %>% round(1)

c(max(abs(res1$targets_pctdiff)),
  max(abs(res2$targets_pctdiff)),
  max(abs(res3$targets_pctdiff)),
  max(abs(res2_2step$targets_pctdiff)),
  max(abs(res3_2step$targets_pctdiff))
  ) %>% round(3)


res3$beta_opt_mat %>% round(4)
res2_2step$beta_opt_mat %>% round(4)
res3_2step$beta_opt_mat %>% round(4)

summary(res3$whs)
summary(res2_2step$whs)

# now do tax data ----



# could we use ipopt?? how fast is gradient?? ----
dw <- get_dweights(p$targets)
betavec <- rep(0, length(p$targets))
# gradient of the sse function
# system.time(a <- numDeriv::grad(sse, x, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw))
system.time(b <- pracma::grad(sse, betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)) # much faster

gf <- function(x, wh, xmat, targets, dweights=NULL){
  pracma::grad(sse, x0=x, wh=wh, xmat=xmat, targets=targets, dweights=dweights)
}

inputs <- list()
inputs$wh <- p$wh
inputs$xmat <- p$xmat
inputs$targets <- p$targets
inputs$dweights <- dw

library(ipoptr)

opts_ipopt <- list("print_level" = 0,
                   "file_print_level" = 5, # integer
                   "print_user_options" = "yes",
                   "max_iter"= 50,
                   "linear_solver" = "ma86", # mumps pardiso ma27 ma57 ma77 ma86 ma97
                   # "mehrotra_algorithm" = "yes", # default no
                   "hessian_approximation" = "exact", # default exact; limited-memory
                   # "limited_memory_update_type" = "sr1", # default bfgs; sr1 docs say "not working well" but this really helps
                   "obj_scaling_factor" = 3, # 1e-3, # default 1; 1e-1 pretty fast to feasible but not to optimal
                   # "nlp_scaling_max_gradient" = 100, # default is 100 - seems good
                   # "accept_every_trial_step" = "yes",
                   # "derivative_test" = "first-order",
                   "output_file" = here::here("temp3.out"))

set.seed(1)
xran <- runif(length(betavec), -10, 10)

# maybe fast hessian?

ip1 <- ipoptr(x0 = betavec,
              eval_f = sse, # arguments: x, inputs; eval_f_xtop eval_f_xm1sq eval_f_wfs
              eval_grad_f = gf, # eval_grad_f_xm1sq, # eval_grad_f_xtop eval_grad_f_xm1sq
              lb=rep(-100, length(betavec)),
              ub=rep(100, length(betavec)),
              # eval_h = eval_h_xm1sq, # the hessian is essential for this problem eval_h_xtop eval_h_xm1sq
              # eval_h_structure = inputs$eval_h_structure,
              opts = opts_ipopt,
              wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)


x2 <- ip1$solution
# res3$beta_opt_mat

iptarg <-  targets_mat(x2, p$wh, p$xmat, p$s)
iptarg


fn <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/data/acs_10krecs_5states.rds"
df <- readRDS(fn)
acs <- df %>%
  mutate(incgroup=ntile(pincp, 10)) %>%
  select(-pop, -nrecs)
usethis::use_data(acs)
quantile(df$pincp)


library(minpack.lm) # nls.lm
library(nleqslv) # nleqslv

p <- make_problem(h=10, s=50, k=10)

p$targets
dw <- get_dweights(p$targets)
betavec <- rep(0, length(p$targets))

v1 <- geoweight(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
        method="Newton")

v2 <- geoweight(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
        method="Broyden")

v3 <- geoweight(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
        method="levmarq")

v1$etime
v2$etime
v3$etime

(diff_vec(v1$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets) / as.vector(p$targets) * 100) %>% round(3)
(diff_vec(v2$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets) / as.vector(p$targets) * 100) %>% round(3)
(diff_vec(v3$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets) / as.vector(p$targets) * 100) %>% round(3)

diff_vec(v1$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets) %>% round(3)
diff_vec(v2$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets) %>% round(3)
diff_vec(v3$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets) %>% round(3)

sse(v1$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets)
sse(v2$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets)
sse(v3$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets)

dt <- t4 - t3
dt

cbind(v1$x, v2$x, v3$x)

first <- list(a="1", b="2")
second <- list(b="3", d="4")

Map(c, first, second)


p <- make_problem(h=1000, s=20, k=10)
dw <- get_dweights(p$targets)

res1 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
                  dweights = get_dweights(p$targets), maxiter=20)

res2 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
                  dweights = get_dweights(p$targets), method = 'Newton', maxiter=10)

res3 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
   dweights = get_dweights(p$targets), method = 'LM', maxiter=20, opts=list(factor=10))


res1
res2
res3
c(res1$sse_unweighted, res2$sse_unweighted, res3$sse_unweighted)
c(res1$sse_weighted, res2$sse_weighted, res3$sse_weighted)

if (!assertthat::is.number(n)) { # make sure n is just one value and numeric
  stop("n should be a number of length 1!")
}
assert_that(is.character(x))

start <- rep(0, length(p$targets))
a <- proc.time()
res4 <- multiroot(start=start,
                  f=diff_vec, verbose=TRUE,
                  wh=p$wh, xmat=p$xmat, targets=p$targets)
b <- proc.time()
b - a
res4

library(BB)
start <- rep(0, length(p$targets))
t1 <- proc.time()
res5 <- BBoptim(par=start, fn=sse,
                method=3,
                wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=get_dweights(p$targets))
t2 <- proc.time()
t2 - t1

res6 <- spg(start, fn=sse, control=list(), quiet=FALSE, alertConvergence=TRUE,
            wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=get_dweights(p$targets))

res1$beta_opt_mat

res1$targets_calc
res2$targets_calc
res3$targets_calc
targets_mat(res5$par, p$wh, p$xmat, p$s)
targets_mat(res6$par, p$wh, p$xmat, p$s)
targets_mat(as.vector(res3$beta_opt_mat), p$wh, p$xmat, p$s)
p$targets

res1$etime; res2$etime; res3$etime; (b - a)
ssew <- sse(betavec=res4$root, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=get_dweights(p$targets))


res4 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=rep(1, length(targets)))
res5 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = "Newton", dweights=rep(1, length(targets)))
res6 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = "LM", dweights=rep(1, length(targets)))

cbind(res1$sse_weighted, res2$sse_weighted, res3$sse_weighted)
cbind(res1$sse_unweighted, res2$sse_unweighted, res3$sse_unweighted)

res1$etime; res2$etime; res3$etime



# res2d <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=rep(1, length(targets)), method = "Newton")

out <- res1$targets_calc
out <- res2$targets_calc
# out <- res2d$targets_calc
out <- res3$targets_calc

out <- res4$targets_calc
out <- res5$targets_calc
out <- res6$targets_calc


out
p$targets

(out - p$targets) %>% round(4)

(out / p$targets * 100 - 100) %>% round(4)


res1$etime
res2$etime
res3$etime

res2$whs %>% round(2) %>% head(10)
whs2 %>% round(2) %>% head(10)

st <- "CA"
st <- "FL"
cor(whs2[, st], res2$whs[, st])


library(purrr)
x <- list(x = 1:10, y = 4, z = list(a = 1, b = 2))
str(x)
# Update values
str(list_modify(x, a = 1))
x
str(list_modify(x, z = 5))

a <- list(x = 1:10, y = 4)
str(a)
b <- list(y=7, m=1:3)
str(b)
str(list_modify(a, !!!b))

str(list_merge(a, b))
str(update_list(a, b))


df <- tibble(
  x = 1:3,
  y = c("a", "d,e,f", "g,h")
)
df %>%
  transform(y = strsplit(y, ",")) %>%
  unnest(y)
