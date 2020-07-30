
# code folding ----
# alt-o, shift-alt-o
# alt-l, shift-alt-l
# alt-r

# notes ----


# libraries ----
source(here::here("include", "libraries.r"))
# remotes::install_github("tidyverse/dplyr") if needed
library(numDeriv)
library(ipoptr)
library(nloptr)

devtools::session_info()
(.packages()) %>% sort

# globals ----
dbox <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/"
(fns <- paste0(c("acs_10krecs_5states", "acs_100krecs_20states", "acs_200krecs_50states", "acs_400krecs_50states"), ".rds"))

# functions ----
source(here::here("include", "functions_prep_dev.r")) # soon we will replace functions_prep.r with the dev version
source(here::here("include", "functions_poisson_model2.r"))


#.. poisson-related functions ----




#..functions to create problem of chosen size ----

make_problem <- function(h, k, s){
  # h: # of households
  # k: # of characteristics per household
  # s: # of states
  
  # returns a list with items created below
  
  # example call: make_problem(8, 2, 3)
  
  set.seed(1234)
  x <- matrix(runif(h * k), nrow=h, byrow=TRUE)
  
  set.seed(1234)
  whs <- matrix(runif(h * s, 10, 20), nrow=h, byrow=TRUE)
  
  wh=rowSums(whs)
  ws=colSums(whs)
  
  targets <- t(whs) %*% x # s x k
  
  keepnames <- c("h", "k", "s", "x", "whs", "wh", "ws", "targets")
  problem <- list()
  for(var in keepnames) problem[[var]] <- get(var)
  problem
}


# choose which file to use ----
samp1 <- readRDS(here::here("data", fns[3])) %>% 
  select(-nrecs, -pop) # note that we no longer need nrecs; pop ordinarily would not be in the data so drop here and create later
glimpse(samp1)
summary(samp1)
count(samp1, mar)
# djb note that we have nrecs and pop variables -- I want to make sure they are not needed for anything ----
# if you want to target the total number of weighted records we need a variable that is 1 for all records ----

# PREPARE DATA  ----
#.. modify the sample (don't think we need a function for this) ----
# - define income groups
# - create an indicator for each income variable as to whether it is nonzero
# _ expand categoricals into dummies as needed
# if we don't have a variable such as pop where all values are 1, we should create it as it makes it easy to get weighted record counts
samp2 <- samp1 %>%
  mutate(pid=row_number(), # pid -- an id variable for each person in the file
         incgroup=ntile(pincp, 10), # divide the data into 10 income ranges
         pop=1, # it's useful to have a variable that is 1 on every record
         # convert categoricals to dummies if we will base targets upon them
         mar1=ifelse(mar==1, 1, 0), # married
         mar5=ifelse(mar==5, 1, 0), # single
         marx15=ifelse(mar %nin% c(1, 5), 1, 0)
  )
summary(samp2)
ht(samp2)

#.. define the kinds of (weighted) targets we want and prepare the file accordingly ----
# sum:    sum of values
# nnz:    number of nonzero values
# sumneg: sum of negative values
# nneg:   number of zero values
# sumpos: sum of positive value
# npos:   number of positive values

# For the PUF the SOI data provide only the first two kinds of targets, but for the ACS we could have any of them.

# TRY TO AVOID DEPENDENT CONSTRAINTS - redundancy - as they can make the problem very hard to solve.
# For example, suppose there are 3 kinds of income (wages, interest, retirement) plus a total (sum of the 3)
#   -- don't create targets for each of the 3 kinds plus a target for the total -- leave one item out
# Another, less obvious example: don't target the total number of returns plus the number for each marital status - leave one out.

nnz_vars <- c("pop", "mar1", "mar5", "pincp", "wagp") # note that I leave the 3rd marital status out -- marx15
sum_vars <- c("pincp", "wagp", "intp", "pap", "retp", "ssip", "ssp") # DO NOT include total plus all sums - leave one out (otherincp)
sumneg_vars <- "otherincp"

# define a vector of variable names for "possible" targets (a superset) -- we may not target all
possible_target_vars <- make_target_names(
  list(nnz_vars, sum_vars, sumneg_vars),
  c("nnz", "sum", "sumneg"))
possible_target_vars

# prepare data by creating variables with those names:
#   nnz, nneg, and npos will have 1 for rows where the respective variable is nz, neg, or pos, respectively, and 0 otherwise
#   sum will have its respective variable's value
#   sumneg and sumpos will have the variable's value if negative or positive, respectively, and 0 otherwise
samp <- prep_data(samp2, possible_target_vars)
glimpse(samp)


summary_vals <- get_summary_vals(samp, .weight=pwgtp, .sum_vars=possible_target_vars, stabbr, incgroup)
summary_vals

# Create a data frame with all targets for all states and income groups ----
# for the PUF, we will create this using information from Historical Table 2
# for the ACS, we construct the targets from the ACS data itself
all_target_values <- summary_vals


# wrap everything we need for a single income group into a function that returns a list ----

# SINGLE INCOME GROUP ----
#.. define target incgroup, target variable names, and target values for each state in the income group ----
target_incgroup <- 2 # define target income group

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

target_vars

# define target values and states, for this income group
targets_wide <- all_target_values %>%
  filter(incgroup==target_incgroup) %>%
  select(stabbr, incgroup, nrecs, all_of(target_vars)) # a small list of variables to target; we have nrecs because we created it in summary_vals
targets_wide # these are the targets we want to hit
summary(targets_wide)

hweights <- samp %>% filter(incgroup==target_incgroup) %>% .$pwgtp
targets <- targets_wide[, target_vars] %>% as.matrix
xmat <- samp %>% filter(incgroup==target_incgroup) %>% .[, target_vars] %>% as.matrix


targets
xmat

pacs <- list()
pacs$h <- nrow(xmat)
pacs$k <- ncol(xmat)
pacs$s <- nrow(targets)
pacs$x <- xmat
pacs$wh <- hweights
pacs$targets <- targets
pacs


# SOLVE in loop ----

p <- make_problem(h=8, k=2, s=3)
p <- make_problem(h=100, k=5, s=10)

res <- solve_poisson(p, step_scale = 100)
res <- solve_poisson(p, scale=TRUE, scale_goal = 100, step_scale=7)

#.. acs 1 ----
res <- solve_poisson(pacs, step_scale=2e3) # acs 1, 4 targets
res <- solve_poisson(pacs, step_scale=.75e3, scale=TRUE, scale_goal=100) # acs 1, 13 targets
res <- solve_poisson(pacs, step_scale=1000, scale=TRUE, scale_goal=100, maxiter = 200) # acs 1, 13 targets

#.. acs 2 ----
res <- solve_poisson(pacs, step_scale=6e3) # acs 2, 4 targets
system.time(res <- solve_poisson(pacs, step_scale=8000, maxiter=20000)) # acs 2, 13 targets
system.time(res <- solve_poisson(pacs, step_scale=8000, maxiter=20000, scale=TRUE, scale_goal=100)) # acs 2, 13 targets

#.. acs 3 ----
res <- solve_poisson(pacs, step_scale=11e3) # acs 3, 4 targets
res <- solve_poisson(pacs, step_scale=11e3, scale=TRUE, scale_goal=100, maxiter=20000) # acs 3, 13 targets


# examine results ----
result <- res
names(result)
# str(result)
result$iter
result$sse
result$sse_vec
result$d %>% round(4)

# check weights
(result$ewh - result$problem_unscaled$wh) %>% round(2) # total household weights

# check targets
result$etargets
result$problem_unscaled$targets
(result$etargets / result$problem_unscaled$targets * 100 - 100) %>% round(2)



# OTHER APPROACHES ----
# nlm approach ----
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
  g <- grad(sse_fn, x=betavec, wh=wh, xmat=xmat, targets=targets)
  
  attr(sse, "gradient") <- g
  sse
}


f_nlmxg <- function(betavec, wh, xmat, targets){
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  sse <- sum(d^2)
  sse
}

g <- grad(sse_fn, x=betavec, wh=wh, xmat=xmat, targets=targets)



system.time(res3 <- nlm(f_nlm, as.vector(beta0), wh=hweights, xmat=sxmat, targets=stargets, iterlim = 500))
# 283 iter = 2 mins
res3 # 462 iter, 7 mins

system.time(res3x <- nlm(f_nlmxg, as.vector(beta0), wh=hweights, xmat=sxmat, targets=stargets, iterlim = 500))
res3x # 43 secs all vars 5k records



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


# try it with ipopt -----

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


inputs <- list()
inputs$wh <- hweights
inputs$xmat <- sxmat
inputs$targets <- stargets

x0 <- as.vector(beta0)

opts <- list("print_level" = 0,
             "file_print_level" = 5, # integer
             "max_iter"= 500,
             "linear_solver" = "ma57", # mumps pardiso ma27 ma57 ma77 ma86 ma97
             "output_file" = here::here("out", "v5.out"))

a <- proc.time()
v1 <- ipoptr(x0=x0,
             lb=rep(-500, length(x0)),
             ub=rep(500, length(x0)),
             eval_f=fwrap,
             eval_grad_f=gwrap,
             opts=opts,
             inputs=inputs)
b <- proc.time()
b - a

str(v1)

bvals <- v1$solution
fbeta <- matrix(bvals, nrow=nrow(targets), byrow=FALSE)
fdelta <- get_delta(hweights, fbeta, sxmat)
fwhs <- get_weights(fbeta, fdelta, sxmat)
round(hweights - rowSums(fwhs), 2)
ftargets <- t(fwhs) %*% xmat
ftargets
targets
(ftargets / targets * 100 - 100) %>% round(2)




opts <- list("print_level" = 0,
             "file_print_level" = 5, # integer
             "max_iter"= 100,
             "linear_solver" = "ma57", # mumps pardiso ma27 ma57 ma77 ma86 ma97
             # "ma57_automatic_scaling" = "yes", # if using ma57
             # "ma57_pre_alloc" = 3, # 1.05 is default; even changed, cannot allocate enough memory, however
             # "ma77_order" = "amd",  # metis; amd -- not clear which is faster
             "mehrotra_algorithm" = "yes",
             "obj_scaling_factor" = 1, # 1e-3, # default 1; 1e-1 pretty fast to feasible but not to optimal
             # "nlp_scaling_method" = "equilibration-based", # NO - use default gradient_based
             "nlp_scaling_max_gradient" = 100, # default is 100 - seems good 
             # "jac_c_constant" = "yes", # does not improve on moderate problems
             # "jac_d_constant" = "yes", # does not improve on  moderate problems
             # "hessian_constant" = "yes", # KEEP default NO - if yes Ipopt asks for Hessian of Lagrangian function only once and reuses; default "no"
             # "hessian_approximation" = "limited-memory", # KEEP default of exact
             # "derivative_test" = "first-order",
             # "derivative_test_print_all" = "yes",
             "output_file" = here::here("out", "test57.out"))

setwd(here::here("temp1"))
getwd()
result <- ipoptr(x0 = inputs$x0,
                 lb = inputs$xlb,
                 ub = inputs$xub,
                 eval_f = eval_f_xm1sq, # arguments: x, inputs; eval_f_xtop eval_f_xm1sq
                 eval_grad_f = eval_grad_f_xm1sq, # eval_grad_f_xtop eval_grad_f_xm1sq
                 eval_g = eval_g, # constraints LHS - a vector of values
                 eval_jac_g = eval_jac_g,
                 eval_jac_g_structure = inputs$eval_jac_g_structure,
                 eval_h = eval_h_xm1sq, # the hessian is essential for this problem eval_h_xtop eval_h_xm1sq
                 eval_h_structure = inputs$eval_h_structure,
                 constraint_lb = inputs$clb,
                 constraint_ub = inputs$cub,
                 opts = opts,
                 inputs = inputs)




d <- 1:5


D(w_h1s2, "b22") # w_h1s2 * x12
exp(b21 * x11 + b22 * x12 + c1) * x12


