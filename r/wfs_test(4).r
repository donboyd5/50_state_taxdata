source(here::here("include", "functions_ipopt.r"))
library(numDeriv)

# targets_wide <- targets_wide %>%
#   rename(stabbr=STATE)
  # rename(STATE=stabbr)

# djb start ----

targets_df2 <- get_targets2(targets_wide, target_vars, .pid=incgroup_data$pid, .weight_total=incgroup_data$weight_total)
ht(targets_df2)

# targets_df <- targets_df %>%
#   filter(targtype != "addup")

# targets_df %>%
#   filter(targtype!="aggregate")

# CONSTRUCT initial weights ----
#.. prepare the data we will reweight to hit (or come close to) the targets ----
# let's stack the data so that each person appears 1 time for each state that is targeted
# create a stub with a record for every person and every targeted state

# define initial weight for each person-state combination. This is important as the
# optimization will try to keep the final weights near these initial weights. The
# more we can improve these initial weights, the better the results are likely to be
iweights <- get_initial_weights(targets_wide, incgroup_data, .popvar=pop_nnz)
head(iweights)


# DEFINE constraint coefficients (cc) ----
# cc_sparse <- get_constraint_coefficients(incgroup_data, target_vars, iweights, targets_df) # can take a while on big problems


cc_sparse2 <- get_constraint_coefficients2(.incgroup_data=incgroup_data,
                                           .target_vars=target_vars,
                                           .iweights=iweights,
                                           .targets_df=targets_df2) # can take a while on big problems
count(cc_sparse2, targtype)


inputs2 <- get_inputs(.targets_df=targets_df2, .iweights=iweights, .cc_sparse=cc_sparse2,
                      .targtol=.05, .adduptol=0, .xub=50, .conscaling=TRUE, scale_goal=10)
iweights

# make what we need for adding-up constraints
cc_addup_sparse <- inputs2$cc_sparse %>%
  select(j, pid, stabbr) %>%
  distinct() %>%
  group_by(pid) %>%
  mutate(ipid=cur_group_id()) %>%
  ungroup %>%
  arrange(j) %>%
  left_join(iweights %>% select(pid, stabbr, weight_total, iweight_state), by = c("pid", "stabbr"))

weight_totals <- cc_addup_sparse %>%
  select(ipid, weight_total) %>%
  distinct %>%
  arrange(ipid)

inputs2$cc_addup_sparse <- cc_addup_sparse
inputs2$weight_totals <- weight_totals
inputs2$eval_jac_g_structure <- define_jac_g_structure_sparse(inputs2$cc_addup_sparse, ivar="ipid", jvar="j")



# also make the dense jacobian - DANGER - costly in time ----
# jacg_dense <- matrix(0, nrow=length(inputs2$eval_jac_g_structure), ncol = inputs2$n_variables)
# f <- function(rownum){
#   colnums <- inputs2$eval_jac_g_structure[[rownum]]
#   vals <- inputs2$cc_addup_sparse$iweight_state[colnums]
#   jacg_dense[rownum, colnums] <<- vals
# }
# l_ply(1:length(inputs2$eval_jac_g_structure), f)
# sum(jacg_dense)
# dim(jacg_dense)
# jacg_dense[1:10, 1:30]
# inputs2$jacg_dense <- jacg_dense
# jacg_dense[20, 1:100]
# 
# quantile(jacg_dense)
# str(jacg_dense)

names(inputs2)
inputs2$n_variables
inputs2$n_constraints
inputs2$n_targets
inputs2$objscale

eval_f_wfs(inputs2$x0, inputs2)
eval_g_addup(inputs2$x0, inputs2)

# djb jacobian checks ----
# jcheck <- eval_jac_g_addup(inputs2$x0, inputs2)
# jcheck[1:12]
# 
# inputs2$eval_jac_g_structure[1:10]

# system.time(jfd <- jacobian(eval_g_addup, inputs2$x0, inputs=inputs2))
# jfd[1:10, 1:12]
# 
# jrows <- 30
# jcols <- inputs2$eval_jac_g_structure[[jrows]]
# start <- (jrows - 1)*5 + 1
# jcheck[start:(start+4)]
# jfd[jrows, jcols]

# end jacobian checks ----


# inputs2$objscale <- 1e4

# a <- proc.time()
# der <- check.derivatives(inputs2$x0, 
#                          func=eval_f_wfs, 
#                          func_grad=eval_grad_f_wfs, 
#                          inputs=inputs2,
#                          check_derivatives_print = "errors")
# b <- proc.time()
# b - a


# system.time(tmp <- check.derivatives(.x=inputs2$x0, func=eval_f_wfs, func_grad=eval_grad_f_wfs,
#                          check_derivatives_print='none', inputs=inputs2))
# str(tmp)
# tb <- tibble(ad=tmp$analytic, fd=tmp$finite_difference, re=tmp$relative_error, flag=tmp$flag_derivative_warning, i=1:length(flag))
# ht(tb)
# tb %>% filter(flag)

xlb <- rep(0, inputs2$n_variables)
xub <- rep(50, inputs2$n_variables)
inputs2$objscale <- 1
a <- proc.time()
sse_opt <- optim(par=sse_opt$par, # inputs2$x0 sse_opt$par
              fn = eval_f_wfs,
              gr = eval_grad_f_wfs,
              inputs=inputs2,
              lower = xlb, upper = xub,
              method = "L-BFGS-B",
              control = list(trace=1, maxit=100))
b <- proc.time()
b - a

eval_f_wfs(inputs2$x0, inputs2)
sse_opt$value

cw <- check_weights(sse_opt$par)
cw %>% ht
sse_weights(cw)
cw %>%
  ggplot(aes(weight_total, weight)) +
  geom_point()

# sse_opt2 ----
library(optimr) # seems older than optimx
library(optimx) # this seems newer than optimr (??)

inputs2$objscale <- 1

xlb <- rep(0, inputs2$n_variables)
xub <- rep(50, inputs2$n_variables)
a <- proc.time()
sse_opt2 <- optimx::optimr(par=inputs2$x0, # inputs2$x0  sse_opt2$par ip1a$solution
                 fn = eval_f_wfs,
                 gr = eval_grad_f_wfs,
                 lower = xlb, upper = xub,
                 method = "L-BFGS-B", # L-BFGS-B Rcgmin BFGS
                 control = list(trace=1, 
                                maxit=2000, 
                                kkt=FALSE,
                                starttests=FALSE),
                inputs=inputs2)
b <- proc.time()
b - a
str(sse_opt2)
# par_bak <- sse_opt2$par


cw <- check_weights(sse_opt2$par)
cw %>% ht
sse_weights(cw)
cw %>%
  ggplot(aes(weight_total, weight)) +
  geom_point()

# how did we do against constraints?
wadjust <- inputs2$cc_addup_sparse %>%
  mutate(x=sse_opt2$par,
         wstate=iweight_state * x) %>%
  group_by(ipid, pid) %>%
  mutate(wstate_adj=wstate / sum(wstate) * weight_total,
         xadj=wstate_adj / iweight_state) %>%
  ungroup

eval_f_wfs(sse_opt2$par, inputs2)
# (v1 <- eval_f_wfs(wadjust$xadj, inputs2))
# v1
eval_f_wfs(wadjust$xadj, inputs2)

cw <- check_weights(wadjust$xadj)
cw %>% ht
sse_weights(cw)
cw %>%
  ggplot(aes(weight_total, weight)) +
  geom_point()


scalef <- inputs2$cc_sparse %>% select(cname, scale) %>% distinct
  
unscaled <- function(x, scale) x * scale
condf <- tibble(cname=inputs2$constraint_names,
                con=inputs2$constraints) %>%
  mutate(conx0=eval_f_obj_elements(inputs2$x0, inputs2),
         conx=eval_f_obj_elements(sse_opt2$par, inputs2),
         conx_adj=eval_f_obj_elements(wadjust$xadj, inputs2)) %>%
  left_join(scalef, by = "cname") %>%
  # mutate(scale=1) %>%
  mutate(across(c(con, conx0, conx, conx_adj), ~unscaled(.x, scale), .names = "{col}_raw")) %>%
  mutate(diff=conx_adj_raw - con_raw,
         pdiff=diff / con_raw * 100)

condf %>%
  select(cname, contains("raw"), diff, pdiff) %>%
  arrange(desc(abs(pdiff)))

condf %>%
  filter(str_detect(cname, "pap")) %>%
  arrange(desc(abs(pdiff)))


# ipopt ----

inputs2$objscale <- 1
eval_f_wfs(inputs2$x0, inputs2)
# eval_f_wfs(sse_opt2$par, inputs2)

opts_ipopt <- list("print_level" = 0,
                   "file_print_level" = 5, # integer
                   "max_iter"= 2000,
                   "linear_solver" = "ma86", # mumps pardiso ma27 ma57 ma77 ma86 ma97
                   # "mehrotra_algorithm" = "yes", # default no
                   "jac_c_constant" = "yes", # default no
                   "hessian_approximation" = "limited-memory", # default exact; limited-memory
                   "limited_memory_update_type" = "sr1", # default bfgs; sr1 docs say "not working well" but this really helps
                   # "hessian_constant" = "yes", # default no; yes
                   "obj_scaling_factor" = 1e7, # 1e-3, # default 1; 1e-1 pretty fast to feasible but not to optimal
                   "nlp_scaling_max_gradient" = 100, # default is 100 - seems good
                   "output_file" = here::here("out", "fn4t13_os1e7_mos1_ma86_xsse.out"))

# setwd(here::here("temp1"))
#getwd()
ip1 <- ipoptr(x0 = wadjust$xadj, # inputs2$x0 wadjust$xadj ip1$solution wfs1$solution, sse_opt$par sse_opt2$par
              lb = inputs2$xlb,
              ub = inputs2$xub,
              eval_f = eval_f_wfs, # arguments: x, inputs; eval_f_xtop eval_f_xm1sq eval_f_wfs
              eval_grad_f = eval_grad_f_wfs, # eval_grad_f_xm1sq, # eval_grad_f_xtop eval_grad_f_xm1sq
              eval_g = eval_g_addup, # constraints LHS - a vector of values
              eval_jac_g = eval_jac_g_addup,
              eval_jac_g_structure = inputs2$eval_jac_g_structure,
              # eval_h = eval_h_xm1sq, # the hessian is essential for this problem eval_h_xtop eval_h_xm1sq
              # eval_h_structure = inputs$eval_h_structure,
              constraint_lb = inputs2$weight_totals$weight_total,
              constraint_ub = inputs2$weight_totals$weight_total,
              opts = opts_ipopt,
              inputs = inputs2)


# length(eval_jac_g_addup(inputs2$x0, inputs2))
# length(unlist(inputs2$eval_jac_g_structure))
ip1a <- ip1

eval_f_wfs(inputs2$x0, inputs2)
eval_f_wfs(ip1a$solution, inputs2)

cw <- check_weights(ip1a$solution)
cw %>% ht
sse_weights(cw)
cw %>%
  ggplot(aes(weight_total, weight)) +
  geom_point()


wcheck <- inputs2$weight_totals %>%
  mutate(#weight1=eval_g_addup(wfs1$solution, inputs2),
    weight2=eval_g_addup(ip1a$solution, inputs2))
wcheck

# how did we do against constraints?
# wadjust <- inputs2$cc_addup_sparse %>%
#   mutate(x=ip1$solution,
#          wstate=iweight_state * x) %>%
#   group_by(ipid, pid) %>%
#   mutate(wstate_adj=wstate / sum(wstate) * weight_total,
#          xadj=wstate_adj / iweight_state) %>%
#   ungroup


scalef <- inputs2$cc_sparse %>% select(cname, scale) %>% distinct

unscaled <- function(x, scale) x * scale
condf <- tibble(cname=inputs2$constraint_names,
                con=inputs2$constraints) %>%
  mutate(conx0=eval_f_obj_elements(inputs2$x0, inputs2),
         conx=eval_f_obj_elements(ip1a$solution, inputs2)) %>%
  left_join(scalef, by = "cname") %>%
  # mutate(scale=1) %>%
  mutate(across(c(con, conx0, conx), ~unscaled(.x, scale), .names = "{col}_raw")) %>%
  mutate(diff=conx_raw - con_raw,
         pdiff=diff / con_raw * 100)

condf %>%
  select(cname, contains("raw"), diff, pdiff) %>%
  arrange(desc(abs(pdiff)))

condf %>%
  filter(str_detect(cname, "pap")) %>%
  arrange(desc(abs(pdiff)))



count(inputs2$cc_sparse, cname, scale)



xlb <- rep(0, inputs2$n_variables)
xub <- rep(50, inputs2$n_variables)
opts_nlopt <- list("algorithm"="NLOPT_LD_LBFGS", # NLOPT_LD_MMA NLOPT_LD_LBFGS BUT NOT!: NLOPT_LD_SLSQP NLOPT_LD_MMA
             "check_derivatives" = FALSE,
             "xtol_rel"=1e-6, # default 1e-4 if 0, ignored
             "maxeval"=10)
a <- proc.time()
wfs1 <- nloptr(inputs2$x0,
                 eval_f=eval_f_wfs,
                 eval_grad_f = eval_grad_f_wfs,
                 lb = xlb, ub = xub,
                 opts = opts_nlopt, inputs=inputs2)
b <- proc.time()
b - a

str(wfs1)
eval_f_wfs(inputs2$x0, inputs2)
wfs1$objective
quantile(wfs1$solution)

cw <- check_weights(tmp$par)
cw %>% ht
sse_weights(cw)
cw %>%
  ggplot(aes(weight_total, weight)) +
  geom_point()

# NLOPT_LD_MMA 100 0.002563917  4.7 secs
# NLOPT_LD_LBFGS 100 2.410433e-12 2.2 secs



# augmented lagrangian 
inputs2$objscale <- 1

xlb <- rep(0, inputs2$n_variables)
xub <- rep(50, inputs2$n_variables)
inputs2$heqjac <- ajac
heqjac <- function(x, inputs) ajac

ajac <- nl.jacobian(inputs2$x0, heq, inputs=inputs2)
ajac2 <- nl.jacobian(runif(length(inputs2$x0)), heq, inputs=inputs2)
dim(ajac)
ajac[1:5, 1:10]
ajac2[1:5, 1:10]
inputs2$cc_addup_sparse$iweight_state[1:10]
nstates <- length(unique(inputs2$cc_addup_sparse$stabbr))
npids <- length(unique(inputs2$cc_addup_sparse$pid))
jac_dense <- matrix(0, nrow=npids, ncol = nstates * npids) # cannot allocate this
names(inputs2)


local_opts <- list("algorithm" = "NLOPT_LD_LBFGS", # NLOPT_LD_MMA NLOPT_LD_LBFGS
                   "xtol_rel"  = 1.0e-6 ) 
opts_nlopt <- list("algorithm" = "NLOPT_LD_AUGLAG",
                   "xtol_rel"  = 1.0e-6,
                   "maxeval"   = 10,
                   #"check_derivatives" = TRUE,
                   #"check_derivatives_print"='errors',
                   "local_opts" = local_opts)
a <- proc.time()
opt_al <- nloptr(inputs2$x0, # inputs2$x0, wfs1$solution
               eval_f=eval_f_wfs,
               eval_grad_f = eval_grad_f_wfs,
               eval_g_eq = heq,
               eval_jac_g_eq = heqjac,
               lb = xlb, ub = xub,
               opts = opts_nlopt, inputs=inputs2)
b <- proc.time()
b - a
str(opt_al)

library(Rsolnp)
a <- solnp(pars=inputs2$x0, fun=eval_f_wfs, eqfun = NULL, 
           eqB = NULL, LB = inputs2$xlb, UB = inputs2$xub, 
           control = list(outer.iter=2,
                          inner.iter=2), inputs=inputs2)



nloptr::auglag(x0=inputs2$x0, 
               fn = eval_f_wfs, 
               gr = eval_grad_f_wfs, 
               lower = xlb, upper = xub, 
               hin = hin,
               localsolver = c("COBYLA"),
               localtol = 1e-06, ineq2local = FALSE,
               nl.info = FALSE, control = list(), inputs=inputs2)

# NLOPT_LD_LBFGS 100 11.94 secs 7312.283
# NLOPT_LD_MMA 100 12.03 secs 12240.71, 200 23 secs 19523

eval_f_wfs(inputs2$x0, inputs2)
eval_f_wfs(wfs2$solution, inputs2)

library(alabama)
heq <- function(x, inputs) {
  tot_weights <- inputs$cc_addup_sparse %>%
    mutate(weight_state=iweight_state * x[j]) %>%
    group_by(ipid) %>%
    summarise(weight_state=sum(weight_state), .groups = "drop")
  tot_weights$weight_state - inputs2$weight_totals$weight_total
}

hin <- function(x, inputs) {
  tot_weights <- inputs$cc_addup_sparse %>%
    mutate(weight_state=iweight_state * x[j]) %>%
    group_by(ipid) %>%
    summarise(weight_state=sum(weight_state), .groups = "drop")
  
  cmin <- tot_weights$weight_state * 1.01 - inputs2$weight_totals$weight_total # all must be > 0
  cmax <- inputs2$weight_totals$weight_total - tot_weights$weight_state * .99 # all must be > 0
  c(cmin, cmax)
}

# hjfn <- function(x, inputs) ajac

# constrOptim.nl
# In optim(par = theta, fn = fun, gr = grad, control = control.optim,  :
#            method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
a <- proc.time()
tmp <- alabama::constrOptim.nl(par=inputs2$x0, 
              fn=eval_f_wfs, 
              gr = eval_grad_f_wfs,
              heq = heq,
              # heq.jac = hjfn,
              control.outer = list(itmax = 10, 
                                   trace=TRUE, 
                                   kkt2.check=FALSE,
                                   mu0=.9, 
                                   sig0=.9, # .9 worked well
                                   eps = 1e-8,
                                   method="L-BFGS-B"), 
              control.optim = list(trace=1, maxit=10), # , factr=1e-10, pgtol=1e-10
              inputs=inputs2)
b <- proc.time()
b - a
str(tmp)
eval_f_wfs(tmp$par, inputs2)
cw <- check_weights(tmp$par)
cw %>% ht
sse_weights(cw)
cw %>%
  ggplot(aes(weight_total, weight)) +
  geom_point()



spg(par=p0, fn=fn, gr=gr, project="projectLinear", 
    projectArgs=list(A=Amat, b=b, meq=meq))




# tmp <- constrOptim(theta=inputs2$x0,  f=eval_f_wfs, grad= eval_grad_f_wfs,
#                    ui, ci, mu = 1e-04, control = list(),
#             method = if(is.null(grad)) "Nelder-Mead" else "BFGS",
#             outer.iterations = 100, outer.eps = 1e-05, …,
#             hessian = FALSE)



auglag(x0, fn, gr = NULL, lower = NULL, upper = NULL, hin = NULL,
       hinjac = NULL, heq = NULL, heqjac = NULL,
       localsolver = c("COBYLA"), localtol = 1e-06, ineq2local = FALSE,
       nl.info = FALSE, control = list(), ...)

auglag(par, fn, gr, hin, hin.jac, heq, heq.jac, 
       control.outer=list(), control.optim = list(), ...)


eval_g_eq(inputs2$x0, inputs2)

eval_jac_g_addup(inputs2$x0, inputs2)

a <- proc.time()
wfs2 <- alabama::auglag(inputs2$x0,
               fn=eval_f_wfs,
               gr = eval_grad_f_wfs,
               heq = heq,
               control.outer = list(itmax = 10, 
                                    trace=TRUE, 
                                    kkt2.check=FALSE,
                                    # mu0=.9, 
                                    # sig0=.9, # .9 worked well
                                    # eps = 1e-8,
                                    method="L-BFGS-B"), 
               control.optim = list(trace=1, maxit=10), # , factr=1e-10, pgtol=1e-10
               inputs=inputs2)
b <- proc.time()
b - a

# ipopt here ----
# eval_g_addup(x, inputs)
# a <- proc.time()
# der <- check.derivatives(wadjust$xadj, func=eval_f_wfs, func_grad=eval_grad_f_wfs, inputs=inputs2)
# b <- proc.time()
# b - a



inputs2$constraints
eval_g(wfs1$solution, inputs2)
eval_g(wfs2$solution, inputs2)
eval_g(wfs3$solution, inputs2)

cor(wfs1$solution, wfs2$solution)
cor(wfs1$solution, wfs3$solution)
tibble(x=wfs1$solution, y=wfs3$solution) %>%
  ggplot(aes(x, y)) +
  geom_point()

inputs2$cc_addup_sparse %>%
  mutate(x=wfs3$solution,
         wnew=iweight_state * x) %>%
  group_by(ipid, pid) %>%
  summarise(weight_total=first(weight_total),
            iweight_state=sum(iweight_state),
            wnew=sum(wnew))



eval_g_eq <- function(x, inputs) {
  eval_g_addup(x, inputs)
}


# eval_g_eq <- function(x, inputs) {
#   constr <- eval_g_addup(x, inputs)
#   grad <- inputs$jacg_dense
#   return( list( "constraints"=constr, "jacobian"=grad ) )
# }

eval_f_wfs <- function(x, inputs){
  # constraints_vec <- calc_constraints(w, inputs$ccoef, inputs$target_names)
  constraints_vec <- eval_f_obj_elements(x, inputs)
  con_diff <- constraints_vec - inputs$constraints
  obj <- sum(con_diff^2) / inputs$objscale # * inputs$priority_weights)
  return(obj)
}

# inputs <- inputs2
# x <- inputs$x0

eval_grad_f_wfs <- function(x, inputs){
  # http://www.derivative-calculator.net/
  # https://www.symbolab.com/solver/partial-derivative-calculator
  # Example where we only have one target element in the objective function, and only have
  #   3 records and therefore 3 weights, and need to return a vector of 3 partial derivatives
  # Notation: w1, w2, and w3 are the first 3 weights in the vector of weights
  #           a1, a2, and a3 are the 3 constants they are multiplied by (e.g., wages1, wages2, wages3)
  #           t is the sum or target we are comparing to
  #           p is the priority weight of this element of the objective function
  #           calc is the calculated value of the target with the new weights
  # The objective function to minimize is p*(w1*a1 + w2*a2 + w3*a3 - s)^2
  # The first partial derivatives of the objective function are:
  #   wrt w1:          2*a1*p*(w1*a1 + w2*a2 + w3*a3 - t) = 2*a1*p*(calc - t)
  #   wrt w2:          2*a2*p*(w1*a1 + w2*a2 + w3*a3 - t) = 2*a2*p*(calc - t)
  #   wrt w3:          2*a3*p*(w1*a1 + w2*a2 + w3*a3 - t) = 2*a3*p*(calc - t)
  # If we have multiple target elements, each vector element will have a more complex sum
  
  constraints_vec <- eval_f_obj_elements(x, inputs)
  con_diff <- constraints_vec - inputs$constraints
  
  grad <- inputs$cc_sparse %>%
    mutate(iweight=inputs$iweight[j],
           x=x[j],
           diff=con_diff[i],
           grad= 2 * nzcc * x * diff) %>%
    group_by(j) %>%
    summarise(grad=sum(grad), .groups="drop")
  
  return(grad$grad / (inputs$objscale))
}


# inputs <- inputs2
# x <- inputs$x0

eval_g_addup <- function(x, inputs){
  tot_weights <- inputs$cc_addup_sparse %>%
    mutate(weight_state=iweight_state * x[j]) %>%
    group_by(ipid) %>%
    summarise(weight_state=sum(weight_state), .groups = "drop")
  tot_weights$weight_state
}

eval_f_obj_elements <- function(x, inputs) {
  # constraints that must hold in the solution - just give the LHS of the expression
  # return a vector where each element evaluates a constraint (i.e., sum of (x * a ccmat column), for each column)
  
  # ipoptr requires that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  
  # inputs$cc_sparse has the fixed constraint coefficients in sparse form in a dataframe that has:
  #   i -- the constraint number (based on those constraints we send to ipoptr, which could be a subset of all)
  #   j -- index into x (i.e., the variable number)
  #   nzcc
  
  constraint_tbl <- inputs$cc_sparse %>%
    group_by(i) %>%
    summarise(constraint_value=sum(nzcc * x[j]),
              .groups="keep")
  # the column constraint_value is a vector, in the order we want
  
  return(constraint_tbl$constraint_value)
}




eval_jac_g_addup <- function(x, inputs){
  # the Jacobian is the matrix of first partial derivatives of constraints (these derivatives may be constants)
  # this function evaluates the Jacobian at point x
  
  # return: a vector where each element gives a NONZERO partial derivative of constraints wrt change in x
  # so that the first m items are the derivs with respect to each element of first column of ccmat
  # and next m items are derivs with respect to 2nd column of ccmat, and so on
  # so that it returns a vector with length=nrows x ncolumns in ccmat
  
  # because constraints in this problem are linear, the derivatives are all constants
  
  # ipoptr requires that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  
  return(inputs$cc_addup_sparse$iweight_state)
}


get_targets2 <- function(.targets_wide, .target_vars, .pid, .weight_total){
  # .targets_wide dataframe that includes STATE, incgroup, nrecs, and the target variables
  # .target_vars  character vector of target variable names
  # cname is a unique constraint name
  
  # first, make a long data frame of the non-adding-up targets
  targets_long <- .targets_wide %>%
    pivot_longer(cols=all_of(.target_vars), names_to = "vname_calctype") %>%
    separate(vname_calctype, c("vname", "calctype"), sep="_", remove=FALSE) %>%
    mutate(cname=paste0(vname_calctype, "_", stabbr),
           targtype="aggregate") %>%
    select(stabbr, vname, calctype, targtype, vname_calctype, cname, value) %>%
    arrange(cname)
  
  # now the adding-up targets
  # we need person id and total weight
  # targets_adding_up <- tibble(pid=.pid, value=.weight_total) %>%
  #   mutate(cname=make_pid_cname(pid),
  #          vname="pid",
  #          targtype="addup") %>%
  #   select(vname, targtype, pid, cname, value) %>%
  #   arrange(cname)
  
  # we also need these as a vector -- make sure they are in the same order
  targets_df <- # bind_rows(targets_long, targets_adding_up) %>%
    targets_long %>%
    mutate(i=row_number()) %>%
    # select(i, STATE, vname, calctype, targtype, vname_calctype, pid, cname, value)
    select(i, stabbr, vname, calctype, targtype, vname_calctype, cname, value)
  
  targets_df
}



get_constraint_coefficients2 <- function(.incgroup_data, .target_vars, .iweights, .targets_df) {
  #.. dense matrix of constraint coefficients, not including the adding-up constraints
  # it has one record per person per state (if there are 50 states, it will have 50 records per person)
  cc_dense1 <- .iweights %>%
    select(pid, stabbr, iweight_state) %>%
    left_join(.incgroup_data %>% select(pid, all_of(.target_vars)), by = "pid") %>%
    mutate_at(vars(all_of(.target_vars)), list(~ . * iweight_state)) %>% # constraint coefficients
    arrange(pid, stabbr) %>%
    mutate(j=row_number()) %>% # j is an index for x, the variables we will solve for
    select(j, pid, stabbr, iweight_state, all_of(.target_vars))
  
  #.. sparse matrix of constraint coefficients, not including adding-up constraints
  cc_sparse1 <- cc_dense1 %>%
    select(j, pid, stabbr, all_of(.target_vars)) %>%
    pivot_longer(cols = all_of(.target_vars),
                 names_to="vname_calctype",
                 values_to = "nzcc") %>% # nzcc for nonzero constraint coefficient
    filter(nzcc!=0) %>%
    mutate(cname=paste0(vname_calctype, "_", stabbr)) %>%
    select(j, nzcc, pid, stabbr, cname)
  
  # Build the second part, for the adding up constraints. It will have 1 row per person per state,
  #   the same as the number of rows in cc_dense1 above
  # Each row will have the following variables:
  #   pid  the identifier of the person in the original data
  #   i  the row number for the constraint in the imaginary dense matrix, which is its position in the constraints vector
  #   j  the column number for the x variable in the imaginary dense matrix, which is its row number in the dense1 matrix
  #   cname the constraint name
  #   nzcc  the nonzero constraint coefficient, which is the amount by which the constraint value will change for a
  #      unit change in the x value, which will simply be the initial weight value
  # we can keep some other variables if we want
  
  # cc_sparse2 <- cc_dense1 %>%
  #   select(j, pid, STATE, nzcc=iweight_state) %>%
  #   mutate(cname=make_pid_cname(pid))
  
  cc_sparse <- # bind_rows(cc_sparse1, cc_sparse2) %>%
    cc_sparse1 %>%
    left_join(.targets_df %>% select(i, cname, targtype), by = c("cname")) %>% # get i and targtype from .targets_df
    arrange(i, j) # this ordering is crucial for the Jacobian
  
  return(cc_sparse)
}



get_inputs2 <- function(.targets_df, .iweights, .cc_sparse, .objscale=1, .p=2, .targtol=.01, .adduptol=0, .xub=20, .conscaling=FALSE, scale_goal=1){
  inputs_unscaled <- list()
  inputs_unscaled$p <- .p
  inputs_unscaled$iweight <- .iweights$iweight_state # the initial weight
  inputs_unscaled$cc_sparse <- .cc_sparse
  inputs_unscaled$constraints <- .targets_df$value
  inputs_unscaled$constraint_names <- .targets_df$cname
  inputs_unscaled$n_variables <- length(inputs_unscaled$iweight)
  inputs_unscaled$n_constraints <- length(inputs_unscaled$constraints)
  inputs_unscaled$n_targets <- nrow(inputs_unscaled$cc_sparse %>%
                                      filter(targtype=="aggregate") %>%
                                      select(cname) %>%
                                      distinct)
  inputs_unscaled$objscale <- .objscale
  
  conbounds <- get_conbounds(.constraints=inputs_unscaled$constraints, .n_targets=inputs_unscaled$n_targets, .targtol=.targtol, .adduptol=.adduptol)
  inputs_unscaled$clb <- conbounds$clb
  inputs_unscaled$cub <- conbounds$cub
  
  if(.conscaling==TRUE) inputs <- scale_inputs(inputs_unscaled, scale_goal) else inputs <- inputs_unscaled
  
  # finally, add xlb, xub, x0, and the relevant structures
  inputs$xlb <- rep(0, inputs$n_variables)
  inputs$xub <- rep(.xub, inputs$n_variables)
  inputs$x0 <- rep(1, inputs$n_variables)
  
  inputs$eval_jac_g_structure <- define_jac_g_structure_sparse(inputs$cc_sparse, ivar="i", jvar="j")
  inputs$eval_h_structure <- lapply(1:inputs$n_variables, function(x) x) # diagonal elements of our Hessian
  
  inputs
}


check_weights <- function(x){
  wcheck <- inputs2$weight_totals %>%
    mutate(weight=eval_g_addup(x, inputs2))
  wcheck
}

sse_weights <- function(wcheck){
  wcheck %>%
    summarise(sse2=sum((weight - weight_total)^2))
}

