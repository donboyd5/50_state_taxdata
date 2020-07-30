source(here::here("include", "functions_ipopt.r"))

targets_wide <- targets_wide %>%
  rename(stabbr=STATE)
  # rename(STATE=stabbr)

# djb start ----
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



cc_sparse2 <- get_constraint_coefficients2(.incgroup_data=incgroup_data,
                                           .target_vars=target_vars,
                                           .iweights=iweights,
                                           .targets_df=targets_df2) # can take a while on big problems
count(cc_sparse2, targtype)


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


inputs2 <- get_inputs(.targets_df=targets_df2, .iweights=iweights, .cc_sparse=cc_sparse2,
                      .targtol=.05, .adduptol=0, .xub=50, .conscaling=TRUE, scale_goal=1)
names(inputs2)
inputs2$n_variables
inputs2$n_constraints
inputs2$n_targets

eval_f_wfs(inputs2$x0, inputs2)


system.time(tmp <- check.derivatives(.x=inputs2$x0, func=eval_f_wfs, func_grad=eval_grad_f_wfs,
                         check_derivatives_print='none', inputs=inputs2))
str(tmp)
tb <- tibble(ad=tmp$analytic, fd=tmp$finite_difference, re=tmp$relative_error, flag=tmp$flag_derivative_warning, i=1:length(flag))
ht(tb)
tb %>% filter(flag)

xlb <- rep(0, inputs2$n_variables)
xub <- rep(50, inputs2$n_variables)
opts_nlopt <- list("algorithm"="NLOPT_LD_LBFGS", # NLOPT_LD_MMA NLOPT_LD_LBFGS BUT NOT!: NLOPT_LD_SLSQP NLOPT_LD_MMA
             "check_derivatives" = FALSE,
             "xtol_rel"=1e-4, # default 1e-4 if 0, ignored
             "maxeval"=50)
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



opts_ipopt <- list("print_level" = 0,
             "file_print_level" = 5, # integer
             "max_iter"= 100,
             "linear_solver" = "ma57", # mumps pardiso ma27 ma57 ma77 ma86 ma97
             # "mehrotra_algorithm" = "yes",
             "obj_scaling_factor" = 1e6, # 1e-3, # default 1; 1e-1 pretty fast to feasible but not to optimal
             "nlp_scaling_max_gradient" = 100, # default is 100 - seems good
             "output_file" = here::here("out", "wfs3os.out"))

wfs3 <- ipoptr(x0 = inputs2$x0,
               lb = inputs2$xlb,
               ub = inputs2$xub,
               eval_f = eval_f_wfs, # arguments: x, inputs; eval_f_xtop eval_f_xm1sq eval_f_wfs
               eval_grad_f = eval_grad_f_wfs, # eval_grad_f_xm1sq, # eval_grad_f_xtop eval_grad_f_xm1sq
               # eval_g = eval_g, # constraints LHS - a vector of values
               # eval_jac_g = eval_jac_g,
               # eval_jac_g_structure = inputs$eval_jac_g_structure,
               # eval_h = eval_h_xm1sq, # the hessian is essential for this problem eval_h_xtop eval_h_xm1sq
               # eval_h_structure = inputs$eval_h_structure,
               # constraint_lb = inputs$clb,
               # constraint_ub = inputs$cub,
               opts = opts_ipopt,
               inputs = inputs2)




inputs2$constraints
eval_g(wfs1$solution, inputs2)
eval_g(wfs3$solution, inputs2)

cor(wfs1$solution, wfs3$solution)

tibble(x=wfs1$solution, y=wfs3$solution) %>%
  ggplot(aes(x, y)) +
  geom_point()


eval_f_wfs <- function(x, inputs){
  # constraints_vec <- calc_constraints(w, inputs$ccoef, inputs$target_names)
  constraints_vec <- eval_g(x, inputs)
  con_diff <- constraints_vec - inputs$constraints
  obj <- sum(con_diff^2) # * inputs$priority_weights)
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
  
  constraints_vec <- eval_g(x, inputs)
  con_diff <- constraints_vec - inputs$constraints
  
  grad <- inputs$cc_sparse %>%
    mutate(iweight=inputs$iweight[j],
           x=x[j],
           diff=con_diff[i],
           grad= 2 * nzcc * x * diff) %>%
    group_by(j) %>%
    summarise(grad=sum(grad), .groups="drop")
  
  return(grad$grad)
}


