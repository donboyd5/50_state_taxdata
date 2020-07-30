# .iweights <- iweights
# .incgroup_data <- igdata
# .target_vars <- def_targs
# .targets_df <- targdf

get_cc_national <- function(.incgroup_data, .target_vars, .iweights, .targets_df) {
  #.. dense matrix of constraint coefficients for the nation
  # it has one record per person
  cc_dense1 <- .iweights %>%
    select(sortid, iweight) %>%
    left_join(.incgroup_data %>% select(sortid, all_of(.target_vars)), by = "sortid") %>%
    mutate(across(all_of(.target_vars),  function(x) x * iweight)) %>% # constraint coefficients
    arrange(sortid)

  # caution: row_number doesn't seem to work in parallel so I do this in 2 steps
  cc_dense1 <- cc_dense1 %>%
    mutate(j=1:nrow(cc_dense1)) %>% # j is an index for x, the variables we will solve for
    select(j, sortid, iweight, all_of(.target_vars))

  #.. sparse matrix of constraint coefficients, not including adding-up constraints
  cc_sparse1 <- cc_dense1 %>%
    pivot_longer(cols = all_of(.target_vars),
                 names_to="cname",
                 values_to = "nzcc") %>% # nzcc for nonzero constraint coefficient
    filter(nzcc!=0) %>%
    select(j, nzcc, sortid, cname)

  cc_sparse <- cc_sparse1 %>%
    left_join(.targets_df %>% select(i, cname), by = c("cname")) %>% # get i from .targets_df
    select(i, j, cname, nzcc, sortid) %>%
    arrange(i, j) # this ordering is crucial for the Jacobian

  return(cc_sparse)
}


get_cc_states <- function(.incgroup_data, .target_vars, .iweights, .targets_df) {
  #.. dense matrix of constraint coefficients, not including the adding-up constraints
  # it has one record per person per state (if there are 50 states, it will have 50 records per person)
  cc_dense1 <- .iweights %>%
    select(j, sortid, STATE, iweight) %>%
    left_join(.incgroup_data %>% select(sortid, all_of(.target_vars)), by = "sortid") %>%
    mutate(across(all_of(.target_vars), function(x) x * iweight)) %>% # constraint coefficients
    ungroup %>% # to be safe
    arrange(j) %>%
    select(j, sortid, STATE, iweight, all_of(.target_vars))

  #.. sparse matrix of constraint coefficients, not including adding-up constraints
  cc_sparse1 <- cc_dense1 %>%
    pivot_longer(cols = all_of(.target_vars),
                 names_to="var_calc",
                 values_to = "nzcc") %>% # nzcc for nonzero constraint coefficient
    filter(nzcc!=0) %>%
    mutate(cname=paste0(STATE, "_", var_calc),
                        targtype="aggregate") %>%
    separate(var_calc, c("var", "calc")) %>%
    select(j, nzcc, sortid, cname, STATE, var, calc, targtype)

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
  cc_sparse2 <- .iweights %>%
    select(j, sortid, STATE, nzcc=iweight) %>%
    mutate(cname=paste0("p", str_pad(sortid, width=8, side="left", pad="0")))

  cc_sparse <- bind_rows(cc_sparse1, cc_sparse2) %>%
    left_join(.targets_df %>% select(i, cname), by = c("cname")) %>% # get i from .targets_df
    arrange(i, j) # this ordering is crucial for the Jacobian

  return(cc_sparse)
}


get_conbounds <- function(.constraints, .n_targets, .targtol=.01){
  # define constraint lower and upper bounds -- for now they will be the same
  tol <- rep(0, length(.constraints))
  tol[1:.n_targets] <- .targtol

  clb <- ifelse(.constraints==0,
                -Inf,
                .constraints - abs(.constraints) * tol)
  cub <- ifelse(.constraints==0,
                +Inf,
                .constraints + abs(.constraints) * tol)
  list(clb=clb, cub=cub)
}

# .incgroup_data <- igdata
# .target_vars <- target_vars
# .iweights <- iweights3
# .targets_df <- targdf
# .cc_sparse <- cc_sparse

get_inputs <- function(.targets_df, .iweights, .cc_sparse, .objscale=1, .p=2, .targtol=.01, .xub=20, .conscaling=FALSE, scale_goal=1){
  inputs_unscaled <- list()
  inputs_unscaled$p <- .p
  inputs_unscaled$iweight <- .iweights$iweight # the initial weight
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

  conbounds <- get_conbounds(.constraints=inputs_unscaled$constraints, .n_targets=inputs_unscaled$n_targets, .targtol=.targtol)
  inputs_unscaled$clb <- conbounds$clb
  inputs_unscaled$cub <- conbounds$cub

  # if(.conscaling==TRUE) inputs <- scale_inputs(inputs_unscaled, scale_goal) else inputs <- inputs_unscaled
  inputs <- inputs_unscaled

  # finally, add xlb, xub, x0, and the relevant structures
  inputs$xlb <- rep(0, inputs$n_variables)
  inputs$xub <- rep(.xub, inputs$n_variables)
  inputs$x0 <- rep(1, inputs$n_variables)

  inputs$eval_jac_g_structure <- define_jac_g_structure_sparse(inputs$cc_sparse, ivar="i", jvar="j")
  inputs$eval_h_structure <- lapply(1:inputs$n_variables, function(x) x) # diagonal elements of our Hessian

  inputs
}

get_inputs_national <- function(.targets_df, .iweights, .cc_sparse, .objscale=1, .p=2, .targtol=.01, .xub=20, .conscaling=FALSE, scale_goal=1){
  inputs_unscaled <- list()
  inputs_unscaled$p <- .p
  inputs_unscaled$iweight <- .iweights$iweight # the initial weight
  inputs_unscaled$cc_sparse <- .cc_sparse
  inputs_unscaled$constraints <- .targets_df$value
  inputs_unscaled$constraint_names <- .targets_df$cname
  inputs_unscaled$n_variables <- length(inputs_unscaled$iweight)
  inputs_unscaled$n_constraints <- length(inputs_unscaled$constraints)

  inputs_unscaled$objscale <- .objscale

  conbounds <- get_conbounds(.constraints=inputs_unscaled$constraints, .n_targets=inputs_unscaled$n_constraints, .targtol=.targtol)
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

scale_inputs <- function(inputs, scale_goal=1){

  ccsum <- inputs$cc_sparse %>%
    group_by(i, cname) %>%
    summarise(nzmin=min(nzcc), nzmdn=median(nzcc), nzmax=max(nzcc),
              .groups="drop")

  scale_factors <- tibble(cname=inputs$constraint_names, cvalue=inputs$constraints) %>%
    left_join(ccsum, by = "cname") %>%
    mutate(scale=abs(nzmax / scale_goal)) %>%
    mutate_at(vars(cvalue, starts_with("nz")), list(~ . / scale))

  inputs_scaled <- inputs

  inputs_scaled$cc_sparse <- inputs$cc_sparse %>%
    left_join(scale_factors %>% select(cname, scale), by = "cname") %>%
    mutate(nzcc = nzcc / scale)

  inputs_scaled$constraints <- scale_factors$cvalue

  inputs_scaled$clb <- inputs$clb / scale_factors$scale
  inputs_scaled$cub <- inputs$cub / scale_factors$scale

  return(inputs_scaled)
}


one_stub <- function(stub){
  # reweight nationally for a single stub
  print(paste0("Reweighting for AGI_STUB: ", stub)) # won't print if parallel

  a <- proc.time()

  # iweights is vector, length n
  # xmat is n x k matrix of data
  # targets is vector, length k
  # target names is vector, length k
  # tol is vector, length k

  iweights <- wts_all %>% filter(AGI_STUB==stub) %>% .$wh_sum

  targs_stub <- ustargs %>%
    filter(AGI_STUB == stub) %>%
    pivot_longer(-AGI_STUB, names_to = "cname", values_to = "value") %>%
    filter(!is.na(value))

  targets <- targs_stub$value
  tol <- abs(targets) * .001
  target_names <- targs_stub$cname

  data_stub <- pufstrip %>%
    filter(AGI_STUB == stub) %>%
    select(AGI_STUB, sortid, all_of(target_names))

  xmat <- data_stub %>%
    select(all_of(target_names)) %>%
    as.matrix

  fname <- here::here("ignore", paste0("stubrw_", stub, ".out"))
  optlist <- list(file_print_level = 5,
                  output_file = fname,
                  linear_solver = "ma57")

  result <- reweight(iweights, targets, target_names, tol, xmat, optlist=optlist)

  b <- proc.time()

  result
}
