
#**********************************************************************************************#
# Data preparation functions ----
#**********************************************************************************************#
make_target_names <- function(lnames, csuffixes){
  # lnames: list of variable-name vectors to add suffixes to
  # csuffixes: character vector of suffixes -- one suffix per variable-name vector
  
  # returns a character vector of variable names concatenated with their respective suffixes
  
  paste_suffix <- function(vnames, suffix){
    paste0(vnames, "_", suffix)
  }
  
  map2(lnames, csuffixes, paste_suffix) %>%
    unlist
}


prep_data <- function(data, target_vars){
  # prepare the FULL data set (ALL income groups)
  #   data: input data frame
  #   target_vars: character vector of variable names with suffixes:
  #     _nnz -- variables for which we want the weighted number of nonzero values
  #     _sum -- variable for which we want the weighted sum
  #     later we can expand this to have _npos (weighted # of positive values), _sumpos, _nneg, _sumneg, if
  #     we ever get to the point where we have targets for those kinds of variables
  #     
  #     For example, if target_vars includes "wagp_nnz" it means we want a variable that can be used to 
  #     get the number nonzero values of wagp (wages). If target_vars includes "ssip_sum", we want a variable
  #     that can be used to get the sum of values for ssip (SSI income).
  
  #     If target_vars contains "wagp_nnz" we will create a variable named wagp_nnz that is 1 for every record
  #     in which wagp is nonzero and 0 otherwise.
  
  #     If target_vars contains "ssip_sum" we will create a variable named ssip_sum that has the value of ssip
  #     for every record.
  
  #     The variables specified in target_values before the suffix, MUST exist in data. Eventually add
  #     error checking.
  
  # return: the data frame, enhanced with created variables, as specified in target_vars
  
  # create variables that when weighted and summed will produce values to compare to targets
  # currently includes sum and nnz, npos, sumpos, nneg, sumneg
  
  # trick mutate_at by naming the vectors of variables to be transformed, so that
  # it will give desired names to the constructed variables
  getbase <- function(suffix) {
    var_backend <- str_extract(target_vars, "_.*$") # gets the ending part of each variable name
    base_vars <- target_vars[which(var_backend==suffix)] %>% str_remove(suffix)
    names(base_vars) <- base_vars
    base_vars
  }
  
  nnz_vars <-getbase("_nnz") # variables for which we want the weighted number of nonzero values
  sum_vars <- getbase("_sum") # variables for which we want the weighted sum
  sumneg_vars <- getbase("_sumneg") # variables for which we want the weighted sum
  
  data2 <- data %>%
    mutate_at(nnz_vars,
              list(nnz = ~ 1 * (. != 0))) %>%
    mutate_at(sum_vars, 
              list(sum = ~ . * (. != 0))) %>%
    mutate_at(sumneg_vars, 
              list(sumneg = ~ . * (. < 0)))
  
  return(data2)
}


prep_incgroup <- function(.data, .incgroup, .weight, .target_vars){
  # prepare data for a single income group for a specific set of target variables
  # .data:         data frame that has been prepared for targeting
  # .incgroup:     integer that identifies the income group we will analyze
  # .weight:       bare name of the weight variable (no quotes around it)
  # .target_vars:  character vector of names of target variables (including suffixes)
  
  # check that data has variables named pid and incgroup
  check_vars <- setdiff(c("pid", "incgroup"), names(.data))
  if(length(check_vars) > 0){ 
    errmsg <- paste0("ERROR: .data must have variables named pid and incgroup")
    print(errmsg)
    return(NULL)
  }
  
  # check that .weight is in the data
  wname <- paste0(enexpr(.weight)) # there must be a better way...
  if(!wname %in% names(.data)) {
    errmsg <- paste0("ERROR: Weight variable ", wname, " not in .data.")
    print(errmsg)
    return(NULL)
  }
  
  # check that all desired targets are in the data
  check_targets <- setdiff(.target_vars, names(.data))
  if(length(check_targets) > 0 ) {
    errmsg <- paste0("ERROR: The following target variables are not in .data: ", paste(check_targets, collapse = ", "))
    print(errmsg)
    return(NULL)
  }
  
  # if we get here, function inputs are good
  base_data <- .data %>%
    filter(incgroup==.incgroup) %>%
    select(pid, incgroup, weight_total={{.weight}}, all_of(.target_vars))
  
  base_data
}



#**********************************************************************************************#
# Data summarization functions ----
#**********************************************************************************************#

count_recs <- function(data){
  # crosstab of number of records by stabbr (long) and incgroup (wide), with margin totals
  # as of May 16, 2020 this requires the dev version of dplyr; soon on CRAN
  #   remotes::install_github("tidyverse/dplyr")
  count(data, incgroup, stabbr) %>%
    pivot_wider(names_from = incgroup, values_from = n) %>%
    add_row(stabbr="all") %>%
    mutate_at(vars(-stabbr), list(~ ifelse(stabbr=="all", sum(., na.rm=TRUE), .))) %>%
    rowwise() %>%
    mutate(sum=sum(c_across(-stabbr)))
}


get_summary_vals <- function(.data, .weight, .sum_vars, ...){
  # get grouped weighted sums for a data frame
  .data %>%
    group_by(...) %>%
    summarise(nrecs=n(),
              across(.cols=(all_of(.sum_vars)), ~ sum(.x * {{.weight}})),
              .groups="keep") %>%
    ungroup
}


get_tabs <- function(summary_vals, .popvar){
  tablist <- list()
  
  # sums of weighted values:
  tablist$wsum <- summary_vals %>%
    select(-ends_with(("_nnz"))) %>%
    kable(digits=0, format.args = list(big.mark=","), format="rst")
  
  # sums of weighted counts:
  tablist$wcount <- summary_vals %>%
    select(stabbr, incgroup, ends_with(("_nnz"))) %>%
    kable(digits=0, format.args = list(big.mark=","), format="rst")
  
  # weighted means
  tablist$wmean <- summary_vals %>%
    select(stabbr, incgroup, nrecs, {{.popvar}}, ends_with("sum")) %>%
    mutate_at(vars(ends_with("sum")), list(~ . / {{.popvar}})) %>%
    kable(digits=0, format.args = list(big.mark=","), format="rst")
  
  # percentages of population that have nonzero values
  tablist$pctnz <- summary_vals %>%
    select(stabbr, incgroup, nrecs, {{.popvar}}, ends_with(("_nnz"))) %>%
    mutate_at(vars(ends_with("_nnz")), list(pct= ~ . / {{.popvar}} * 100)) %>%
    select(-ends_with("_nnz")) %>%
    kable(digits=1, format.args = list(big.mark=","), format="rst")
  
  return(tablist)
}


#**********************************************************************************************#
# Target preparation functions ----
#**********************************************************************************************#
make_pid_cname <- function(pid) {paste0("p", str_pad(pid, width=8, pad="0"))}


get_targets <- function(.targets_wide, .target_vars, .pid, .weight_total){
  # .targets_wide dataframe that includes stabbr, incgroup, nrecs, and the target variables
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
  targets_adding_up <- tibble(pid=.pid, value=.weight_total) %>%
    mutate(cname=make_pid_cname(pid),
           vname="pid",
           targtype="addup") %>%
    select(vname, targtype, pid, cname, value) %>%
    arrange(cname)
  
  # we also need these as a vector -- make sure they are in the same order
  targets_df <- bind_rows(targets_long, targets_adding_up) %>%
    mutate(i=row_number()) %>%
    select(i, stabbr, vname, calctype, targtype, vname_calctype, pid, cname, value)
  
  targets_df
}


#**********************************************************************************************#
# Initial weights ----
#**********************************************************************************************#
get_initial_weights <- function(.targets_wide, .base_data, .popvar){
  # define state_shares: from some source (e.g., SOI HT2) we will know the 
  # fraction of weighted records in an income group that come from each state
  state_shares <- .targets_wide %>%
    select(stabbr, incgroup, {{.popvar}}) %>%
    mutate(st_share={{.popvar}} / sum({{.popvar}})) %>%
    select(-{{.popvar}})
  
  iweights <- expand_grid(pid=unique(.base_data$pid), stabbr=.targets_wide$stabbr) %>%
    arrange(pid, stabbr) %>%
    left_join(.base_data %>% select(pid, incgroup, weight_total), by = "pid") %>%
    left_join(state_shares, by = c("stabbr", "incgroup")) %>%
    mutate(iweight_state=weight_total * st_share) %>%
    select(-st_share)
  iweights
}


#**********************************************************************************************#
# Constraint coefficients ----
#**********************************************************************************************#
get_constraint_coefficients <- function(.incgroup_data, .target_vars, .iweights, .targets_df) {
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
  
  cc_sparse2 <- cc_dense1 %>%
    select(j, pid, stabbr, nzcc=iweight_state) %>%
    mutate(cname=make_pid_cname(pid))
  
  cc_sparse <- bind_rows(cc_sparse1, cc_sparse2) %>%
    left_join(.targets_df %>% select(i, cname, targtype), by = c("cname")) %>% # get i and targtype from .targets_df
    arrange(i, j) # this ordering is crucial for the Jacobian
  
  return(cc_sparse)
}


#**********************************************************************************************#
# Build inputs list ----
#**********************************************************************************************#

get_conbounds_old <- function(.constraints, .n_targets, .targtol=.01){
  # define constraint lower and upper bounds -- for now they will be the same
  tol <- rep(0, length(.constraints))
  
  if(length(.constraints) > .n_targets){
    tol[(.n_targets + 1):length(.constraints)] <- .adduptol
  }
  
  clb <- ifelse(.constraints==0,
                -Inf,
                .constraints - abs(.constraints) * tol)
  cub <- ifelse(.constraints==0,
                +Inf, 
                .constraints + abs(.constraints) * tol)
  list(clb=clb, cub=cub)
}

get_conbounds <- function(.constraints, .n_targets, .targtol=.01, .adduptol=0){
  # define constraint lower and upper bounds -- for now they will be the same
  tol <- rep(0, length(.constraints))
  tol[1:.n_targets] <- .targtol
  
  if(length(.constraints) > .n_targets){
    tol[(.n_targets + 1):length(.constraints)] <- .adduptol
  }
  
  clb <- ifelse(.constraints==0,
                -Inf,
                .constraints - abs(.constraints) * tol)
  
  cub <- ifelse(.constraints==0,
                +Inf,
                .constraints + abs(.constraints) * tol)
  
  list(clb=clb, cub=cub)
}

get_inputs <- function(.targets_df, .iweights, .cc_sparse, .objscale=1, .p=2, .targtol=.01, .adduptol=0, .xub=20, .conscaling=FALSE, scale_goal=1){
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
  
  if(.conscaling==TRUE) inputs <- scale_inputs(inputs_unscaled, scale_goal) else{
    inputs <- inputs_unscaled
  } 
  
  # finally, add xlb, xub, x0, and the relevant structures
  inputs$xlb <- rep(0, inputs$n_variables)
  inputs$xub <- rep(.xub, inputs$n_variables)
  inputs$x0 <- rep(1, inputs$n_variables)
  
  inputs$eval_jac_g_structure <- define_jac_g_structure_sparse(inputs$cc_sparse, ivar="i", jvar="j")
  inputs$eval_h_structure <- lapply(1:inputs$n_variables, function(x) x) # diagonal elements of our Hessian
  
  inputs
}


get_inputs_old <- function(.targets_df, .iweights, .cc_sparse, .objscale=1, .p=2, .targtol=.01, .xub=20, .conscaling=FALSE, scale_goal=1){
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
  
  # conbounds <- get_conbounds(.constraints=inputs_unscaled$constraints, .n_targets=inputs_unscaled$n_targets, .targtol=.01)
  
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

