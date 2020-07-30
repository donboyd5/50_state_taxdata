
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
    mutate(nrecs=1 / {{.weight}}) %>%
    group_by(...) %>%
    summarise_at(vars(nrecs, all_of(.sum_vars)),
                 list(~ sum(. * {{.weight}})),
                 .groups="keep") %>%
    ungroup
}


get_tabs <- function(summary_vals){
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
    select(stabbr, incgroup, nrecs, pop_nnz, ends_with("sum")) %>%
    mutate_at(vars(ends_with("sum")), list(~ . / pop_nnz)) %>%
    kable(digits=0, format.args = list(big.mark=","), format="rst")
  
  # percentages of population that have nonzero values
  tablist$pctnz <- summary_vals %>%
    select(stabbr, incgroup, nrecs, pop_nnz, ends_with(("_nnz"))) %>%
    mutate_at(vars(ends_with("_nnz")), list(pct= ~ . / pop_nnz * 100)) %>%
    select(-ends_with("_nnz")) %>%
    kable(digits=1, format.args = list(big.mark=","), format="rst")
  
  return(tablist)
}


#**********************************************************************************************#
# Target preparation functions ----
#**********************************************************************************************#

get_targets <- function(.targets_wide, .target_vars, .pid, .weight_total){
  # cname is a unique constraint name
  
  # first, make a long data frame of the non-adding-up targets
  targets_long <- .targets_wide %>%
    pivot_longer(cols=all_of(.target_vars), names_to = "vname_targtype") %>%
    separate(vname_targtype, c("vname", "targtype"), sep="_", remove=FALSE) %>%
    mutate(cname=paste0(vname_targtype, "_", stabbr)) %>%
    select(stabbr, incgroup, vname, targtype, vname_targtype, cname, nrecs, value) %>%
    arrange(cname)
  
  # now the adding-up targets
  # we need person id and total weight
  targets_adding_up <- tibble(pid=.pid, value=.weight_total) %>%
    mutate(cname=paste0("p", str_pad(pid, width=8, pad="0"))) %>%
    select(cname, pid, value) %>%
    arrange(cname)
  
  # we also need these as a vector -- make sure they are in the same order
  constraint_names <- c(targets_long$cname, targets_adding_up$cname)
  constraints <- c(targets_long$value, targets_adding_up$value)
  
  targlist <- list()
  targlist$targets_long <- targets_long
  targlist$targets_adding_up <- targets_adding_up
  targlist$constraint_names <- constraint_names
  targlist$constraints <- constraints
  targlist
}


#**********************************************************************************************#
# Initial weights ----
#**********************************************************************************************#
get_initial_weights <- function(.targets_wide, .base_data){
  # define state_shares: from some source (e.g., SOI HT2) we will know the 
  # fraction of weighted records in an income group that come from each state
  state_shares <- .targets_wide %>%
    select(stabbr, incgroup, pop_nnz) %>%
    mutate(st_share=pop_nnz / sum(pop_nnz)) %>%
    select(-pop_nnz)
  
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
get_constraint_coefficients <- function(.base_data, .target_vars, .iweights, .targlist) {
  #.. dense matrix of constraint coefficients, not including the adding-up constraints
  cc_dense1 <- .iweights %>%
    select(pid, stabbr, iweight_state) %>%
    left_join(.base_data %>% select(pid, all_of(.target_vars)), by = "pid") %>%
    mutate_at(vars(all_of(.target_vars)), list(~ . * iweight_state)) %>% # constraint coefficients
    arrange(pid, stabbr) %>%
    mutate(j=row_number()) %>% # j is an index for x, the variables we will solve for
    select(j, pid, stabbr, iweight_state, all_of(.target_vars))
  
  #.. sparse matrix of constraint coefficients, not including adding-up constraints
  cc_sparse1 <- cc_dense1 %>%
    select(j, pid, stabbr, all_of(.target_vars)) %>%
    pivot_longer(cols = all_of(.target_vars), 
                 names_to="vname_targtype",
                 values_to = "nzcc") %>% # nzcc for nonzero constraint coefficient
    filter(nzcc!=0) %>%
    separate(vname_targtype, c("vname", "targtype"), sep="_", remove=FALSE) %>%
    mutate(cname=paste0(vname_targtype, "_", stabbr),
           i=match(cname, .targlist$constraint_names),
           ctype="target") %>%
    select(i, j, nzcc, pid, ctype, cname, vname, vname_targtype, stabbr)
  
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
    left_join(.targlist$targets_adding_up %>% select(cname, pid), by="pid") %>%
    mutate(i=match(cname, .targlist$constraint_names), # i is index for the constraint column
           ctype="addup")
  
  cc_sparse <- bind_rows(cc_sparse1, cc_sparse2) %>%
    arrange(i, j) # this ordering is crucial for the Jacobian
  
  return(cc_sparse)
}


#**********************************************************************************************#
# Build inputs list ----
#**********************************************************************************************#

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


get_inputs <- function(.targlist, .iweights, .cc_sparse, .objscale=1, .p=2, .targtol=.01, .xub=20, .conscaling=FALSE, scale_goal=1){
  inputs_unscaled <- list()
  inputs_unscaled$p <- .p
  inputs_unscaled$iweight <- .iweights$iweight_state # the initial weight
  inputs_unscaled$cc_sparse <- .cc_sparse
  inputs_unscaled$constraints <- .targlist$constraints
  inputs_unscaled$constraint_names <- .targlist$constraint_names
  inputs_unscaled$n_variables <- length(inputs_unscaled$iweight)
  inputs_unscaled$n_constraints <- length(inputs_unscaled$constraints)
  inputs_unscaled$n_targets <- nrow(inputs_unscaled$cc_sparse %>% 
                                      filter(ctype=="target") %>%
                                      select(cname) %>% 
                                      distinct)
  inputs_unscaled$objscale <- .objscale
  
  conbounds <- get_conbounds(.constraints=inputs_unscaled$constraints, .n_targets=inputs_unscaled$n_targets, .targtol=.01)
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

