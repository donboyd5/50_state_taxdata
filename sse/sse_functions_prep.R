
get_targname <- function(vname_calctype, stabbr) paste0(vname_calctype, "_", stabbr)

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



scale_inputs_sse <- function(inputs_unscaled, scale_goal=1, scale_type="max"){
  
  coef_sum <- inputs_unscaled$ofe_sparse %>%
    group_by(targnum, targname) %>%
    summarise(nzmin=min(nz_coef), 
              nzmdn=median(nz_coef), 
              nzmax=max(nz_coef),
              .groups="drop")
  
  target_scale_factors <- inputs_unscaled$targets_df %>%
    left_join(coef_sum, by = c("targnum", "targname")) %>%
    mutate(target_unscaled = target,
           # scale_type = scale_type_in,
           scale_numerator = case_when(
             scale_type != "max" ~ nzmdn,
             scale_type == "max" & target >= 0 ~ nzmax,
             scale_type == "max" & target < 0 ~ nzmin,
             TRUE ~ NA_real_),
           scale=abs(scale_numerator) / scale_goal) %>%
    mutate_at(vars(target, starts_with("nz")), list(~ . / scale))
  
  # create the scaled inputs
  inputs_scaled <- inputs_unscaled
  inputs_scaled$ofe_sparse <- inputs_unscaled$ofe_sparse %>%
    left_join(target_scale_factors %>% 
                select(targname, scale), 
              by = "targname") %>%
    mutate(nz_coef = nz_coef / scale)
  
  inputs_scaled$targets_df <- target_scale_factors
  inputs_scaled$targets <- inputs_scaled$targets_df$target
  return(inputs_scaled)
}


get_inputs_sse <- function(targets_df, iweights, ofe_sparse, target_scaling=FALSE, scale_goal=1, scale_type="max", xlb=0, xub=50){
  inputs_unscaled <- list()
  inputs_unscaled$iweights_df <- iweights
  inputs_unscaled$iweight <- iweights$iweight_state # the initial weight
  inputs_unscaled$ofe_sparse <- ofe_sparse
  inputs_unscaled$targets_df <- targets_df
  inputs_unscaled$targets <- targets_df$target
  inputs_unscaled$n_variables <- length(inputs_unscaled$iweight)
  inputs_unscaled$n_targets_type <- count(targets_df, vname, calctype)
  inputs_unscaled$n_targets_all <- length(inputs_unscaled$targets)
  
  inputs_unscaled$constraints <- iweights %>%
    group_by(i, pid) %>%
    summarise(weight_total=first(weight_total), .groups="drop") %>%
    .$weight_total
  
  inputs_unscaled$cc_sparse <- iweights %>%
    select(i, pid, j, xname, stabbr, iweight_state)
  
  # finally, add xlb, xub, x0, and the relevant structures
  inputs_unscaled$xlb <- rep(xlb, inputs_unscaled$n_variables)
  inputs_unscaled$xub <- rep(xub, inputs_unscaled$n_variables)
  inputs_unscaled$x0 <- rep(1, inputs_unscaled$n_variables)
  
  inputs_unscaled$eval_jac_g_structure <- define_jac_g_structure_sparse(inputs_unscaled$cc_sparse, ivar="i", jvar="j")
  inputs_unscaled$objscale <- 1
  
  if(target_scaling==TRUE) inputs <- scale_inputs_sse(inputs_unscaled, scale_goal, scale_type) else{
    inputs_unscaled$targets_df$scale <- 1
    inputs <- inputs_unscaled
  } 
  
  inputs
}

