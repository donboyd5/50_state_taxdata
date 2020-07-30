calc_targets_mat <- function(stub, targetvec, wtsdf=wts_all, pufdf=pufstrip){
  # wtsdf is a df that has state abbreviations for columns and tax returns
  #   for rows, for one or more subs; each cell is a weight
  # pufdf is tax return data
  whs <- get_whs_stub(wts_all, stub)
  xmat <-  get_xmat_stub(pufdf, stub, targetvec)
  t(whs) %*% xmat
}


get_whs_stub <- function(wts, stub){
  wts %>%
    filter(AGI_STUB == stub) %>%
    select(-c(AGI_STUB, sortid, wh, wh_sum)) %>% # keep just the state weights
    as.matrix
}

get_xmat_stub <- function(pufdf, stub, targetvec){
  pufdf %>%
    filter(AGI_STUB == stub) %>%
    select(all_of(targetvec)) %>%
    as.matrix
}


calc_targets <- function(wtsdf=wts_all, pufdf=pufstrip){
  stfull <- c(state.abb, "DC", "OA")
  stnames <- intersect(stfull, names(wtsdf)) %>% sort

  targetvec <- setdiff(names(pufdf), c("sortid", "AGI_STUB", "s006"))

  whslong <- wtsdf %>%
    pivot_longer(cols = all_of(stnames), names_to = "STATE", values_to = "whs")

  targs_wide <- whslong %>%
    left_join(pufdf, by=c("sortid", "AGI_STUB")) %>%
    group_by(AGI_STUB, STATE) %>%
    summarise(across(.cols = all_of(targetvec),
                     function(x) sum(x * whs, na.rm=TRUE)),
                     .groups="drop")

  targs_wide
}
