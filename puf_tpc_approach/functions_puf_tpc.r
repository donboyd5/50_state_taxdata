
remove_suffix <- function(suffix, target_vars) {
  # return a named vector with the base names of target variables
  var_backend <- str_extract(target_vars, "_.*$") # gets the ending part of each variable name
  base_vars <- target_vars[which(var_backend==suffix)] %>% str_remove(suffix)
  names(base_vars) <- base_vars
  base_vars
}


check_xmat <- function(xmat){
  xmat_rank <- qr(xmat[,])$rank

  rankifremoved <- sapply(1:ncol(xmat), function (x) qr(xmat[,-x])$rank)

  lindep_cols <- which(rankifremoved == max(rankifremoved)) # indices of columns that are linearly dependent
  if(length(lindep_cols) == xmat_rank) lindep_cols <- 0

  check <- list()
  check$rank <- xmat_rank
  check$rankifremoved <- rankifremoved
  check$lindep_cols <- lindep_cols
  check$remove <- max(lindep_cols)
  check
}


one_state <- function(st, incgroup, target_vars_df, quiet=TRUE, method='LM', maxiter=50){
  print(sprintf("State: %s, Income group: %i", st, incgroup))

  # grab the targets for this STATE-AGI_STUB combination
  target_vars <- target_vars_df %>%
    filter(STATE==st, AGI_STUB==incgroup) %>%
    .$target_vec %>%
    unlist

  xmat <- pufstrip %>%
    filter(AGI_STUB==incgroup) %>%
    select(all_of(target_vars)) %>%
    as.matrix

  # only check once!!??
  check <- check_xmat(xmat)
  # print(check)
  if(check$remove > 0) {
    print(paste0("WARNING: Removing linearly dependent column: ", target_vars[check$remove]))
    print("CAUTION: Not checking for further linear dependence...")
    xmat <- xmat[, -check$remove]
    target_vars <- target_vars[-check$remove]
  }

  targs <- puf_targ %>%
    filter(AGI_STUB==incgroup) %>%
    select(AGI_STUB, STATE, all_of(target_vars)) %>%
    mutate(STATE=ifelse(STATE==st, STATE, "ZZ")) %>%
    group_by(AGI_STUB, STATE) %>%
    summarise(across(all_of(target_vars), sum), .groups="drop")

  targmat <- targs %>%
    select(all_of(target_vars)) %>%
    as.matrix
  rownames(targmat) <- targs$STATE
  targmat

  sortid <- pufstrip %>%
    filter(AGI_STUB==incgroup) %>%
    .$sortid

  wh <- pufstrip %>%
    filter(AGI_STUB==incgroup) %>%
    .$s006

  result <- geoweight(wh=wh, xmat=xmat, targets=targmat, method = method, quiet=quiet, maxiter=maxiter)
  result$sortid <- sortid

  return(result)
}

