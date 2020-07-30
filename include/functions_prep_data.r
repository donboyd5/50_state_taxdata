
prep_ACS_sample <- function(fn){
  samp1 <- readRDS(here::here("data", fn)) %>% 
    select(-nrecs, -pop) # note that we no longer need nrecs; pop ordinarily would not be in the data so drop here and create later
  
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
  
  # prepare data by creating variables with those names:
  #   nnz, nneg, and npos will have 1 for rows where the respective variable is nz, neg, or pos, respectively, and 0 otherwise
  #   sum will have its respective variable's value
  #   sumneg and sumpos will have the variable's value if negative or positive, respectively, and 0 otherwise
  samp <- prep_data(samp2, possible_target_vars)
  
  # Create a data frame with all targets for all states and income groups ----
  # for the PUF, we will create this using information from Historical Table 2
  # for the ACS, we construct the targets from the ACS data itself
  all_target_values <- get_summary_vals(samp, .weight=pwgtp, .sum_vars=possible_target_vars, stabbr, incgroup)

  acslist <- list()
  acslist$data <- samp
  acslist$possible_target_vars <- possible_target_vars
  acslist$all_target_values <- all_target_values
  acslist
}


make_acs_problem <- function(acslist, target_incgroup=2, target_vars){
  samp <- acslist$data %>%
    filter(incgroup==target_incgroup)
  
  possible_target_vars <- acslist$possible_target_vars
  all_target_values <- acslist$all_target_values
  
  # wrap everything we need for a single income group into a function that returns a list ----
  
  # define target values and states, for this income group
  targets_wide <- all_target_values %>%
    filter(incgroup==target_incgroup) %>%
    select(stabbr, incgroup, nrecs, all_of(target_vars)) # a small list of variables to target; we have nrecs because we created it in summary_vals
  
  wh <- samp %>% .$pwgtp
  targets <- targets_wide[, target_vars] %>% as.matrix
  xmat <- samp %>% filter(incgroup==target_incgroup) %>% .[, target_vars] %>% as.matrix
  
  pacs <- list()
  pacs$h <- nrow(xmat)
  pacs$k <- ncol(xmat)
  pacs$s <- nrow(targets)
  pacs$xmat <- xmat
  pacs$wh <- wh
  pacs$targets <- targets
  pacs
}




