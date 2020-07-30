
#****************************************************************************************************
#                NEW: constraint setup ####
#****************************************************************************************************

get_ccdf <- function(pufdf, wtvar, voi){
  # get the constraint coefficients for all "n" and all "value" constraints for the variables of interest
  # for each state and income range
  ccdf <- pufdf %>%
    select(stabbr, incrange, one_of(c(wtvar, voi))) %>%
    mutate(uid=row_number(), wt_toadj=get(wtvar)) %>%
    gather(variable, value, -uid, -stabbr, -incrange, -wt_toadj) %>%
    filter(value!=0) %>%
    mutate(n=wt_toadj, value=wt_toadj * value) %>%
    gather(cctype, value, value, n) %>%
    mutate(cname=paste(stabbr, incrange, variable, cctype, sep="_"))
  return(ccdf)
}

get_ccoef <- function(ccrules, ccdf){
  ccoef <- ccrules %>%
    left_join(ccdf %>% select(cname, j=uid, wt_toadj, value), by="cname") %>%
    arrange(cname) %>%
    mutate(i=group_indices(., cname)) %>%
    select(cname, i, j, value, stabbr, incrange, soivar, cctype) %>%
    arrange(i, j)
}


apply_tol_rules <- function(constraint_rhs, tol_rules){
  get_clb_cub <- function(constraint_rhs){
    constraint_rhs <- constraint_rhs %>%
      mutate(clb.unscaled=ifelse(is.infinite(tol), -Inf, target.unscaled - abs(tol * target.unscaled)),
             cub.unscaled=ifelse(is.infinite(tol), +Inf, target.unscaled + abs(tol * target.unscaled)))
    return(constraint_rhs)
  }
  
  constraint_rhs$tol <- 0
  
  for(i in 1:nrow(tol_rules)){
    rule <- parse(text=tol_rules$rule[i])
    tolval <- tol_rules$tol[i] %>% as.numeric
    
    if(length(rule)==1){ # this test makes sure that commented-out rules don't cause an error
      constraint_rhs <- constraint_rhs %>%
        mutate(toltemp=ifelse(with(., eval(rule)), tolval, 0),
               tol=pmax(tol, toltemp)) %>%
        select(-toltemp)
    }
  }
  constraint_rhs <- get_clb_cub(constraint_rhs)
  
  return(constraint_rhs)
}



#****************************************************************************************************
#                Expand constraint rules ####
#****************************************************************************************************
expand_incranges <- function(rulesdf){
  if(rulesdf$inc_ranges=="all"){
    rulesdf2 <- rulesdf %>% 
      right_join(rangemap %>% 
                   filter(incrange %in% c("lt50k", "ge50lt75k", "ge75lt100k", "ge100lt200k", "ge200k")) %>% 
                   mutate(inc_ranges="all"), 
                 by="inc_ranges")
  } else rulesdf2 <- rulesdf %>% 
      mutate(incrange=inc_ranges) %>%
      left_join(rangemap) # get inc_rule
  return(rulesdf2)
}


expand_states <- function(rulesdf){
  if(rulesdf$states=="all"){
    rulesdf2 <- rulesdf %>% 
      right_join(tibble(stabbr=stdc, states="all"), by="states")
  } else rulesdf2 <- rulesdf %>% mutate(stabbr=states)
  return(rulesdf2)
}


get_constraint_rules_all <- function(rulesdf){
  # expand the "all" rules
  # input: rulesdf with columns:
  #   soivar, states, inc_ranges, value_calc
  rules2 <- rulesdf %>%
    filter(!is.na(soivar)) %>%
    group_by(recnum=row_number()) %>%
    do(expand_incranges(.)) %>%
    ungroup %>% # we must reset recnum
    group_by(recnum=row_number()) %>%
    do(expand_states(.)) %>%
    ungroup %>%
    select(-recnum)
  
  # now make the final constraint rules
  constraint_rules_all <- rules2 %>%
    mutate(state_rule=paste0("(stabbr=='", stabbr, "')"),
           cname=paste(soivar, stabbr, incrange, value_calc, sep="_"),
           full_rule=paste(inc_rule, state_rule, sep=" & ")) %>%
    select(cname, full_rule, value_calc, stabbr, soivar, incrange, inc_ranges, states)
  return(constraint_rules_all)
}


#****************************************************************************************************
#                # General function for creating constraints ####
#****************************************************************************************************

con.fg <- function(st, vname, puf.df, sts, ids){
  # This function fills part of an imaginary matrix with constraint coefficients, where a coefficient is
  # the amount by which a constraint value changes if we change the x for a given person in a given state.
  # Rows are constraints and columns are person-state combinations.
  # Indexes:
  #   i for rows -- 1 per constraint
  #   j for columns (people-state combinations)
  
  # If a problem has 1,000 people and 6 states, and 5 constraints per state, the imaginary
  # matix will have 6,000 columns and 1,030 rows (1 adding-up constraint per person where the sum of their x values
  # must equal 1, plus 5 constraints per state)
  
  # The columns are ordered so that we have all of the columns for person 1, then all for person 2, and so on.
  # Thus, if the 6 states are AK, AL, AR, AZ, CA, and CO, then column 3 is has the coefficients for person 1 in AR,
  # and column 8 has the coefficients for person 2 in AL (person 2 starts in column 7).
  
  # Example: if person 2's agi is $10,000 and their weight is 20, and if constraint # 1,003 is
  # a target for agi in AL, then the constraint coefficient for person 2 in AL for this constraint is 200e3. If
  # we increase person 2's x in AL by 1, the amount of agi in AL increases by $200k. We would store this
  # coefficient in row 1,003, column 8 of the imaginary matrix (i.e., its indexes are 1003, 8).
  
  # Most of the cells in this matrix will be zero. In the example above, person 2's coefficients start in column 7 and
  # run through column 12, storing coefficients for AK, AL, AR, AZ, CA, and CO. If row 1,003 is the agi constraint for AL,
  # then the value in row 1003, column 8 is 200e3 but the values in this row and columns 7 and 9-12 are zero - if we change
  # the x value for person 2 in AL, it does not change person 2's agi in AK, AR, or the other states; it only changes in AL.
  
  # We store this imaginary matrix in a compact way, only keeping the non-zero coefficients, and their row and column 
  # indexes. We also keep a name for the constraint and the record id of the person, for convenience
  
  # This function calculates the nonzero elements for one state column for each person, for a particular constraint. 
  # It is stored in a single row of this imaginary matrix. We get the column indexes now. We will determine the 
  # row index later, when we know which constraint is in which row.
  
  # Normally, we would call this function one time for each state-column.
  
  stnum <- which(sts == st) # what is the number of the state to which this state-specific index applies?
  
  colidx <- seq(from=stnum, by=length(sts), length.out=length(ids)) # get each person's column index for this state
  
  df <- tibble(cname=paste0(vname, ".", st),
               RECID=puf.df$RECID,
               j=as.integer(colidx)) %>%
    mutate(value=puf.df$wt * puf.df[[vname]]) %>%
    filter(value!=0) # in case a coeff is zero, drop it
  return(df)
}

  
get_a_constraint_for_each_state <- function(cname, puf.df, sts, ids){
  # get the constraint vname for each state in sts by calling
  # con.fg2 repeatedly (one time for each state)
  print(paste0("Getting constraint: ", cname))
  df <- ldply(sts, con.fg, cname, puf.df, sts, ids)
  return(df)
}

get_all_state_constraints <- function(cnames, puf.df, sts, ids){
  # get all of the constraints named in cnames by calling get_a_constraint_for_each_state
  # repeatedly, once for each cname in cnames
  print(paste("Getting constraints: ", paste(cnames, collapse=", ")))
  df <- ldply(cnames, get_a_constraint_for_each_state, puf.df, sts, ids)
  return(df)
}


#****************************************************************************************************
#                # Constraint coefficients, rhs, and scaling ####
#****************************************************************************************************
get_ccmatdf <- function(puf.df, cnames, sts){
  # Create sparse constraint coefficients matrix, in data frame form
  # create needed vectors
  ids <- puf.df$RECID %>% sort
  
  # adding-up constraints
  print("....getting the adding-up constraints for each return (sum of shares allocated to states must equal 1)")
  # xij sum for each person i across states j must be 1 -- this is dense
  ccmat.ids.df <- expand.grid(RECID=ids, 
                              stabbr=sts) %>%
    mutate(stabbr=as.character(stabbr)) %>%
    arrange(RECID, stabbr) %>%
    mutate(cname=paste0("sum.", RECID),
           j=row_number(),
           value=1)
  
  state_constraints <- get_all_state_constraints(cnames, puf.df, sts, ids)

  # put the constraints together
  ccmat.df <- bind_rows(ccmat.ids.df, 
                        state_constraints) %>%
    mutate(i=match(cname, unique(cname))) %>%
    select(RECID, cname, i, j, value)
  
  return(ccmat.df)
}


get_rhs <- function(nrecs, pufsumdf, cnames, rhs_shares=NULL){
  #  Create RHS of constraints
  if(is.null(rhs_shares)) rhs.sharetotals <- rep(1, nrecs) else
    rhs.sharetotals <- rhs_shares
  
  # rhsvars <- ifelse(cnames=="wt_one", "wt", cnames)
  rhsvars <- cnames
  
  targets <- c(rhs.sharetotals,
               stack(pufsumdf[, rhsvars])$values)
  return(targets)
}

scale_median <- function(ccmat.df, scale, mdnval=50){
  # scale so that the median of the nonzero constraint coefficients for each constraint is 50
  # get the median value of each raw constraint coefficient
  mdn <- ccmat.df %>% 
    group_by(i) %>% # group by each constraint
    filter(value!=0) %>% # value is the constraint coefficient for each j (variable), for this constraint (i)
    summarise(value=median(value)) %>%
    .[["value"]] # result is a vector
  
  scale <- ifelse(mdn!=0, (1 / mdn) * mdnval, scale)
  return(scale)
}

scale_2norm <- function(ccmat.df, scale){
  # Modeling and Optimization in Space Engineering, edited by Giorgio Fasano, János D. Pintér, pp.51-52
  # scale so that the 2-norm  of each constraint is 1, where 2-norm is sqrt of sum of squared coefficients
  norm_vec <- function(vec) sqrt(sum(vec^2))
  
  norm2 <- ccmat.df %>% 
    group_by(i) %>% # group by each constraint
    filter(value!=0) %>% # value is the constraint coefficient for each j (variable), for this constraint (i)
    summarise(value=norm_vec(value)) %>%
    .[["value"]] # result is a vector
  
  scale <- ifelse(norm2!=0, (1 / norm2), scale)
  return(scale)
}


scale_constraints <- function(ccmat.df, targets, scale_type="none"){
  # define a scale value for each constraint and then apply
  
  scale <- rep(1, length(targets)) # default scaling is 1
  
  if(scale_type=="none"){
    # do nothing
  } else if(scale_type=="median50") {
    scale <- scale_median(ccmat.df, scale)
    
  } else if(scale_type=="2norm") {
    scale <- scale_2norm(ccmat.df, scale)
  }
  
  # print(ht(scale))
  
  targets.scaled <- targets * scale
  
  ccmat.scaled <- ccmat.df %>%
    mutate(value = value * scale[i]) # use i in ccmat.df as index into scale
  
  constraints.scaled <- list()
  constraints.scaled$targets.scaled <- targets.scaled
  constraints.scaled$ccmat.scaled <- ccmat.scaled
  
  return(constraints.scaled)
}


#****************************************************************************************************
#                # More constraint items ####
#****************************************************************************************************

getcc <- function(constraint_rules, pufdf, progress=FALSE){
  # the puf data frame MUST have the variable uid -- a UNIQUE id for each RECORD
  
  get_a_cc <- function(i){
    # get a single constraint coefficient for constraint i
    full_rule <- parse(text=constraint_rules$full_rule[i])
    value_calc <- parse(text=constraint_rules$value_calc[i])
    cname <- constraint_rules$cname[i]
    
    df <- pufdf %>% filter(eval(with(., full_rule))) %>%
      mutate(cname=cname,
             value=eval(with(., value_calc))) %>%
      filter(value!=0) %>%
      select(cname, j=uid, stabbr, incrange, value.unscaled=value)
    return(df)
  }
  
  # get coefficients for each of the constraints in constraint_rules
  progtext <- "none"
  if(progress) progtext <- "text"
  concoef <- ldply(1:nrow(constraint_rules), get_a_cc, .progress=progtext) %>%
    mutate(i=group_indices(., cname)) %>%
    arrange(i, j)
  return(concoef)
}


pufsum <- function(constraint_rules, pufdf) {
  getcon <- function(i){
    full_rule <- parse(text=constraint_rules$full_rule[i])
    value_calc <- parse(text=constraint_rules$value_calc[i])
    cname <- constraint_rules$cname[i]
    df2 <- tibble(cname=cname, initvalue=sum(with(pufdf, eval(full_rule) * eval(value_calc))))
    return(df2)
  }
  
  consums <- ldply(1:nrow(constraint_rules), getcon) # , .progress="text"
  return(consums)
}


calcsums <- function(df, wtvar){
  dfsum <- df %>%
    mutate(wt_n=1,
           mars2_n=((MARS==2)*1),
           netcg_n=(netcg!=0)*1,
           E18425_n=(E18425!=0)*1, # state-local income tax
           E18500_n=(E18500!=0)*1, # state-local real estate tax
           wt_calc=.[[wtvar]]) %>%
    group_by(stabbr, incrange) %>%
    summarise_at(vars(wt_n, mars2_n, E00100, E00200,  # weight, married joint, agi, wages
                      netcg, netcg_n, 
                      E18425, E18425_n, # state-local income tax
                      E18500, E18500_n, # state-local real estate tax
                      E04800, E06500), # taxable income, fed income tax
                 funs(sum(. * wt_calc)))
  return(dfsum)
}


#****************************************************************************************************
#                # Starting points ####
#****************************************************************************************************
x0_wtshare <- function(inputs){
  wtshares_vec <- inputs$pufsum$wt_one / sum(inputs$pufsum$wt_one)
  names(wtshares_vec) <- inputs$pufsum$stabbr
  x0m <- matrix(data=wtshares_vec, nrow=inputs$nrecs, ncol=length(wtshares_vec), 
                byrow=TRUE, dimnames=list(NULL, inputs$sts)) # recycle wtshares_vec
  x0 <- xij_to_x(x0m) # convert back to vector
  return(x0)
}


x0_nstates <- function(inputs) {
  x0 <- rep(1 / length(inputs$sts), inputs$num.vars)
  return(x0)
}



#****************************************************************************************************
#                # Problem definition ####
#****************************************************************************************************

get_problem <- function(puf.df, cnames, rhs_targets, rhs_shares=NULL){
  print("Defining problem list")
  problem <- list()
  problem$sts <- rhs_targets$stabbr
  problem$cnames <- cnames
  problem$nrecs <- nrow(puf.df)
  problem$pufsum <- rhs_targets
  
  print("..getting unscaled RHS targets")
  problem$targets.unscaled <- get_rhs(problem$nrecs, problem$pufsum, problem$cnames, rhs_shares)
  
  print("..getting unscaled constraint coefficients, may take some time")
  problem$ccmat.df.unscaled <- get_ccmatdf(puf.df, problem$cnames, problem$sts)
  
  problem$num.vars <- max(problem$ccmat.df.unscaled$j)
  problem$num.constraints <- max(problem$ccmat.df.unscaled$i)
  problem$nnz <- nrow(problem$ccmat.df.unscaled)
  
  print("..getting indexes of nonzero constraint coefficients, to define sparsity structure")
  problem$ccmat.indexes <- dlply(problem$ccmat.df.unscaled, "i", function(x) return(x$j))
  
  print("Done defining problem")
  cat("\n")
  
  return(problem)
}


get_inputs_list <- function(problem, scale_type="none"){
  # scale_type: none, median50, 2norm
  print("Creating inputs list for Ipopt")
  inputs <- list()
  inputs$eps <- .0001
  inputs$fmult <- 1
  inputs$sts <- problem$sts
  inputs$nrecs <- problem$nrecs
  inputs$num.vars <- problem$num.vars
  inputs$num.constraints <- problem$num.constraints
  inputs$nnz <- problem$nnz
  inputs$pufsum <- problem$pufsum
  
  # scale constraint coefficients and targets even though we may not use them
  print("..scaling constraint coefficients and targets if requested")
  inputs$scaling <- scale_type
  all.constraints.scaled <- scale_constraints(problem$ccmat.df.unscaled,
                                              problem$targets.unscaled,
                                              scale_type=inputs$scaling)
  inputs$targets  <- all.constraints.scaled$targets.scaled
  inputs$ccmat.df <- all.constraints.scaled$ccmat.scaled

  inputs$eval_jac_g_structure <- problem$ccmat.indexes
  
  print("..Defining sparsity structure of the Hessian")
  inputs$eval_h_structure <- lapply(1:problem$num.vars, function(x) x) # diagonal elements of our Hessian
  
  print("Done creating inputs list")
  cat("\n")
  return(inputs)
}



#****************************************************************************************************
#                # Get constraints-comparison data frame ####
#****************************************************************************************************

get_constraints_df <- function(inputs, x0){
  dftarg <- tibble(cname=unique(inputs$ccmat.df$cname),
                   con_target=inputs$targets,
                   con_x0=eval_g(x0, inputs),
                   diff=con_x0 - con_target,
                   pdiff=diff / con_target * 100)
  return(dftarg)
}


#****************************************************************************************************
#                # Problem solution - entropy maximization ####
#****************************************************************************************************
solve_problem_entmax <- function(inputs, opts, x0=NULL, xlb=NULL, xub=NULL, ctol=0){
  clb <- ifelse(inputs$targets >= 0, 
                inputs$targets * (1 - ctol),
                inputs$targets * (1 + ctol))
  cub <- ifelse(inputs$targets >= 0, 
                inputs$targets * (1 + ctol),
                inputs$targets * (1 - ctol))
  
  # Create starting values
  if(is.null(x0)) x0 <- rep(1 / length(inputs$sts), inputs$num.vars)

  if(is.null(xlb)) xlb <- rep(0, length(x0))
  if(is.null(xub)) xub <- rep(1, length(x0))
  
  a <- proc.time()
  res <- ipoptr(x0 = x0,
                lb = xlb,
                ub = xub,
                eval_f = eval_f_entmax, 
                eval_grad_f = eval_grad_f_entmax, 
                eval_g = eval_g, 
                eval_jac_g = eval_jac_g,
                eval_jac_g_structure = inputs$eval_jac_g_structure,
                eval_h = eval_h_entmax,
                eval_h_structure = inputs$eval_h_structure,
                constraint_lb = clb,
                constraint_ub = cub,
                opts = opts,
                inputs = inputs)
  b <- proc.time()
  print((b - a) / 60)
  return(res)
}


# inputs <- inputs0; ctol <- .02
solve_problem_xmx2 <- function(inputs, opts, x0=NULL, xlb=NULL, xub=NULL, ctol=0){
  clb <- ifelse(inputs$targets >= 0, 
                inputs$targets * (1 - ctol),
                inputs$targets * (1 + ctol))
  cub <- ifelse(inputs$targets >= 0, 
                inputs$targets * (1 + ctol),
                inputs$targets * (1 - ctol))
  
  # Create starting values
  if(is.null(x0)) x0 <- rep(0, inputs$num.vars)
  
  if(is.null(xlb)) xlb <- rep(0, inputs$num.vars)
  if(is.null(xub)) xub <- rep(1, inputs$num.vars)
  
  a <- proc.time()
  res <- ipoptr(x0 = x0,
                lb = xlb,
                ub = xub,
                eval_f = eval_f_xmx2, 
                eval_grad_f = eval_grad_f_xmx2, 
                eval_g = eval_g, 
                eval_jac_g = eval_jac_g,
                eval_jac_g_structure = inputs$eval_jac_g_structure,
                #eval_h = eval_h_xmx2,
                #eval_h_structure = inputs$eval_h_structure,
                constraint_lb = clb,
                constraint_ub = cub,
                opts = opts,
                inputs = inputs)
  b <- proc.time()
  print((b - a) / 60)
  return(res)
}

# length(eval_h_lown(x0, 1, rep(1, inputs$num.constraints)))
# length(unlist(inputs$eval_h_structure)) # 100,000
# names(inputs)



