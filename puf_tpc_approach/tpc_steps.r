
# from start of step 8 in puf_tpc.r ----
ig <- 2

wh <- wts_adj %>%
  filter(AGI_STUB==ig) %>%
  .$wh_adj # the new national weight -- use either wh_sum or wh_sumadj

# target_vars <- default_targets
target_vars <- target_vars_df %>%
  filter(STATE=="AL", AGI_STUB==ig) %>%
  .$target_vec %>%
  unlist
target_vars <- setdiff(target_vars, "posAGI_nnz")

targets <- puf_targ %>%
  filter(AGI_STUB==ig) %>%
  select(AGI_STUB, STATE, all_of(target_vars))

# adjust these targets to get rid of any zero or NA values - for now, calculate an adjusted
# target for any zero or NA target so that its per-return value is equal to 50% of
# the lowest nonzero per-return value for other states for that target
# N1_nnz
per_return <- puf_targ %>%
  filter(AGI_STUB==ig) %>%
  select(AGI_STUB, STATE, all_of(target_vars)) %>%
  pivot_longer(-c(AGI_STUB, STATE, N1_nnz)) %>%
  mutate(per_ret=value / N1_nnz) %>%
  group_by(name) %>%
  mutate(isgood=!is.na(per_ret) & per_ret!=0,
         abspr=ifelse(isgood, abs(per_ret), NA),
         min_isgood=min(abspr, na.rm=TRUE),
         whichmin=which.min(abspr),
         min_nzval_per_ret=min_isgood * sign(value[whichmin]),
         value_adj=ifelse(is.na(value) | value==0,
                          N1_nnz * min_nzval_per_ret * .5,
                          value))
per_return %>%
  filter(is.na(value) | value==0)

per_return %>%
  filter(name=="A09600_sum")

per_return %>%
  filter(value != value_adj)

targets_adj <- per_return %>%
  select(AGI_STUB, STATE, N1_nnz, name, value_adj) %>%
  pivot_wider(values_from = value_adj) %>%
  select(AGI_STUB, STATE, all_of(target_vars)) # ensure that columns are in the original order

# after reviewing, use targets_adj instead of targets
targmat <- targets %>%
  select(all_of(target_vars)) %>%
  as.matrix
rownames(targmat) <- targets$STATE

targmat_adj <- targets_adj %>%
  select(all_of(target_vars)) %>%
  as.matrix
rownames(targmat_adj) <- targets$STATE

# compare the two
targmat
targmat_adj
targmat_adj - targmat


targets <- targmat_adj


xmat <- pufstrip %>%
  filter(AGI_STUB==ig) %>%
  select(all_of(target_vars)) %>%
  as.matrix


# now start the tpc approach ----
xpx <- t(xmat) %*% xmat
invxpx <- solve(xpx) # TODO: add error check and exit if not invertible


maxiter <- 100
beta0 <- matrix(0, nrow=nrow(targets), ncol=ncol(targets))
delta0 <- get_delta(wh, beta0, xmat) # tpc uses initial delta based on initial beta

ebeta <- beta0 # tpc uses 0 as beta starting point
edelta <- delta0 # tpc uses initial delta based on initial beta

sse_vec <- rep(NA_real_, maxiter)

step_scale <- 1
step_scale <- 1000
step_scale <- 10e3

step_inputs <- list()
step_inputs$targets <- targets
step_inputs$step_scale <- step_scale
step_inputs$xmat <- xmat
step_inputs$invxpx <- invxpx
step_inputs$wh <- wh

iter <- 1


iter <- iter + 1
# for(iter in 1:maxiter){
  # iter <- iter + 1
  edelta <- get_delta(wh, ebeta, xmat)
  ewhs <- get_weights(ebeta, edelta, xmat)
  ews <- colSums(ewhs)
  ewh <- rowSums(ewhs)
  step_inputs$ews <- ews

  etargets <- t(ewhs) %*% xmat
  d <- targets - etargets
  step_inputs$d <- d

  rel_err <- ifelse(targets==0, NA, abs(d / targets) * 100)
  max_rel_err <- max(rel_err, na.rm=TRUE)
  sse <- sum(d^2)
  # if(is.na(sse)) break # bad result, end it now, we have already saved the prior best result

  sse_vec[iter] <- sse
  sse_vec
  step_scale <- 4000
  # sse_vec <- c(seq(200, 100, -1), NA, NA)
  sse_rel_change <- sse_vec / lag(sse_vec) - 1
  # iter <- 5
  # test2 <- ifelse(iter >= 5, !any(sse_rel_change[iter - 0:2] < -.01), FALSE)
  # test2
  # any(sse_rel_change[c(4, 5, 6)] < -.01)

  best_sse <- min(sse_vec, na.rm=TRUE)
  if(sse==best_sse) best_ebeta <- ebeta
  prior_sse <- sse

  #if(iter <=20 | iter %% 20 ==0) print(sprintf("iteration: %i, sse: %.5e, max_rel_err: %.5e", iter, sse, max_rel_err))
  print(sprintf("iteration: %i, sse: %.5e, max_rel_err: %.5e", iter, sse, max_rel_err))

  #.. stopping criteria ---- iter <- 5
  #test1 <- max_rel_err < tol # every distance from target is within our desired error tolerance
  # test2: none the of last 3 iterations had sse improvement of 0.1% or more
  #test2 <- ifelse(iter >= 5, !any(sse_rel_change[iter - 0:2] < -.001), FALSE)

  # if(test1 | test2) {
  #   # exit if good
  #   print(sprintf("exit at iteration: %i, sse: %.5e, max_rel_err: %.5e", iter, sse, max_rel_err))
  #   break
  # }

  # if sse is > prior sse, adjust step scale downward
  # if(step_method=="adhoc" & (sse > best_sse)){
  #   step_scale <- step_scale * .5
  #   ebeta <- best_ebeta # reset and try again
  # }

  prior_ebeta <- ebeta

  # ad hoc step
  step <- -(1 / ews) * d %*% invxpx * step_scale
  #step_inputs$step_scale <- step_scale
  #step <- step_fn(ebeta, step_inputs) #  * (1 - iter /maxiter) # * step_scale # * (1 - iter /maxiter)
  # print(step)

  ebeta <- ebeta - step
}

best_edelta <- get_delta(ewh, best_ebeta, xmat)
ewhs <- get_weights(best_ebeta, best_edelta, xmat)
ewh <- rowSums(ewhs)
if(scale==TRUE) etargets <- sweep(etargets, 2, problem$scale_factor, "/")
final_step_scale <- step_scale

t2 <- proc.time()
total_seconds <- as.numeric((t2 - t1)[3])
















# solve poisson problem ----
solve_poisson <- function(problem, maxiter=100, scale=FALSE, scale_goal=1000, step_method="adhoc", step_scale=1, tol=1e-3, start=NULL){
  t1 <- proc.time()

  if(step_method=="adhoc") step_fn <- step_adhoc else{
    if(step_method=="finite_diff") step_fn <- step_fd
  }

  init_step_scale <- step_scale

  problem_unscaled <- problem
  if(scale==TRUE) problem <- scale_problem(problem, scale_goal)

  # unbundle the problem list and create additional variables needed
  targets <- problem$targets
  wh <- problem$wh
  xmat <- problem$xmat

  xpx <- t(xmat) %*% xmat
  invxpx <- solve(xpx) # TODO: add error check and exit if not invertible

  if(is.null(start)) beta0 <- matrix(0, nrow=nrow(targets), ncol=ncol(targets)) else # tpc uses 0 as beta starting point
    beta0 <- start
  delta0 <- get_delta(wh, beta0, xmat) # tpc uses initial delta based on initial beta

  ebeta <- beta0 # tpc uses 0 as beta starting point
  edelta <- delta0 # tpc uses initial delta based on initial beta

  sse_vec <- rep(NA_real_, maxiter)

  step_inputs <- list()
  step_inputs$targets <- targets
  step_inputs$step_scale <- step_scale
  step_inputs$xmat <- xmat
  step_inputs$invxpx <- invxpx
  step_inputs$wh <- wh

  for(iter in 1:maxiter){
    # iter <- iter + 1
    edelta <- get_delta(wh, ebeta, xmat)
    ewhs <- get_weights(ebeta, edelta, xmat)
    ews <- colSums(ewhs)
    ewh <- rowSums(ewhs)
    step_inputs$ews <- ews

    etargets <- t(ewhs) %*% xmat
    d <- targets - etargets
    step_inputs$d <- d

    rel_err <- ifelse(targets==0, NA, abs(d / targets))
    max_rel_err <- max(rel_err, na.rm=TRUE)
    sse <- sum(d^2)
    if(is.na(sse)) break # bad result, end it now, we have already saved the prior best result

    sse_vec[iter] <- sse
    # sse_vec <- c(seq(200, 100, -1), NA, NA)
    sse_rel_change <- sse_vec / lag(sse_vec) - 1
    # iter <- 5
    # test2 <- ifelse(iter >= 5, !any(sse_rel_change[iter - 0:2] < -.01), FALSE)
    # test2
    # any(sse_rel_change[c(4, 5, 6)] < -.01)

    best_sse <- min(sse_vec, na.rm=TRUE)
    if(sse==best_sse) best_ebeta <- ebeta
    prior_sse <- sse

    if(iter <=20 | iter %% 20 ==0) print(sprintf("iteration: %i, sse: %.5e, max_rel_err: %.5e", iter, sse, max_rel_err))

    #.. stopping criteria ---- iter <- 5
    test1 <- max_rel_err < tol # every distance from target is within our desired error tolerance
    # test2: none the of last 3 iterations had sse improvement of 0.1% or more
    test2 <- ifelse(iter >= 5, !any(sse_rel_change[iter - 0:2] < -.001), FALSE)

    if(test1 | test2) {
      # exit if good
      print(sprintf("exit at iteration: %i, sse: %.5e, max_rel_err: %.5e", iter, sse, max_rel_err))
      break
    }

    # if sse is > prior sse, adjust step scale downward
    if(step_method=="adhoc" & (sse > best_sse)){
      step_scale <- step_scale * .5
      ebeta <- best_ebeta # reset and try again
    }

    prior_ebeta <- ebeta

    # ad hoc step
    # step <- -(1 / ews) * d %*% invxpx * step_scale
    step_inputs$step_scale <- step_scale
    step <- step_fn(ebeta, step_inputs) #  * (1 - iter /maxiter) # * step_scale # * (1 - iter /maxiter)
    # print(step)

    ebeta <- ebeta - step
  }

  best_edelta <- get_delta(ewh, best_ebeta, xmat)
  ewhs <- get_weights(best_ebeta, best_edelta, xmat)
  ewh <- rowSums(ewhs)
  if(scale==TRUE) etargets <- sweep(etargets, 2, problem$scale_factor, "/")
  final_step_scale <- step_scale

  t2 <- proc.time()
  total_seconds <- as.numeric((t2 - t1)[3])

  keepnames <- c("total_seconds", "maxiter", "iter", "max_rel_err", "sse", "sse_vec", "d", "best_ebeta", "best_edelta", "ewh", "ewhs", "etargets",
                 "problem_unscaled", "scale", "scale_goal", "init_step_scale", "final_step_scale")
  result <- list()
  for(var in keepnames) result[[var]] <- get(var)
  print("all done")
  result
  # end solve_poisson
}

