
# delta, weights and vtom ----

get_delta <- function(wh, beta, xmat){
  # we cannot let beta %*% xmat get too large!! or exp will be Inf and problem will bomb
  # it will get large when a beta element times an xmat element is large, so either
  # beta or xmat can be the problem
  beta_x <- exp(beta %*% t(xmat))
  log(wh / colSums(beta_x)) # denominator is sum for each person
}


get_weights <- function(beta, delta, xmat){
  # get whs: state weights for households, given beta matrix, delta vector, and x matrix
  beta_x <- beta %*% t(xmat)
  # add delta vector to every row of beta_x and transpose
  beta_xd <- apply(beta_x, 1 , function(mat) mat + delta) 
  exp(beta_xd)
}


scale_problem <- function(problem, scale_goal){
  # problem is a list with at least the following:
  #  targets
  #  xmat
  # return:
  #   list with the scaled problem, including all of the original elements, plus
  #   scaled versions of x and targets
  #   plus new items scale_goal and scale_factor

  max_targets <- apply(problem$targets, 2, max) # find max target in each row of the target matrix
  scale_factor <- scale_goal / max_targets
  
  scaled_targets <- sweep(problem$targets, 2, scale_factor, "*")
  scaled_x <- sweep(problem$xmat, 2, scale_factor, "*")
  
  scaled_problem <- problem
  scaled_problem$targets <- scaled_targets
  scaled_problem$xmat <- scaled_x
  scaled_problem$scale_factor <- scale_factor
  
  scaled_problem
}


scale_problem_mdn <- function(problem, scale_goal){
  # problem is a list with at least the following:
  #  targets
  #  xmat
  # return:
  #   list with the scaled problem, including all of the original elements, plus
  #   scaled versions of x and targets
  #   plus new items scale_goal and scale_factor
  
  mdn_targets <- apply(problem$targets, 2, median) # find median target in each row of the target matrix
  scale_factor <- scale_goal / mdn_targets
  
  scaled_targets <- sweep(problem$targets, 2, scale_factor, "*")
  scaled_x <- sweep(problem$xmat, 2, scale_factor, "*")
  
  scaled_problem <- problem
  scaled_problem$targets <- scaled_targets
  scaled_problem$xmat <- scaled_x
  scaled_problem$scale_factor <- scale_factor
  
  scaled_problem
}


vtom <- function(vec, nrows){
  # vector to matrix in the same ordering as a beta matrix
  matrix(vec, nrow=nrows, byrow=FALSE)
}


# differences and sse ----

diff_vec <- function(betavec, wh, xmat, targets){
  # return a vector of differences between targets and corresponding
  # values calculated given a beta vector, household weights, and x matrix
  beta <- vtom(betavec, nrows=nrow(targets))
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  as.vector(d)
}


sse_fn <- function(betavec, wh, xmat, targets){
  # return a single value - sse (sum of squared errors)
  sse <- sum(diff_vec(betavec, wh, xmat, targets)^2)
  sse
}


# step functions ----
step_fd <- function(ebeta, step_inputs){
  # finite differences
  bvec <- as.vector(ebeta)
  gbeta <- numDeriv::grad(sse_fn, bvec,  wh=step_inputs$wh, xmat=step_inputs$xmat, targets=step_inputs$targets)
  hbeta <- numDeriv::hessian(sse_fn, bvec,  wh=step_inputs$wh, xmat=step_inputs$xmat, targets=step_inputs$targets)
  ihbeta <- solve(hbeta)
  stepvec <- t(ihbeta %*% gbeta)
  step <- vtom(stepvec, nrows=nrow(ebeta))
  step
}

step_fd <- function(ebeta, step_inputs){
  # finite differences -- version to print time
  bvec <- as.vector(ebeta)
  t1 <- proc.time()
  gbeta <- numDeriv::grad(sse_fn, bvec,  wh=step_inputs$wh, xmat=step_inputs$xmat, targets=step_inputs$targets)
  t2 <- proc.time()
  print(sprintf("gradient time in seconds: %.1e", (t2-t1)[3]))
  hbeta <- numDeriv::hessian(sse_fn, bvec,  wh=step_inputs$wh, xmat=step_inputs$xmat, targets=step_inputs$targets)
  t3 <- proc.time()
  print(sprintf("hessian time in seconds: %.1e", (t3-t2)[3]))
  ihbeta <- solve(hbeta)
  t4 <- proc.time()
  print(sprintf("inverse time in seconds: %.1e", (t4-t3)[3]))
  stepvec <- t(ihbeta %*% gbeta)
  step <- vtom(stepvec, nrows=nrow(ebeta))
  step
}


step_adhoc <- function(ebeta, step_inputs){
  -(1 / step_inputs$ews) * step_inputs$d %*% step_inputs$invxpx * step_inputs$step_scale
}

get_step <- function(step_method, ebeta, step_inputs){
  step <- case_when(step_method=="adhoc" ~ step_adhoc(ebeta, step_inputs),
                    step_method=="finite_diff" ~ step_fd(ebeta, step_inputs),
                    TRUE ~ step_adhoc(ebeta, step_inputs))
  step
}

jac <- function(ewhs, xmatrix){
  # jacobian of distance vector relative to beta vector, IGNORING delta
  x2 <- xmatrix * xmatrix
  ddiag <- - t(ewhs) %*% x2 # note the minus sign in front
  diag(as.vector(ddiag)) 
}


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
}
