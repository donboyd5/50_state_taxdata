
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

etargs_mat <- function(betavec, wh, xmat, s){
  # return a vector of calculated targets and corresponding
  # values calculated given a beta vector, household weights, and x matrix
  beta <- vtom(betavec, nrows=s)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  etargets
}


etargs_vec <- function(betavec, wh, xmat, s){
  # return a vector of calculated targets and corresponding
  # values calculated given a beta vector, household weights, and x matrix
  beta <- vtom(betavec, nrows=s)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  as.vector(etargets)
}


diff_vec <- function(betavec, wh, xmat, targets, dweights=rep(1, length(targets))){
  # return a vector of differences between targets and corresponding
  # values calculated given a beta vector, household weights, and x matrix
  beta <- vtom(betavec, nrows=nrow(targets))
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  d <- targets - etargets
  as.vector(d) * dweights
}


sse_fn <- function(betavec, wh, xmat, targets, dweights=NULL){
  # return a single value - sse (sum of squared errors)
  sse <- sum(diff_vec(betavec, wh, xmat, targets, dweights)^2)
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

jac_tpc <- function(ewhs, xmatrix){
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
  # end solve_poisson
}


grad_sse <- function(betavec, wh, xmat, targets){
  # return gradient of the sse function wrt each beta
  
  # get the deltas as we will need them
  beta <- vtom(betavec, nrows=nrow(targets))
  delta <- get_delta(wh, beta, xmat)
  
  # make a data frame for each relevant variable, with h, k, and/or s indexes as needed
  h <- nrow(xmat)
  k <- ncol(xmat)
  s <- nrow(targets)
  
  hstub <- tibble(h=1:h)
  skstub <- expand_grid(s=1:s, k=1:k) %>% arrange(k, s)
  hkstub <- expand_grid(h=1:h, k=1:k) %>% arrange(k, h)
  
  diffs <- diff_vec(betavec, wh, xmat, targets)
  diffsdf <- skstub %>% mutate(diff=diffs)
  
  whdf <- hstub %>% mutate(wh=wh)
  
  xdf <- hkstub %>% mutate(x=as.vector(xmat))
  
  etargets <- etargs_vec(betavec, wh, xmat, s)
  
  targdf <- skstub %>% mutate(target=as.vector(targets))
  etargdf <- skstub %>% mutate(etarget=etargets)
  betadf <- skstub %>% mutate(beta=betavec)
  deltadf <- hstub %>% mutate(delta=delta) 
  
  # now that the data are set up we are ready to calculate the gradient of the sse function
  # break the calculation into pieces using first the chain rule and then the product rule

  # sse = f(beta) = sum over targets [s,k] of (target - g(beta))^2
  #   where g(beta[s,k]) = sum over h(ws[h] * x[h,k]) and ws[h] is the TPC formula
  
  # chain rule for grad, for each beta[s,k] (where gprime is the partial of g wrt beta[s,k]):
  # = 2 * (target - g(beta)) * gprime(beta)
  # = 2 * diffs * gprime(beta[s,k])
  
  # for a single target[s,k]:
  #   g(beta)= weighted sum of X over hh, or sum[h] of X * exp(beta %*% X + delta[h]) where delta[h] is a function of all beta[s,k]
  #     Re-express:
  #          = sum[h] of X * exp(beta*X) * exp(delta[h]) # involving just the beta and x needed for this target
  #      
  #          = sum[h]  of X * a * b where a=exp(beta *X) and b=exp(delta[h]) and delta is a function of beta
  
  # product rule, still for a single target, gives:
  #     gprime(beta)=sum[h] of X * (a * bprime + b * aprime)
  
  # a = exp(beta * x)
  adf_base <- xdf %>%
    left_join(betadf, by="k") %>%
    mutate(a_exponent=beta * x) %>%
    select(h, s, k, x, beta, a_exponent)
  # adf_base %>% filter(h==1)
  
  adf <- adf_base %>%
    group_by(s, h) %>%
    summarise(a=exp(sum(a_exponent)), .groups="drop") %>% # these are the state weights for each hh BEFORE delta impact
    select(h, s, a)
  # adf %>% filter(h==1)
  
  #    aprime:
  #      deriv of a=exp(beta * x) wrt beta is = x * a
  # this is how much each hh's state weight will change if a beta changes, all else equal, BEFORE delta impact
  aprimedf <- adf %>%
    left_join(xdf, by="h") %>%
    mutate(aprime=x * a) %>%
    select(h, s, k, x, a, aprime) %>%
    arrange(h, s, k)
  # aprimedf %>% filter(h==1)
  
  # b = exp(delta[h])
  bdf <- deltadf %>%
    mutate(b=exp(delta))
  
  # check b: -- good
  # bcheck <- adf %>%
  #   left_join(whdf, by = "h") %>%
  #   group_by(h) %>%
  #   summarise(wh=first(wh), a=sum(a), .groups="drop") %>%
  #   mutate(bcheck=wh / a)
  
  # bprimedf # do this for each hh for each target I think
  # this is how much the delta impact will change if we change a beta - thus we have 1 per h, s, k
  # delta =log(wh/log(sum[s] exp(betas*X))
  # b=exp(delta(h))
  # which is just b = wh / log(sum[s] exp(betas*X))
  # bprime= for each h, for each beta (ie each s-k combination): (according to symbolic differentiation checks):
  #   bprime =  - (wh * xk *exp(bs) / sum[s] exp(bs))^2  where bs is exp(BX) for just that S and just that h
  # note that this bs is the same as a above: the sum, for an s-h combo, of exp(BX)
  
  # for each h, get the sum of their exp(beta * x) as it is a denominator; this is in adf
  # adf %>% filter(h==1)
  asums <- adf %>%
    group_by(h) %>%
    summarise(asum=sum(a), .groups="drop")
  
  bprimedf_base <- adf %>%
    left_join(whdf, by = "h") %>%
    left_join(xdf, by="h") %>%
    left_join(asums, by="h") %>%
    select(h, s, k, wh, x, a, asum)
  # bprimedf_base %>% filter(h==1)
  
  bprimedf <- bprimedf_base %>%
    mutate(bprime= -(wh * x * a / (asum^2)))
  # bprimedf %>% filter(h==1)
  
  # now get gprime:
  #     gprime(beta)=sum[h] of X * (a * bprime + b * aprime)
  #  bprimedf has most of what we need
  gprime_h <- bprimedf %>%
    select(-a) %>% # drop a as it is also in aprime, below
    left_join(bdf, by="h") %>%
    left_join(aprimedf %>% select(-x), by = c("h", "s", "k")) %>% # drop x as it is in bprime
    mutate(gprime_h=x * (a * bprime + b * aprime))
  # gprime_h %>% filter(h==1)
  
  gprime <- gprime_h %>%
    group_by(s, k) %>%
    summarise(gprime=sum(gprime_h), .groups="drop") # sum over the households h
  
  # put it all together to get the gradient by s, k
  # 2 * (target - g(beta)) * gprime(beta)    
  graddf <- diffsdf %>%
    left_join(gprime, by = c("s", "k")) %>%
    mutate(grad=2 * diff * gprime) %>%
    arrange(k, s)
  # graddf <- targdf %>%
  #   left_join(etargdf, by = c("s", "k")) %>%
  #   left_join(gprime, by = c("s", "k")) %>%
  #   mutate(term1=2 * target * gprime,
  #          term2=2 * etarget * gprime)
  
  
  graddf$grad
}


grad_sse_v2 <- function(betavec, wh, xmat, targets){
  # return gradient of the sse function wrt each beta
  
  # get the deltas as we will need them
  beta <- vtom(betavec, nrows=nrow(targets))
  delta <- get_delta(wh, beta, xmat)
  
  # make a data frame for each relevant variable, with h, k, and/or s indexes as needed
  h <- nrow(xmat)
  k <- ncol(xmat)
  s <- nrow(targets)
  
  hstub <- tibble(h=1:h)
  skstub <- expand_grid(s=1:s, k=1:k) %>% arrange(k, s)
  hkstub <- expand_grid(h=1:h, k=1:k) %>% arrange(k, h)
  
  diffs <- diff_vec(betavec, wh, xmat, targets)
  etargets <- etargs_vec(betavec, wh, xmat, s)
  diffsdf <- skstub %>% 
    mutate(target=as.vector(targets),
           etarget=etargets,
           diff=diffs)
  
  whdf <- hstub %>% mutate(wh=wh)
  
  xdf <- hkstub %>% mutate(x=as.vector(xmat))
  
  betadf <- skstub %>% mutate(beta=betavec)
  deltadf <- hstub %>% mutate(delta=delta) 
  
  # now that the data are set up we are ready to calculate the gradient of the sse function
  # break the calculation into pieces using first the chain rule and then the product rule
  
  # sse = f(beta) = sum over targets [s,k] of (target - g(beta))^2
  #   where g(beta[s,k]) = sum over h(ws[h] * x[h,k]) and ws[h] is the TPC formula
  
  # for each target, chain rule for grad, 
  # for each beta[s,k] (where gprime is the partial of g wrt beta[s,k]):
  # = - 2 * (target - g(beta)) * gprime(beta)
  # = - 2 * diffs * gprime(beta[s,k])
  
  # for a single target[s,k]:
  #   g(beta)= weighted sum of X over hh, or sum[h] of X * exp(beta %*% X + delta[h]) where delta[h] is a function of all beta[s,k]
  
  #     Re-express g(beta):
  #          = sum[h] of X * exp(beta*X) * exp(delta[h]) # involving just the beta and x needed for this target
  #      
  #          = sum[h]  of X * a * b where a=exp(beta *X) and b=exp(delta[h]) and delta is a function of beta
  
  # gprime(beta), still for a single target -- product rule gives:
  #     gprime(beta)=sum[h] of X * (a * bprime + b * aprime)
  
  # a = exp(beta * x)  - this is the same for all
  adf_base <- xdf %>%
    left_join(betadf, by="k") %>%
    mutate(a_exponent=beta * x) %>% # we will raise e to this power
    select(h, s, k, x, beta, a_exponent)
  # adf_base %>% filter(h==1)
  
  adf <- adf_base %>%
    group_by(s, h) %>% # weights are specific to state and household
    # we sum the k elements of the exponent and then raise e to that power
    summarise(a=exp(sum(a_exponent)), .groups="drop") %>% # these are the state weights for each hh BEFORE delta impact
    select(h, s, a)
  # adf %>% filter(h==1)
  
  #    aprime:
  #      deriv of a=exp(beta * x) wrt beta is = x * a
  # this is how much each hh's state weight will change if a beta changes, all else equal, BEFORE delta impact
  # since there is a beta for each s, k combination this will vary with different x[h, k] values
  aprimedf <- adf %>%
    left_join(xdf, by="h") %>%
    mutate(aprime=x * a) %>%
    select(h, s, k, x, a, aprime) %>%
    arrange(h, s, k)
  # aprimedf %>% filter(h==1)
  
  # b = exp(delta[h])
  bdf <- deltadf %>%
    mutate(b=exp(delta))
  
  # bprime - the hardest part -- how much does delta for an h change wrt a change in any beta
  # this is how much the delta impact will change if we change a beta - thus we have 1 per h, s, k
  # delta =log(wh/log(sum[s] exp(betas*X))
  # b=exp(delta(h))
  # which is just b = wh / log(sum-over-s]: exp(beta-for-given-s * x))
  
  # bprime= for each h, for each beta (ie each s-k combination):
  
  # from symbolic differentiation we have:
  # .e3 <- exp(b.s1k1 * x1 + b.s1k2 * x2) # which is a for a specific state -- s1 in this case
  # .e9 <- .e3 + exp(b.s2k1 * x1 + b.s2k2 * x2) + exp(b.s3k1 * x1 + b.s3k2 * x2)
  # so e9 <- exp(b.s1k1 * x1 + b.s1k2 * x2) + exp(b.s2k1 * x1 + b.s2k2 * x2) + exp(b.s3k1 * x1 + b.s3k2 * x2)
  # which is just asum as calculated below
  #  -(wh * x1 * .e3/(.e9 * log(.e9)^2))
  
  #  .e3 <- exp(b.s1k1 * x1 + b.s1k2 * x2)
  #   -(wh * x1 * .e3 / (.e3 + exp(b.s2k1 * x1 + b.s2k2 * x2) + exp(b.s3k1 * x1 + b.s3k2 * x2))^2)
  
  # which in a, asum notation is:
  #  -(wh * x1 * a / (asum)^2)
  
  # in my a, asum notation below, we have
  # deriv wrt b.s1k1 =
  #    -(wh * x[k=1] * a[s=1] / {asum^2})
  
  # for each h, get the sum of their exp(beta * x) as it is a denominator; this is in adf
  # log(wh / colSums(beta_x)) # denominator is sum for each person
  
  # asum is, for each h, exp(beta-s1k1 +...+ beta-s1kn) + ...+ exp(beta-smk1 + ...+ beta-smkn)
  # each a that we start with is exp(.) for one of the states so this the sum of exp(.) over the states
  asums <- adf %>%
    group_by(h) %>%
    summarise(asum=sum(a), .groups="drop")
  
  bprimedf_base <- adf %>%
    left_join(whdf, by = "h") %>%
    left_join(xdf, by="h") %>%
    left_join(asums, by="h") %>%
    select(h, s, k, wh, x, a, asum)
  # bprimedf_base %>% filter(h==1)
  
  # deriv wrt b.s1k1 =
  #    -(wh * x[k=1] * a[s=1] / {asum^2})
  
  # bprime -- how much does delta for an h change wrt a change in any beta
  bprimedf <- bprimedf_base %>%
    mutate(bprime= -(wh * x * a / {asum^2})) %>%
    arrange(h, k, s)
  # bprimedf %>% filter(h==1)
  
  # now get gprime:
  #     gprime(beta)=sum[h] of X * (a * bprime + b * aprime)
  #  bprimedf has most of what we need
  gprime_h <- bprimedf %>%
    select(-a) %>% # drop a as it is also in aprime, below
    left_join(bdf, by="h") %>%
    left_join(aprimedf %>% select(-x), by = c("h", "s", "k")) %>% # drop x as it is in bprime
    mutate(gprime_h= x * (a * bprime + b * aprime))
  # gprime_h %>% filter(h==1)
  
  gprime <- gprime_h %>%
    group_by(s, k) %>%
    summarise(gprime=sum(gprime_h), .groups="drop") %>% # sum over the households h
    arrange(k, s)

  # put it all together to get the gradient by s, k
  # - 2 * (target - g(beta)) * gprime(beta) FOR EACH TARGET AND ADD THEM UP
  # diffs; gprime now we need to cross each gprime with all distances
  grad_base <- expand_grid(s.d=1:s, s.k=1:k, s.gp=1:s, k.gp=1:k) %>%
    left_join(diffsdf %>% select(s, k, diff) %>% rename(s.d=s, s.k=k), by = c("s.d", "s.k")) %>%
    left_join(gprime %>% select(s, k, gprime) %>% rename(s.gp=s, k.gp=k), by = c("s.gp", "k.gp")) %>%
    mutate(grad=-2 * diff * gprime)
  
  graddf <- grad_base %>%
    group_by(s.gp, k.gp) %>%
    summarise(grad=sum(grad), .groups="drop")
  
  graddf$grad
}



# DON'T GO ABOVE HERE ----
idiff <- 1; jbeta <- 1
idiff <- 1; jbeta <- 2
pd3 <- function(idiff, jbeta, ijsk, wh, xmat, beta){
  # determine one element of the Jacobian matrix -- the partial derivative of:
  #   difference in row i of Jacobian (difference between target and calculated target)
  #     with respect to
  #   beta in column j of Jacobian
  
  # ijsk is a 4-column matrix that maps the row number of the difference (idiff) or
  #   column number of the beta (jbeta), which corresponds to the first column of ijsk, named "ij",
  #   to the state for the idiff or jbeta, in column 3 of the matrix, named "s" and to the
  #   characteristic for the idiff or jbeta, in column 4, named "k"
  
  # get the s and k values for the idiff passed to this function
  i.s <- ijsk[idiff, "s"]
  i.k <- ijsk[idiff, "k"]
  
  # get the s and k values for the jbeta passed to this function
  j.s <- ijsk[jbeta, "s"]
  j.k <- ijsk[jbeta, "k"]
  
  # i.s; j.s; i.k; j.k
  
  
  #  Let each difference between a target and its calculated value be:
  #    (target[s, k] - g(beta[s, k]))
  #  
  # For a single target difference, dropping the subscripts for now, but still looking at just one difference
  
  #   g(beta)= weighted sum of X over hh, or sum[h] of X * exp(beta %*% X + delta[h]) where delta[h] is a function of all beta[s,k]
  
  #     Re-express g(beta):
  #          = sum[h] of X * exp(beta*X) * exp(delta[h]) # involving just the beta and x needed for this target
  #      
  #          = sum[h]  of X * A * B where
  #              A=exp(beta * X) i.e., A=a household's weight for a given state before considering delta
  #              B=exp(delta[h]) and delta is a household-specific value and is a function of beta
  
  # we need the NEGATIVE OF THE partial derivative, gprime, wrt a particular beta
  # gprime(beta), still for a single target -- product rule gives:
  #     gprime(beta)=sum[h] of X * (A * Bprime + B * Aprime)
  
  # make a dataframe of households, and calculate each household's A, Aprime, B, and Bprime
  hstub <- tibble(h.df=1:h)
  hsstub <- expand_grid(h.df=1:h, s.df=1:s)
  
  exponents_sum <- function(xrow, sidx){
    # element by element multiplication of an x row by the corresponding beta values
    # this is the sum of the exponents for a given state for a given household
    as.numeric(xrow %*% beta[sidx, ])
  }
  
  Adf <- hstub %>%
    mutate(xh_ki=xmat[h.df, i.k]) %>% # get the x column involved in this target
    rowwise() %>%
    mutate(bx_sum= # sum of the x[, k] values for this person times the beta[i.s, k] coeffs for this target's state
        exponents_sum(xmat[h.df,], i.s), # djb ---- i.s we need this for the state that's in the target
      A=exp(bx_sum),
      # Aprime is the derivative wrt beta-j, which is the x value for that beta
      Aprime=A * xmat[h.df, j.k]) %>%
    ungroup
  
  # for each state, for this person, we need: B_exponent <- beta %*% t(xmat)
  # we also need the sums across states
  # I think? this is the same for all i, j for a given beta so could move out of here and pass it in
  Bsx <- hsstub %>%
    rowwise() %>%
    mutate(bsx=exponents_sum(xmat[h.df, ], s.df), # MAYBE HERE?? ---- s.df here we want beta exponent-sum for each given state
           ebsx=exp(bsx)) %>%
    ungroup()
  
  # make this a matrix -- h x s
  mBsx <- Bsx %>%
    select(h.df, s.df, ebsx) %>%
    pivot_wider(names_from=s.df, values_from=ebsx) %>%
    select(-h.df) %>%
    as.matrix()
  
  mBsx_sum <- rowSums(mBsx)

  # now we are ready to get B and B prime  
  # -(w[h] * x[h, j.k] * es[j.s]  / sum(es)^2)
  
  Bdf <- hstub %>%
    mutate(wh=wh[h.df],
           delta=delta[h.df],
           B=exp(delta),
           xh_kj=xmat[h.df, j.k],
           state_betax =mBsx[h.df, j.s],
           # state_betax =mBsx[h.df, i.s],
           betax_sum=mBsx_sum[h.df],
           Bprime= -wh * xh_kj * state_betax / betax_sum^2)
  
  # gprimeh <- Adf %>% 
  #   left_join(Bdf, by = "h.df") %>%
  #   mutate(gprimeh=- xh_ki * (A * Bprime + B * Aprime))
  
  gprimeh <- Adf %>% 
    left_join(Bdf, by = "h.df") %>%
    # djb fudge ----
    mutate(ABprime = A * Bprime,
           BAprime = B * Aprime,
           xABprime=xh_ki * ABprime) %>%
    # end fudge ----
    mutate(gprimeh=- xh_ki * (ABprime + BAprime))
  
  # gprimeh <- Adf %>% 
  #   left_join(Bdf, by = "h.df") %>%
  #   # djb fudge ----
  #   mutate(ABprime = A * Bprime,
  #          BAprime = B * Aprime,
  #          BAprime=ifelse(idiff != jbeta, 0, BAprime)) %>% # djb ????? ----
  #   mutate(gprimeh=- xh_ki * (ABprime + BAprime))
    # - xmat[, i.k] * (A * Bprime + B * Aprime)
  # gprimeh
  element <- sum(gprimeh$gprimeh)
  element2 <- -sum(gprimeh$xABprime)
  print(gprimeh)
  pdval <- ifelse(idiff==jbeta, element, element2)
  pdval

  # calulate A and Aprime
  # A <- exp(b.s1k1*xhk1 + b.s1k2*xhk2)
  # first get the exponent and then compute the exponentiation
  # we want it for the state of the target in question i.e., for which we have the difference -- i.s

  
  # we get the x values that correspond to the characteristic for the beta, j.k, because
  # we are differentiating wrt that beta
  
  # calculate B and Bprime
  # B is exp(delta), where delta=ln(wh / sum[s]exp(beta[s]X))
  # B has 1 element per household
  #   this simplifies to:
  #   B <- wh / sum[s]exp(beta[s]x)
  # we need each exp(beta[s]X) -- h households, s states
  
  # create B_exponent: a matrix with 1 row per state and 1 column per household
  # for each column it has the exponent for a given state in the expression above
  # that is, it has beta for that state times the k characteristics for the household
  # now gprime
}



# djb test ----
#.. make problem ----
p <- make_problem(h=2, k=1, s=2) # good but for sign
p <- make_problem(h=2, k=1, s=3) # good for 1, 1 but bad for 1, 2 -- pd3(1,2)=pd3(1,1) but should differ

p <- make_problem(h=2, k=2, s=2) # good
p <- make_problem(h=2, k=2, s=3) # not good djb ----
# work on the above -- when we add state #3 it stops working - why? 
# the gprime does not change moving from i1, j1 (correct) to i1, j2 (not correct) - why?

p <- make_problem(h=3, k=1, s=2)
p <- make_problem(h=3, k=2, s=2) # good
p <- make_problem(h=3, k=2, s=3) # not good djb ----
p <- make_problem(h=6, k=5, s=4) # not good djb ----

p <- make_problem(h=5, k=2, s=2) # good through here djb ----
p <- make_problem(h=5, k=2, s=3) # bad results maybe bad data in the function??

p <- make_problem(h=10, k=4, s=8)
p <- make_problem(h=100, k=6, s=20)

p <- make_problem(h=50, k=2, s=2)
p <- make_problem(h=10, k=5, s=2)

p <- make_problem(h=10, k=3, s=4)

p <- make_problem(h=1000, k=3, s=4)
p <- make_problem(h=1000, k=30, s=50)

#.. define indexes ----
ijsk <- expand_grid(s=1:p$s, k=1:p$k) %>% 
  arrange(k, s) %>%
  mutate(ij=row_number()) %>%
  select(ij, s, k) %>%
  as.matrix

#.. extract variables ----
h <- nrow(p$xmat)
s <- nrow(p$targets)
k <- ncol(p$xmat)
h; s; k

targets <- p$targets
xmat <- p$xmat
wh <- p$wh
# whs <- p$whs


#.. define betavec one way or another ----
# betavec <- rep(0, p$s * p$k)
set.seed(2345); betavec <- runif(p$s * p$k)

#.. get beta-dependent values ----
beta <- vtom(betavec, p$s)
delta <- get_delta(wh, beta, xmat)
whs <- get_weights(beta, delta, xmat)

#.. adjust targets if desired ----
etargets <- t(whs) %*% xmat
targets <- etargets
row <- 1; col <- 1
targets[row, col] <- etargets[row, col] + 1

targets; etargets
diff_vec(betavec, wh, xmat, targets)

#.. jacobian finite differences ----
jacobian(diff_vec, x=betavec, wh=p$wh, xmat=p$xmat, targets=targets)
pd3(1, 1, ijsk, p$wh, p$xmat, beta=beta)
pd3(14, 28, ijsk, p$wh, p$xmat, beta=beta)
pd3(30, 31, ijsk, p$wh, p$xmat, beta=beta)


#.. run pd ----
pd3(1, 1, ijsk, p$wh, p$xmat, beta=beta)
pd3(1, 2, ijsk, p$wh, p$xmat, beta=beta)
pd3(1, 3, ijsk, p$wh, p$xmat, beta=beta)
pd3(1, 4, ijsk, p$wh, p$xmat, beta=beta)

pd3(2, 1, ijsk, p$wh, p$xmat, beta=beta)
pd3(3, 1, ijsk, p$wh, p$xmat, beta=beta)
pd3(2, 3, ijsk, p$wh, p$xmat, beta=beta)
pd3(3, 2, ijsk, p$wh, p$xmat, beta=beta)
pd3(4, 1, ijsk, p$wh, p$xmat, beta=beta)

pd3(2, 2, ijsk, p$wh, p$xmat, beta=beta)
pd3(3, 3, ijsk, p$wh, p$xmat, beta=beta)



pd <- function(idiff, jbeta, ijsk, wh, xmat, beta){
  # determine one element of the Jacobian matrix -- the partial derivative of:
  #   difference in row i of Jacobian (difference between target and calculated target)
  #     with respect to
  #   beta in column j of Jacobian
  
  # ijsk is a 4-column matrix that maps the row number of the difference (idiff) or
  #   column number of the beta (jbeta), which corresponds to the first column of ijsk, named "ij",
  #   to the state for the idiff or jbeta, in column 3 of the matrix, named "s" and to the
  #   characteristic for the idiff or jbeta, in column 4, named "k"
  
  # get the s and k values for the idiff passed to this function
  i.s <- ijsk[idiff, "s"]
  i.k <- ijsk[idiff, "k"]
  
  # get the s and k values for the jbeta passed to this function
  j.s <- ijsk[jbeta, "s"]
  j.k <- ijsk[jbeta, "k"]
  
  # i.s; j.s; i.k; j.k
  
  
  #  Let each difference between a target and its calculated value be:
  #    (target[s, k] - g(beta[s, k]))
  #  
  # For a single target difference, dropping the subscripts for now, but still looking at just one difference
  
  #   g(beta)= weighted sum of X over hh, or sum[h] of X * exp(beta %*% X + delta[h]) where delta[h] is a function of all beta[s,k]
  
  #     Re-express g(beta):
  #          = sum[h] of X * exp(beta*X) * exp(delta[h]) # involving just the beta and x needed for this target
  #      
  #          = sum[h]  of X * A * B where
  #              A=exp(beta * X) i.e., A=a household's weight for a given state before considering delta
  #              B=exp(delta[h]) and delta is a household-specific value and is a function of beta
  
  # we need the NEGATIVE OF THE partial derivative, gprime, wrt a particular beta
  # gprime(beta), still for a single target -- product rule gives:
  #     gprime(beta)=sum[h] of X * (A * Bprime + B * Aprime)
  
  # make a dataframe of households, and calculate each household's A, Aprime, B, and Bprime
  hstub <- tibble(h.df=1:h)
  hsstub <- expand_grid(h.df=1:h, s.df=1:s)
  
  exponents_sum <- function(xrow, sidx){
    # element by element multiplication of an x row by the corresponding beta values
    # this is the sum of the exponents for a given state for a given household
    as.numeric(xrow %*% beta[sidx, ])
  }
  
  Adf <- hstub %>%
    mutate(xh_ki=xmat[h.df, i.k]) %>% # get the x column involved in this target
    rowwise() %>%
    mutate(bx_sum= # sum of the x[, k] values for this person times the beta[i.s, k] coeffs for this target's state
             exponents_sum(xmat[h.df,], i.s), # djb ---- i.s we need this for the state that's in the target
           A=exp(bx_sum),
           # Aprime is the derivative wrt beta-j, which is the x value for that beta
           Aprime=A * xmat[h.df, j.k]) %>%
    ungroup
  
  # for each state, for this person, we need: B_exponent <- beta %*% t(xmat)
  # we also need the sums across states
  # I think? this is the same for all i, j for a given beta so could move out of here and pass it in
  Bsx <- hsstub %>%
    rowwise() %>%
    mutate(bsx=exponents_sum(xmat[h.df, ], s.df), # MAYBE HERE?? ---- s.df here we want beta exponent-sum for each given state
           ebsx=exp(bsx)) %>%
    ungroup()
  
  # make this a matrix -- h x s
  mBsx <- Bsx %>%
    select(h.df, s.df, ebsx) %>%
    pivot_wider(names_from=s.df, values_from=ebsx) %>%
    select(-h.df) %>%
    as.matrix()
  
  mBsx_sum <- rowSums(mBsx)
  
  # now we are ready to get B and B prime  
  # -(w[h] * x[h, j.k] * es[j.s]  / sum(es)^2)
  
  Bdf <- hstub %>%
    mutate(wh=wh[h.df],
           delta=delta[h.df],
           B=exp(delta),
           xh_kj=xmat[h.df, j.k],
           state_betax =mBsx[h.df, j.s],
           # state_betax =mBsx[h.df, i.s],
           betax_sum=mBsx_sum[h.df],
           Bprime= -wh * xh_kj * state_betax / betax_sum^2)
  
  # gprimeh <- Adf %>% 
  #   left_join(Bdf, by = "h.df") %>%
  #   mutate(gprimeh=- xh_ki * (A * Bprime + B * Aprime))
  
  gprimeh <- Adf %>% 
    left_join(Bdf, by = "h.df") %>%
    # djb fudge ----
  mutate(ABprime = A * Bprime,
         BAprime = B * Aprime,
         xABprime=xh_ki * ABprime,
         xBAprime=xh_ki * BAprime) %>%
    # end fudge ----
  mutate(gprimeh=- xh_ki * (ABprime + BAprime))
  
  # gprimeh <- Adf %>% 
  #   left_join(Bdf, by = "h.df") %>%
  #   # djb fudge ----
  #   mutate(ABprime = A * Bprime,
  #          BAprime = B * Aprime,
  #          BAprime=ifelse(idiff != jbeta, 0, BAprime)) %>% # djb ????? ----
  #   mutate(gprimeh=- xh_ki * (ABprime + BAprime))
  # - xmat[, i.k] * (A * Bprime + B * Aprime)
  # gprimeh
  element <- sum(gprimeh$gprimeh)
  element2 <- -sum(gprimeh$xABprime)
  # print(gprimeh)
  pdval <- ifelse((idiff==jbeta) | (i.s==j.s), element, element2)
  # pdval <- ifelse(i.s==j.s, element, element2)
  pdval
  
  # calulate A and Aprime
  # A <- exp(b.s1k1*xhk1 + b.s1k2*xhk2)
  # first get the exponent and then compute the exponentiation
  # we want it for the state of the target in question i.e., for which we have the difference -- i.s
  
  
  # we get the x values that correspond to the characteristic for the beta, j.k, because
  # we are differentiating wrt that beta
  
  # calculate B and Bprime
  # B is exp(delta), where delta=ln(wh / sum[s]exp(beta[s]X))
  # B has 1 element per household
  #   this simplifies to:
  #   B <- wh / sum[s]exp(beta[s]x)
  # we need each exp(beta[s]X) -- h households, s states
  
  # create B_exponent: a matrix with 1 row per state and 1 column per household
  # for each column it has the exponent for a given state in the expression above
  # that is, it has beta for that state times the k characteristics for the household
  # now gprime
}

jac <- function(beta, wh, xmat){
  h_n <- nrow(xmat)
  s_n <- nrow(beta)
  k_n <- ncol(beta)
  ij_n <- s_n * k_n
  
  ijsk <- expand_grid(s=1:s_n, k=1:k_n) %>% 
    arrange(k, s) %>%
    mutate(ij=row_number()) %>%
    select(ij, s, k) %>%
    as.matrix
  
  # f <- function(id, jb, ijsk, xmat){
  #   fk <- ijsk[id, "k"]
  #   print(fk)
  #   ijsk[id, "s"] + xmat[id, fk]
  # }
  # 
  # 
  # 
  # jdf <- expand_grid(idiff=1:ij_n, jbeta=1:ij_n) %>%
  #   rowwise() %>%
  #   mutate(jvalue=pd(idiff, jbeta, ijsk, wh, xmat, beta))
  # 
  # jd2 <-jdf %>% 
  #   pivot_wider(names_from = jbeta, values_from=jvalue) %>%
  #   select(-idiff) %>%
  #   as.matrix()
  
  
  jmat <- matrix(0, nrow=ij_n, ncol=ij_n)
  for(idiff in 1:ij_n){
    print(idiff)
     for(jbeta in 1:idiff){
       jmat[idiff, jbeta] <- pd(idiff, jbeta, ijsk, wh, xmat, beta)
     }
  }
  jmat[upper.tri(jmat)] <- t(jmat)[upper.tri(jmat)]
  jmat
}

t1 <- proc.time()
j1 <- jacobian(diff_vec, x=betavec, wh=p$wh, xmat=p$xmat, targets=targets)
t2 <- proc.time()
t2 - t1

t3 <- proc.time()
j2 <- jac(beta, wh, xmat)
t4 <- proc.time()
t4 - t3

j1; j2
(j1 - j2) %>% round(2)
d <- (j1 - j2)

r <- 1:6
c <- 1:6
d <- (j1 - j2)
j1[r, c]; j2[r, c]
d[r, c] %>% round(2)

j1[5, 1]
j2[5, 1]
d[5, 1]
sum(d)
sum(d^2)

pd(5, 1, ijsk, p$wh, p$xmat, beta=beta)

dim(j1)
d <- (j1 - j2)
bad <- which(abs(d) > 0.1, arr.ind = TRUE) %>%
  as_tibble() %>%
  left_join(as_tibble(ijsk) %>% rename(row=ij, s.row=s, k.row=k)) %>%
  left_join(as_tibble(ijsk) %>% rename(col=ij, s.col=s, k.col=k)) %>% 
  arrange(row, col)
bad

guessbad <- expand_grid(sr=ijsk[, "s"], sc=ijsk[, "s"])

ijsk

# fast jacobian ----

pdfast <- function(idiff, jbeta){
  # return the jacobian element for a full 2 vectors: idiff, jbeta
  idiff^2 + jbeta^3
}


xmat <- p$xmat
wh <- p$wh
beta <- rbeta


jacfast <- function(wh, xmat, beta, parallel=TRUE){
  h_n <- nrow(xmat)
  s_n <- nrow(beta)
  k_n <- ncol(beta)
  ij_n <- s_n * k_n
  
  # make the ijsk matrix that maps ij of the jacobian to s and k
  ijsk <- expand_grid(s=1:s_n, k=1:k_n) %>% 
    arrange(k, s) %>%
    mutate(ij=row_number()) %>%
    select(ij, s, k) %>%
    as.matrix
  
  # create a df stub with just the ij values we need, for the lower triangle and diagonal
  ijstub <- expand_grid(idiff=1:ij_n, jbeta=1:ij_n) %>%
    filter(jbeta <= idiff) %>%
    mutate(i_s=ijsk[idiff, "s"],
           i_k=ijsk[idiff, "k"],
           j_s=ijsk[jbeta, "s"],
           j_k=ijsk[jbeta, "k"],
           # define which partial derivatives will include the full gprime, vs just xABprime
           pdfull=ifelse(idiff==jbeta | (i_s==j_s), TRUE, FALSE))
  
  # compute the partial derivatives element by element, vectorizing wherever possible and reusing data
  # wherever possible; this will give us the lower triangle plus diagonal
  
  # pre-define any data and functions that will be reused across different elements
  # create matrices where rows are households and columns are values of interest
  
  f_hh_state_exponents <- function(xmat, beta){
    # get the sum of the exponents for each state for each household
    # obtained by, for each household, multiplying its k characteristics by the corresponding
    # beta coefficients for each state, and summing those k exponents for each state
    # the result is a matrix with 1 row per household, 1 column per state
    # each cell is the relevant exponent for that state
    xmat %*% t(beta)
  }
  
  f_Amat <- function(xmat, beta){
    # the A matrix - note that it does not vary by j_k so we can reuse it for a given xmat, beta
    exp(f_hh_state_exponents(xmat, beta))
  }
  
  f_Aprime <- function(Amat, xmat, j_k, parallel=TRUE){
    # Aprime is the derivative of A wrt beta-j, which is A multiplied by the x value for that beta
    # j_k is a vector of k indexes
    # it will be different for each j_k element
    # thus, return a list that has a matrix for each element of the vector j_k
    f <- function(k) {  
      Amat * xmat[ , k]
    }
    
    lAprime <- llply(j_k, f, .parallel = parallel)
    lAprime
  }
  
  # ahs <- matrix(1:(6*3), ncol=3)
  # ahk <- matrix(1:(6*2), ncol=2)
  # delvec <- seq(10, by=10, length.out = 6)
  # ahs * ahk[, 2]
  # ahk * delvec
  
  # set up Amat and lAprime in advance (a list of Aprime matrices, 1 per each unique j_k) in advance
  # Amat has a row for each household and a column for each state, with exp(Bx) in each column
  # lAprime is a list of Aprime matrices, one per unique k, each matrix has 1 row per hh, 1 col per state
  Amat <- f_Amat(xmat, beta)
  
  print("getting lAprime...")
  lAprime <- f_Aprime(Amat, xmat, 1:k_n, parallel=parallel) # this can be done in parallel
  llply(lAprime, "*" , Bvec, .parallel=parallel)
  llply(1:k_n, function(k) Amat * xmat[, k], .parallel=parallel)
  
  # now prepare B and lBprime
  # B is simply a vector: exp(delta) and is the same for all idiff, jbeta
  delta <- get_delta(wh, beta, xmat)
  Bvec <- exp(delta)
  
  # we need the sum, for each household, of the state exponents, which are in A
  print("getting Bprime_mat...")
  A_hh_sums <- rowSums(Amat) # vector with 1 element per household
  A_hh_sums_squared <- A_hh_sums^2
  # Bprime: # -(w[h] * x[h, j.k] * es[j.s]  / sum(es)^2)
  # for each j, j.k
  # for each hh, we need its state exponents divided (sum of exponents)^2
  # this can be (?) one matrix ?? h x s -- we will draw the j.s column we need
  Bxs_div_BXsum2 <- Amat / A_hh_sums_squared # maybe risk of roundoff here?
  
  # now we need a list with Bprime for each j_k, j_s combination -- i.e., each jbeta
  f_Bprime <- function(wh, xmat, Bxs_div_BXsum2, jbeta, parallel=TRUE){
    # Bprime is the derivative of B wrt beta-j
    # jbeta is a vector
    # it will be different for each element of jbeta
    # thus, return a list that has a vector for each element of the vector jbeta
    f <- function(j) {
      neg_whxkj <- - wh * xmat[ , ijsk[j, "k"]]  # we need the k column for this jbeta
      Bprime <- neg_whxkj * Bxs_div_BXsum2[, ijsk[j, "s"]] # we need the s column for this jbeta
      Bprime # a vector with 1 element per hh
    }
    
    # lBprime <- llply(jbeta, f, .parallel = TRUE)
    Bprime <- laply(jbeta, f, .parallel = parallel)
    t(Bprime) # return 1 row per hh, 1 column per jbeta
  }
  Bprime_mat <- f_Bprime(wh, xmat, Bxs_div_BXsum2, 1:ij_n, parallel=parallel) # h x ij this can be done in parallel
  
  
  print("getting lABprime ...")
  # ABprime: Amat * Bprime_mat -- we need a list of matrices, 1 per column of Amat
  # Amat is fixed, h x s
  # Bprime is an h x ij matrix
  f_ABprime <- function(Amat, Bprime_mat, j_s, parallel=TRUE){
    # return a list with 1 matrix per state
    f <- function(s){
      Amat[, s] * Bprime_mat
    }
    lABprime <- llply(j_s, f, .parallel=parallel)
    lABprime
  }
  lABprime <- f_ABprime(Amat, Bprime_mat, 1:s_n, parallel=parallel)
  
  
  print("getting lBAprime ...")
  # lBAprime is a list of k_n Aprime matrices, one per k, each matrix has 1 row per hh, 1 col per state
  lBAprime <- llply(lAprime, "*" , Bvec, .parallel=parallel) # multiply each matrix in Aprime by the vector Bvec
  
  print("getting lxABprime ...") # djb RETURN this fails on memory write to tempfile?? ----
  # for the non-full elements we need xABprime =xh_ki * ABprime
  # for the full pd elements we need gprimeh=- xh_ki * (ABprime + BAprime) or -xh_ki *ABprime -xh_ki * BAprime
  f_xABprime <- function(k){
    # we have s matrices in lABprime each of which has ij_n columns and h_n rows
    # For each of these we want k_n matrices where each ABprime is multiplied
    # by an xh_ki
    llply(lABprime, "*", xmat[, k], .parallel=parallel)
    
  }
  # lxABprime has k_n elements, each of which has s matrices, with x[, k] multiplied by
  # the ABprime matrix for a given state; it has h_n rows and ij_n columns
  lxABprime <- llply(1:k_n, f_xABprime, .parallel=parallel)
  # djb don't save this full matrix -- get hsums (colsums) ----
  
  print("getting lxBAprime ...")
  f_xBAprime <- function(k){
    # we have k matrices in lBAprime each of which has s_n columns and h_n rows
    # For each of these we want k_n matrices where each BAprime is multiplied
    # by an xh_ki
    llply(lBAprime, "*", xmat[, k], .parallel=parallel)
  }
  # lxABprime has k_n elements, each of which has k_n matrices, with x[, k] multiplied by
  # the BAprime matrix for a given k; it has h_n rows and k_n columns
  lxBAprime <- llply(1:k_n, f_xBAprime, .parallel = parallel)

  # for the non-full elements we need xABprime =xh_ki * ABprime
  # for the full pd elements we need gprimeh=- xh_ki * (ABprime + BAprime) or -xh_ki *ABprime -xh_ki * BAprime
  # collapse each matrix in each list by summing values across households
  # lxBAprime <- llply()
  
  # ijstub %>% filter(idiff==2, jbeta==1)
  
  f <- function(lxABprime){
    # k_n elements, each with s_n matrices, get colsums of each, which will be a vector
    g <- function(matlist){
      laply(matlist, colSums, .parallel=parallel)
    }
    llply(lxABprime, g, .parallel=parallel)
  }
  lxABprime_sums <- f(lxABprime)
  
  f <- function(lxBAprime){
    # k_n elements, each with s_n matrices, get colsums of each, which will be a vector
    g <- function(matlist){
      laply(matlist, colSums, .parallel=parallel)
    }
    llply(lxBAprime, g, .parallel=parallel)
  }
  lxBAprime_sums <- f(lxBAprime)
  
  getpd2 <- function(idiff, jbeta){
    df <- ijstub[(ijstub$idiff==idiff) & (ijstub$jbeta==jbeta), ]
    pdfull <- df$pdfull[1]
    i_k <- df$i_k[1]
    j_k <- df$j_k[1]
    i_s <- df$i_s[1]
    j_s <- df$j_s[1]
    idiff <- df$idiff[1]
    
    if(!pdfull) {
      return(-lxABprime_sums[[j_k]][j_s, idiff])
    } else {
      return(-lxABprime_sums[[j_k]][j_s, idiff] - lxBAprime_sums[[i_k]][j_k, i_s])
      # return(c(lxABprime_sums[[j_k]][j_s, idiff], lxBAprime_sums[[i_k]][j_k, i_s]))
    } 
  }
  
  getpd3 <- function(idiff, jbeta, i_s, i_k, j_s, j_k, pdfull){
    # df <- ijstub[(ijstub$idiff==idiff) & (ijstub$jbeta==jbeta), ]
    # pdfull <- df$pdfull[1]
    # i_k <- df$i_k[1]
    # j_k <- df$j_k[1]
    # i_s <- df$i_s[1]
    # j_s <- df$j_s[1]
    # idiff <- df$idiff[1]
    
    if(!pdfull) {
      return(-lxABprime_sums[[j_k]][j_s, idiff])
    } else {
      return(-lxABprime_sums[[j_k]][j_s, idiff] - lxBAprime_sums[[i_k]][j_k, i_s])
      # return(c(lxABprime_sums[[j_k]][j_s, idiff], lxBAprime_sums[[i_k]][j_k, i_s]))
    } 
  }
  
  
  getpd <- function(idiff, jbeta){
    df <- ijstub[(ijstub$idiff==idiff) & (ijstub$jbeta==jbeta), ]
    pdfull <- df$pdfull[1]
    i_k <- df$i_k[1]
    j_k <- df$j_k[1]
    i_s <- df$i_s[1]
    j_s <- df$j_s[1]
    idiff <- df$idiff[1]
    # print(df)
    # print(j_k); print(j_s); print(idiff)
    
    if(!pdfull) {
      return(-lxABprime[[j_k]][[j_s]][, idiff])
      } else {
        -lxABprime[[j_k]][[j_s]][, idiff] - lxBAprime[[i_k]][[j_k]][, i_s]
    } 
  }
  # getpd(2, 1, 2, 1, 1, 1, FALSE) %>% sum()
  # ijstub %>% filter(idiff==6, jbeta==2)
  # getpd(5, 1)  %>% sum()
  # getpd(6, 2)  %>% sum()
  # getpd(4, 3)  %>% sum()
  # getpd(3, 2)  %>% sum()

  
  
  # now get the pd elements of the jacobian
  ijpd <- ijstub %>%
    group_by(idiff, jbeta) %>%
    mutate(jvalue=getpd3(idiff, jbeta, i_s, i_k, j_s, j_k, pdfull)) %>%
    ungroup
  
  # convert to a matrix with the lower triangle and diagonal
  jmat <-ijpd %>% 
    select(idiff, jbeta, jvalue) %>%
    pivot_wider(names_from = jbeta, values_from=jvalue) %>%
    select(-idiff) %>%
    as.matrix()
  
  jmat[upper.tri(jmat)] <- t(jmat)[upper.tri(jmat)]
  jmat
  return(jmat)
}

j1 <- jacobian(diff_vec, x=as.vector(beta), wh=p$wh, xmat=p$xmat, targets=targets)
j2 <- jacfast(wh, xmat, beta, parallel=TRUE)


p <- make_problem(h=4, s=3, k=2)
p <- make_problem(h=30, s=5, k=4)
p <- make_problem(h=100, s=20, k=8)
p <- make_problem(h=1000, s=20, k=8)
p <- make_problem(h=2000, s=25, k=10) # findiff 23 secs, serial 42, par8 42
p <- make_problem(h=4000, s=30, k=10) # findiff 62, serial 101, par8 80
p <- make_problem(h=6000, s=50, k=20) # par8 821 secs


set.seed(2345); rbetavec <- runif(p$s * p$k)
rbeta <- vtom(rbetavec, nrow(p$targets))

t1a <- proc.time()
j1 <- jacobian(diff_vec, x=rbetavec, wh=p$wh, xmat=p$xmat, targets=p$targets)
t1b <- proc.time()

t2a <- proc.time()
j2 <- jacfast_mda(rbeta, p$wh, p$xmat, ncpar=8)
t2b <- proc.time()
# closeAllConnections()

t1b - t1a
t2b - t2a

sum((j2 - j1)^2)

(j2 - j1) %>% round(2)


jacfast_mda <- function(beta, wh, xmat, ncpar=NULL){
  # ncpar is the number of cores to use in parallel when using Apply; NULL means serial
  # multi dimensional arrays
  # name the dimensions with h, s, and k to reduce risk of error
  dim_rename <- c(h=unname(nrow(xmat)), k=unname(ncol(xmat)))
  dim(xmat) <- dim_rename
  
  dim_rename <- c(h=unname(nrow(beta)), k=unname(ncol(beta)))
  dim(beta) <- c(s=nrow(beta), k=ncol(beta))
  
  h_n <- nrow(xmat)
  s_n <- nrow(beta)
  k_n <- ncol(beta)
  ij_n <- s_n * k_n
  
  # make the ijsk matrix that maps ij of the jacobian to s and k
  ijsk <- expand_grid(s=1:s_n, k=1:k_n) %>% 
    arrange(k, s) %>%
    mutate(ij=row_number()) %>%
    select(ij, s, k) %>%
    as.matrix
  
  # create a df stub with just the ij values we need, for the lower triangle and diagonal
  ijstub <- expand_grid(idiff=1:ij_n, jbeta=1:ij_n) %>%
    filter(jbeta <= idiff) %>%
    mutate(i_s=ijsk[idiff, "s"],
           i_k=ijsk[idiff, "k"],
           j_s=ijsk[jbeta, "s"],
           j_k=ijsk[jbeta, "k"],
           # define which partial derivatives will include the full gprime, vs just xABprime
           pdfull=ifelse(idiff==jbeta | (i_s==j_s), TRUE, FALSE))
  
  # compute the partial derivatives element by element, vectorizing wherever possible and reusing data
  # wherever possible; this will give us the lower triangle plus diagonal
  
  # pre-define any data and functions that will be reused across different elements
  # create matrices where rows are households and columns are values of interest
  
  
  
  # set up Amat and lAprime in advance (a list of Aprime matrices, 1 per each unique j_k) in advance
  # Amat has a row for each household and a column for each state, with exp(Bx) in each column
  # lAprime is a list of Aprime matrices, one per unique k, each matrix has 1 row per hh, 1 col per state
  Amat <- exp(xmat %*% t(beta))
  dim(Amat) <- c(h=nrow(Amat), s=ncol(Amat)) # continue to name dimensions to avoid confusion
  
  
  print("getting aAprime...")
  # get aAprime -- an array of k matrices, each with dimension h x s, where each matrix
  # is Amat multiplied by a column of xmat
  # multiply the "s" dimension (columns) in Amat by the "k" dimension (columns) in xmat
  aAprime <- Apply(data = list(Amat, xmat), 
                   margins=list("s", "k"), 
                   fun = "*",
                   ncores = ncpar)$output1
  # dim(aAprime)
  
  # now prepare B and lBprime
  # B is simply a vector: exp(delta) and is the same for all idiff, jbeta
  delta <- get_delta(wh, beta, xmat)
  # Bvec <- exp(delta)
  Bmat1 <- matrix(exp(delta), ncol=1)
  dim(Bmat1) <- c(h=nrow(Bmat1)) # name the column to be safe
  
  Bxs_div_BXsum2 <- Amat / rowSums(Amat)^2
  
  f_aBprime <- function(wh, xmat, Bxs_div_BXsum2, .ncpar=ncpar){
    neg_wh_x <- -wh * xmat
    
    # create an array of k matrices, each with h rows and s columns, where each matrix is
    # Bxs_div_BXsum2 multiplied by a k column of neg_wh_x
    # is multiplied by a column of neg_wh_x
    # each element of a matrix is deriv of B wrt beta of for a combination of k, s -- B is exp(delta)
    Apply(data = list(Bxs_div_BXsum2, neg_wh_x),
          target_dims="h",
          fun = "*",
          ncores = .ncpar)$output1
  }
  print("getting aBprime...")
  aBprime <- f_aBprime(wh, xmat, Bxs_div_BXsum2)
  
  print("getting aABprime...")
  # caution: aABprime has matrices transposed compared to lABprime
  # dimensions are s, h, k, s -- i.e, row, col, dim3, dim4
  aABprime <- Apply(data = list(aBprime, Amat), # h x s x k, h x s
                    margins=list(c(1, 3), 1:2),
                    fun ="*",
                    ncores=ncpar)$output1
  # we need to rename the dims of aABprime because right now they are s h k s, which repeats
  dim_rename <- dim(aABprime)
  names(dim_rename) <- c("is", "h", "k", "js")
  dim(aABprime) <- dim_rename
  # note that it is symmetric in is, js
  # is <- 1; js <- 2; ih <- 4; ik <- 1
  # aABprime[is, ih, ik, js]
  # aABprime[js, ih, ik, is]
  # Amat, Bvec, aBprime, aABprime so far
  
  print("getting aBAprime...")
  aBAprime <- Apply(data = list(aAprime, Bmat1),
                    target_dims="h",
                    fun = "*",
                    ncores = ncpar)$output1
  # dim(aBAprime) h x s x k
  
  
  # now axABprime_sums - we want every xmat_k times every BAprime_s, then get colSums
  # dim(aABprime)
  # is  h  k js 
  # 3  4  2  3 
  
  # lxABprime_sums # lxABprime -- sums over each person
  print("getting axABprime_sums...")
  axABprime_sums <- Apply(data = list(xmat, aABprime),
                          target_dims=list(c("k", "h"), c("h", "js")),
                          fun = function(m1, m2) m1 %*% m2,
                          output_dims = c("ik", "js"), # names of dimensions returned from function
                          ncores = ncpar)$output1
  # dim(axABprime_sums)
  
  
  print("getting axBAprime_sums...")
  # we have k matrices in aBAprime each of which has h_n rows and s_n columns
  # For each of these we want k_n matrices where each BAprime is multiplied
  # by an xh_ki. Then we want all of the colSums.
  # aBAprime # h, s, k -- multiply each matrix by each k
  # dim(aBAprime)
  # h s k 
  # 4 3 2 
  axBAprime_sums <- Apply(data = list(xmat, aBAprime),
                          target_dims=list(c("k", "h"), c("h", "s")),
                          fun = function(m1, m2) m1 %*% m2,
                          #fun = f,
                          output_dims = c("ik", "s"), # name the dimensions returned from function
                          ncores = ncpar)$output1
  # dim(axBAprime_sums)
  
  # now we have everything we need to create the jacobian ----
  print("building jacobian matrix...")
  # ijsk
  # ijstub
  ijpd <- ijstub %>%
    group_by(idiff, jbeta) %>%
    mutate(jvalue=if_else(!pdfull, 
                          -axABprime_sums[j_k, j_s, i_s, i_k],
                          -axABprime_sums[j_k, j_s, i_s, i_k] - axBAprime_sums[i_k, i_s, j_k])) %>%
    ungroup
  # ijpd
  
  # convert to a matrix with the lower triangle and diagonal
  jmat <-ijpd %>% 
    select(idiff, jbeta, jvalue) %>%
    pivot_wider(names_from = jbeta, values_from=jvalue) %>%
    select(-idiff) %>%
    as.matrix()
  
  jmat[upper.tri(jmat)] <- t(jmat)[upper.tri(jmat)]
  return(jmat)
}







library(doParallel)
library(foreach)
cores <- detectCores()
cores
# registerDoParallel(cores=cores)
cl = makeCluster(cores)
registerDoParallel(cl)

t2a <- proc.time()
j2 <- jacfast(p$wh, p$xmat, rbeta, parallel=FALSE)
t2b <- proc.time()

stopCluster(cl)

t1b - t1a
t2b - t2a

sum((j2 - j1)^2)
(j2 - j1) %>% round(2)

getpd2(1, 1)

getpd2(2, 1)
getpd2(2, 2)

getpd2(3, 1)
getpd2(3, 2)
getpd2(3, 3)

getpd2(4, 1)
getpd2(4, 2)
getpd2(4, 3)

getpd2(5, 1)
getpd2(5, 2)
getpd2(5, 3)
getpd2(5, 4)

getpd2(6, 1)
getpd2(6, 2)
getpd2(6, 3)
getpd2(6, 4)
getpd2(6, 5)


# for each state, for this person, we need: B_exponent <- beta %*% t(xmat)
# we also need the sums across states
# I think? this is the same for all i, j for a given beta so could move out of here and pass it in
Bsx <- hsstub %>%
  rowwise() %>%
  mutate(bsx=exponents_sum(xmat[h.df, ], s.df), # MAYBE HERE?? ---- s.df here we want beta exponent-sum for each given state
         ebsx=exp(bsx)) %>%
  ungroup()

# make this a matrix -- h x s
mBsx <- Bsx %>%
  select(h.df, s.df, ebsx) %>%
  pivot_wider(names_from=s.df, values_from=ebsx) %>%
  select(-h.df) %>%
  as.matrix()

mBsx_sum <- rowSums(mBsx)

# now we are ready to get B and B prime  
# -(w[h] * x[h, j.k] * es[j.s]  / sum(es)^2)

Bdf <- hstub %>%
  mutate(wh=wh[h.df],
         delta=delta[h.df],
         B=exp(delta),
         xh_kj=xmat[h.df, j.k],
         state_betax =mBsx[h.df, j.s],
         # state_betax =mBsx[h.df, i.s],
         betax_sum=mBsx_sum[h.df],
         Bprime= -wh * xh_kj * state_betax / betax_sum^2)



# djb test ----
p <- make_problem(h=4, s=3, k=2)


#.. make problem ----
p <- make_problem(h=2, k=1, s=2)
p <- make_problem(h=2, k=1, s=3)

p <- make_problem(h=2, k=2, s=2)
p <- make_problem(h=2, k=2, s=3)

p <- make_problem(h=3, k=1, s=2)
p <- make_problem(h=3, k=2, s=2) 
p <- make_problem(h=3, k=2, s=3) 
p <- make_problem(h=6, k=5, s=4) # not good djb ----

p <- make_problem(h=5, k=2, s=2) # good through here djb ----
p <- make_problem(h=5, k=2, s=3) # bad results maybe bad data in the function??

p <- make_problem(h=10, k=4, s=8)
p <- make_problem(h=100, k=6, s=20)

p <- make_problem(h=50, k=2, s=2)
p <- make_problem(h=10, k=5, s=2)

p <- make_problem(h=10, k=3, s=4)

p <- make_problem(h=1000, k=3, s=4)
p <- make_problem(h=1000, k=30, s=50)

#.. define indexes ----
#.. extract variables ----
targets <- p$targets
xmat <- p$xmat
wh <- p$wh
# whs <- p$whs


#.. define betavec one way or another ----
# betavec <- rep(0, p$s * p$k)
set.seed(2345); betavec <- runif(p$s * p$k)

#.. get beta-dependent values ----
beta <- vtom(betavec, p$s)
delta <- get_delta(wh, beta, xmat)
whs <- get_weights(beta, delta, xmat)

#.. adjust targets if desired ----
etargets <- t(whs) %*% xmat
targets <- etargets
row <- 1; col <- 1
targets[row, col] <- etargets[row, col] + 1

targets; etargets
diff_vec(betavec, wh, xmat, targets)





jactest <- function(a, b, c){ 
  f <- function(i, j){
    # ijjac$jvalue[match(i, ijjac$idiff), match(j, ijjac$jbeta)]
    i <- 1; j <- 1
    # ijjac$jvalue[match(i, ijjac$idiff) & match(j, ijjac$jbeta)]
  }
  jac <- matrix(nrow=ij_n, ncol=ij_n)
  
  # if(!pdfull) {
  #   return(-lxABprime_sums[[j_k]][j_s, idiff])
  # } else {
  #   return(-lxABprime_sums[[j_k]][j_s, idiff] - lxBAprime_sums[[i_k]][j_k, i_s])
  outer(X, Y, FUN = "*", ...)
  
  outer(1:ij_n, 1:ij_n, FUN = f)
  
  # >   jmat
  # 1          2          3          4          5          6
  # [1,] -8.741518   3.784066   4.957452  -8.043498   3.543698   4.499800
  # [2,]  3.784066 -10.737036   6.952971   3.543698  -9.820095   6.276398
  # [3,]  4.957452   6.952971 -11.910423   4.499800   6.276398 -10.776198
  # [4,] -8.043498   3.543698   4.499800 -10.453183   4.865392   5.587791
  # [5,]  3.543698  -9.820095   6.276398   4.865392 -12.479800   7.614408
  # [6,]  4.499800   6.276398 -10.776198   5.587791   7.614408 -13.202198
  
  # jacobian(diff_vec, x=as.vector(beta), wh=p$wh, xmat=p$xmat, targets=p$targets) gives same
  
  getpd <- function(idiff, jbeta){
    df <- ijstub[(ijstub$idiff==idiff) & (ijstub$jbeta==jbeta), ]
    pdfull <- df$pdfull[1]
    i_k <- df$i_k[1]
    j_k <- df$j_k[1]
    i_s <- df$i_s[1]
    j_s <- df$j_s[1]
    idiff <- df$idiff[1]
    # print(df)
    # print(j_k); print(j_s); print(idiff)
    
    if(!pdfull) {
      return(-lxABprime[[j_k]][[j_s]][, idiff])
    } else {
      -lxABprime[[j_k]][[j_s]][, idiff] - lxBAprime[[i_k]][[j_k]][, i_s]
    } 
  }
  
  
  # list-based work below here ----
  f_hh_state_exponents <- function(xmat, beta){
    # get the sum of the exponents for each state for each household
    # obtained by, for each household, multiplying its k characteristics by the corresponding
    # beta coefficients for each state, and summing those k exponents for each state
    # the result is a matrix with 1 row per household, 1 column per state
    # each cell is the relevant exponent for that state
    xmat %*% t(beta)
  }
  
  f_Amat <- function(xmat, beta){
    # the A matrix - note that it does not vary by j_k so we can reuse it for a given xmat, beta
    exp(f_hh_state_exponents(xmat, beta))
  }
  
  f_Aprime <- function(Amat, xmat, j_k, parallel=TRUE){
    # Aprime is the derivative of A wrt beta-j, which is A multiplied by the x value for that beta
    # j_k is a vector of k indexes
    # it will be different for each j_k element
    # thus, return a list that has a matrix for each element of the vector j_k
    f <- function(k) {  
      Amat * xmat[ , k]
    }
    
    lAprime <- llply(j_k, f, .parallel = parallel)
    lAprime
  }
  
  
  print("getting lxBAprime ...")
  f_xBAprime <- function(k){
    # we have k matrices in lBAprime each of which has s_n columns and h_n rows
    # For each of these we want k_n matrices where each BAprime is multiplied
    # by an xh_ki
    llply(lBAprime, "*", xmat[, k], .parallel=parallel)
  }
  # lxABprime has k_n elements, each of which has k_n matrices, with x[, k] multiplied by
  # the BAprime matrix for a given k; it has h_n rows and k_n columns
  lxBAprime <- llply(1:k_n, f_xBAprime, .parallel = parallel)
  
  # for the non-full elements we need xABprime =xh_ki * ABprime
  # for the full pd elements we need gprimeh=- xh_ki * (ABprime + BAprime) or -xh_ki *ABprime -xh_ki * BAprime
  # collapse each matrix in each list by summing values across households
  # lxBAprime <- llply()
  
  # ijstub %>% filter(idiff==2, jbeta==1)
  
  f <- function(lxABprime){
    # k_n elements, each with s_n matrices, get colsums of each, which will be a vector
    g <- function(matlist){
      laply(matlist, colSums, .parallel=parallel)
    }
    llply(lxABprime, g, .parallel=parallel)
  }
  lxABprime_sums <- f(lxABprime)
  
  f <- function(lxBAprime){
    # k_n elements, each with s_n matrices, get colsums of each, which will be a vector
    g <- function(matlist){
      laply(matlist, colSums, .parallel=parallel)
    }
    llply(lxBAprime, g, .parallel=parallel)
  }
  lxBAprime_sums <- f(lxBAprime)
  
  
  f_xABprime <- function(k){
    # we have s matrices in lABprime each of which has ij_n columns and h_n rows
    # For each of these we want k_n matrices where each ABprime is multiplied
    # by an xh_ki
    llply(lABprime, "*", xmat[, k], .parallel=parallel)
    
  }
  # lxABprime has k_n elements, each of which has s matrices, with x[, k] multiplied by
  # the ABprime matrix for a given state; it has h_n rows and ij_n columns
  lxABprime <- llply(1:k_n, f_xABprime, .parallel=parallel)
  
  
  
  
  # we need the sum, for each household, of the state exponents, which are in A
  print("getting Bprime_mat...")
  A_hh_sums <- rowSums(Amat) # vector with 1 element per household
  A_hh_sums_squared <- rowSums(Amat)^2 # vector with 1 element per household
  
  # Bprime: # -(w[h] * x[h, j.k] * es[j.s]  / sum(es)^2)
  # for each j, j.k
  # for each hh, we need its state exponents divided (sum of exponents)^2
  # this can be (?) one matrix ?? h x s -- we will draw the j.s column we need
  Bxs_div_BXsum2 <- Amat / rowSums(Amat)^2 # maybe risk of roundoff here?
  # dim(Bxs_div_BXsum2)
  
  # now we need a list with Bprime for each j_k, j_s combination -- i.e., each jbeta
  f_Bprime <- function(wh, xmat, Bxs_div_BXsum2, jbeta, parallel=TRUE){
    # Bprime is the derivative of B wrt beta-j
    # jbeta is a vector
    # it will be different for each element of jbeta
    # thus, return a list that has a vector for each element of the vector jbeta
    f <- function(j) {
      neg_whxkj <- - wh * xmat[ , ijsk[j, "k"]]  # we need the k column for this jbeta
      Bprime <- neg_whxkj * Bxs_div_BXsum2[, ijsk[j, "s"]] # we need the s column for this jbeta
      Bprime # a vector with 1 element per hh
    }
    
    # lBprime <- llply(jbeta, f, .parallel = TRUE)
    Bprime <- laply(jbeta, f, .parallel = parallel)
    t(Bprime) # return 1 row per hh, 1 column per jbeta
  }
  Bprime_mat <- f_Bprime(wh, xmat, Bxs_div_BXsum2, 1:ij_n, parallel=parallel) # h x ij this can be done in parallel
  
  f_aBprime <- function(wh, xmat, Bxs_div_BXsum2, .ncpar=ncpar){
    neg_wh_x <- -wh * xmat
    
    # create an array of k matrices, each with h rows and s columns, where each matrix is
    # Bxs_div_BXsum2 multiplied by a k column of neg_wh_x
    # is multiplied by a column of neg_wh_x
    # each element of a matrix is deriv of B wrt beta of for a combination of k, s -- B is exp(delta)
    Apply(data = list(Bxs_div_BXsum2, neg_wh_x),
          target_dims="h",
          # margins=2,
          fun = "*",
          ncores = .ncpar)$output1
  }

  print("getting lABprime ...")
  # ABprime: Amat * Bprime_mat -- we need a list of matrices, 1 per column of Amat
  # Amat is fixed, h x s
  # Bprime is an h x ij matrix
  f_ABprime <- function(Amat, Bprime_mat, j_s, parallel=TRUE){
    # return a list with 1 matrix per state
    f <- function(s){
      Amat[, s] * Bprime_mat
    }
    lABprime <- llply(j_s, f, .parallel=parallel)
    lABprime
  }
  lABprime <- f_ABprime(Amat, Bprime_mat, 1:s_n, parallel=parallel)
  
  # caution: aABprime has matrices transposed compared to lABprime
  
  
  print("getting lBAprime ...")
  # lBAprime is a list of k_n Aprime matrices, one per k, each matrix has 1 row per hh, 1 col per state
  lBAprime <- llply(lAprime, "*" , Bvec, .parallel=parallel) # multiply each matrix in Aprime by the vector Bvec
  
  print("getting lxABprime ...") # djb RETURN this fails on memory write to tempfile?? ----
  # for the non-full elements we need xABprime =xh_ki * ABprime
  # for the full pd elements we need gprimeh=- xh_ki * (ABprime + BAprime) or -xh_ki *ABprime -xh_ki * BAprime
  f_xABprime <- function(k){
    # we have s matrices in lABprime each of which has ij_n columns and h_n rows
    # For each of these we want k_n matrices where each ABprime is multiplied
    # by an xh_ki
    llply(lABprime, "*", xmat[, k], .parallel=parallel)
    
  }
  # lxABprime has k_n elements, each of which has s matrices, with x[, k] multiplied by
  # the ABprime matrix for a given state; it has h_n rows and ij_n columns
  lxABprime <- llply(1:k_n, f_xABprime, .parallel=parallel)
  # djb don't save this full matrix -- get hsums (colsums) ----
  
  print("getting lxBAprime ...")
  f_xBAprime <- function(k){
    # we have k matrices in lBAprime each of which has s_n columns and h_n rows
    # For each of these we want k_n matrices where each BAprime is multiplied
    # by an xh_ki
    llply(lBAprime, "*", xmat[, k], .parallel=parallel)
  }
  # lxABprime has k_n elements, each of which has k_n matrices, with x[, k] multiplied by
  # the BAprime matrix for a given k; it has h_n rows and k_n columns
  lxBAprime <- llply(1:k_n, f_xBAprime, .parallel = parallel)
  
  # for the non-full elements we need xABprime =xh_ki * ABprime
  # for the full pd elements we need gprimeh=- xh_ki * (ABprime + BAprime) or -xh_ki *ABprime -xh_ki * BAprime
  # collapse each matrix in each list by summing values across households
  # lxBAprime <- llply()
  
  # ijstub %>% filter(idiff==2, jbeta==1)
  
  f <- function(lxABprime){
    # k_n elements, each with s_n matrices, get colsums of each, which will be a vector
    g <- function(matlist){
      laply(matlist, colSums, .parallel=parallel)
    }
    llply(lxABprime, g, .parallel=parallel)
  }
  lxABprime_sums <- f(lxABprime)
  
  f <- function(lxBAprime){
    # k_n elements, each with s_n matrices, get colsums of each, which will be a vector
    g <- function(matlist){
      laply(matlist, colSums, .parallel=parallel)
    }
    llply(lxBAprime, g, .parallel=parallel)
  }
  lxBAprime_sums <- f(lxBAprime)
  
  getpd2 <- function(idiff, jbeta){
    df <- ijstub[(ijstub$idiff==idiff) & (ijstub$jbeta==jbeta), ]
    pdfull <- df$pdfull[1]
    i_k <- df$i_k[1]
    j_k <- df$j_k[1]
    i_s <- df$i_s[1]
    j_s <- df$j_s[1]
    idiff <- df$idiff[1]
    
    if(!pdfull) {
      return(-lxABprime_sums[[j_k]][j_s, idiff])
    } else {
      return(-lxABprime_sums[[j_k]][j_s, idiff] - lxBAprime_sums[[i_k]][j_k, i_s])
      # return(c(lxABprime_sums[[j_k]][j_s, idiff], lxBAprime_sums[[i_k]][j_k, i_s]))
    } 
  }
  
  getpd3 <- function(idiff, jbeta, i_s, i_k, j_s, j_k, pdfull){
    # df <- ijstub[(ijstub$idiff==idiff) & (ijstub$jbeta==jbeta), ]
    # pdfull <- df$pdfull[1]
    # i_k <- df$i_k[1]
    # j_k <- df$j_k[1]
    # i_s <- df$i_s[1]
    # j_s <- df$j_s[1]
    # idiff <- df$idiff[1]
    
    if(!pdfull) {
      return(-lxABprime_sums[[j_k]][j_s, idiff])
    } else {
      return(-lxABprime_sums[[j_k]][j_s, idiff] - lxBAprime_sums[[i_k]][j_k, i_s])
      # return(c(lxABprime_sums[[j_k]][j_s, idiff], lxBAprime_sums[[i_k]][j_k, i_s]))
    } 
  }
  
  
  getpd <- function(idiff, jbeta){
    df <- ijstub[(ijstub$idiff==idiff) & (ijstub$jbeta==jbeta), ]
    pdfull <- df$pdfull[1]
    i_k <- df$i_k[1]
    j_k <- df$j_k[1]
    i_s <- df$i_s[1]
    j_s <- df$j_s[1]
    idiff <- df$idiff[1]
    # print(df)
    # print(j_k); print(j_s); print(idiff)
    
    if(!pdfull) {
      return(-lxABprime[[j_k]][[j_s]][, idiff])
    } else {
      -lxABprime[[j_k]][[j_s]][, idiff] - lxBAprime[[i_k]][[j_k]][, i_s]
    } 
  }
  # getpd(2, 1, 2, 1, 1, 1, FALSE) %>% sum()
  # ijstub %>% filter(idiff==6, jbeta==2)
  # getpd(5, 1)  %>% sum()
  # getpd(6, 2)  %>% sum()
  # getpd(4, 3)  %>% sum()
  # getpd(3, 2)  %>% sum()
  
  
  
  # now get the pd elements of the jacobian
  ijpd <- ijstub %>%
    group_by(idiff, jbeta) %>%
    mutate(jvalue=getpd3(idiff, jbeta, i_s, i_k, j_s, j_k, pdfull)) %>%
    ungroup
  
  # convert to a matrix with the lower triangle and diagonal
  jmat <-ijpd %>% 
    select(idiff, jbeta, jvalue) %>%
    pivot_wider(names_from = jbeta, values_from=jvalue) %>%
    select(-idiff) %>%
    as.matrix()
  
  jmat[upper.tri(jmat)] <- t(jmat)[upper.tri(jmat)]
  jmat
  return(jmat)
}

