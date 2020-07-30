
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

etargs_vec <- function(betavec, wh, xmat, s){
  # return a vector of calculated targets and corresponding
  # values calculated given a beta vector, household weights, and x matrix
  beta <- vtom(betavec, nrows=s)
  delta <- get_delta(wh, beta, xmat)
  whs <- get_weights(beta, delta, xmat)
  etargets <- t(whs) %*% xmat
  as.vector(etargets)
}


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


p <- make_problem(h=2, k=1, s=2)
p <- make_problem(h=4, k=2, s=3) # use this
p <- make_problem(h=20, k=4, s=8)

sval <- rep(0, p$s * p$k)
# sval <- runif(p$s*p$k)
betavec <- rep(0, p$s * p$k)
targets <- p$targets
xmat <- p$xmat
wh <- p$wh
beta <- vtom(betavec, p$s)
delta <- get_delta(wh, beta, xmat)
whs <- get_weights(beta, delta, xmat)

# make targets that are all hit, except 1
etargets <- t(whs) %*% xmat
targets <- etargets
row <- 1; col <- 1
targets[row, col] <- etargets[row, col] + 1
diff_vec(betavec, wh, xmat, targets)

grad(sse_fn, x=sval, wh=p$wh, xmat=p$xmat, targets=targets) %>% round(4)
grad_sse_v2(sval, wh=p$wh, xmat=p$xmat, targets=targets) %>% round(4)

jacobian(diff_vec, x=sval, wh=p$wh, xmat=p$xmat, targets=targets)
# gprime$gprime
# the diagonal of the jacobian is right...except for the sign



jac <- function(betavec, wh, xmat, targets){
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
  
  # we need functions that get the s or k associated with a given index
  get_s <- function(idx) {skstub$s[idx]}
  get_k <- function(idx) {skstub$k[idx]}
  
  xvec <- as.vector(xmat) # this is in the proper order 
  xdf <- hkstub %>% mutate(x=as.vector(xmat))
  xmat
  diffs <- diff_vec(betavec, wh, xmat, targets)
  
  # create a long form of the jacobian with i indexing rows [differences],
  # and j indexing columns [betas]
  # each element is a partial derivative of d[i] wrt beta[j] and will depend on the h's and k's also
  
  irows_diff <- expand_grid(s.i=1:s, k.i=1:k) %>% 
    arrange(k.i, s.i) %>%
    mutate(i=row_number(), diff=diffs, beta.diff=betavec) %>%
    select(i, s.i, k.i, diff, beta.diff)
  
  jcols_beta <- expand_grid(s.j=1:s, k.j=1:k) %>% 
    arrange(k.j, s.j) %>%
    mutate(j=row_number(), beta.pd=betavec) %>% # we'll need the beta from this iteration
    select(j, s.j, k.j, beta.pd)
  
  # the long version of the jacobian
  jlong_stub <- expand_grid(i=1:(s * k), j=1:(s * k)) %>%
    left_join(irows_diff, by = "i") %>%
    left_join(jcols_beta, by = "j")
  
  
  get_pd <- function(df){
    # get minus gprime(beta) partial derivative of diff-i (diff) wrt beta-j (beta.pd, for beta used in the partial derivative)
    
    # where
    #   g(beta)= weighted sum of X over hh, or sum[h] of X * exp(beta %*% X + delta[h]) where delta[h] is a function of all beta[s,k]
    
    #     Re-express g(beta):
    #          = sum[h] of X * exp(beta*X) * exp(delta[h]) # involving just the beta and x needed for this target
    #      
    #          = sum[h]  of X * a * b where a=exp(beta *X) and b=exp(delta[h]) and delta is a function of beta
    
    # gprime(beta), still for a single target -- product rule gives:
    #     gprime(beta)=sum[h] of X[k.d] * (a * bprime + b * aprime)
    # get the x values for all households for this target
    
    
    x_k <- tibble(x=xmat[, df$k.i]) # the x values for this target's k, one value per household
    beta_s <- tibble(k=1:k, beta_s.i=beta[df$s.i, ]) # the beta values for this target's state, one value per k
    # df <- df %>%
    #   mutate(xsum=sum(xmat[, get_k(.$i)])) # we use x[, k] for this target thus we get i, not j
    
    # # a = exp(beta * x)  - this is the same for all
    # adf_base <- xdf %>%
    #   left_join(betadf, by="k") %>%
    #   mutate(a_exponent=beta * x) %>% # we will raise e to this power
    #   select(h, s, k, x, beta, a_exponent)
    # adf_base %>% filter(h==1)
    # 
    # adf <- adf_base %>%
    #   group_by(s, h) %>% # weights are specific to state and household
    #   # we sum the k elements of the exponent and then raise e to that power
    #   summarise(a=exp(sum(a_exponent)), .groups="drop") %>% # these are the state weights for each hh BEFORE delta impact
    #   select(h, s, a)
    # adf %>% filter(h==1)
    amat <- xmat %*% t(beta) # each row has the a_exponent values, 1 per state, for a household
    a <- exp(amat)
    asums <- rowSums(a)
    
    # a <- exp(xmat %*% beta[df$s.i, ]) # one value per h, weight[h, s] calculated without delta
    # we also need asum - the sum of a over all states for this person
    
    # asums <- adf %>%
    #   group_by(h) %>%
    #   summarise(asum=sum(a), .groups="drop")
    # 
    # bprimedf_base <- adf %>%
    #   left_join(whdf, by = "h") %>%
    #   left_join(xdf, by="h") %>%
    #   left_join(asums, by="h") %>%
    #   select(h, s, k, wh, x, a, asum)
    # bprimedf_base %>% filter(h==1)
    
    # deriv wrt b.s1k1 =
    #    -(wh * x[k=1] * a[s=1] / {asum^2})
    
    # bprime -- how much does delta for an h change wrt a change in any beta
    # bprimedf <- bprimedf_base %>%
    #   mutate(bprime= -(wh * x * a / {asum^2})) %>%
    #   arrange(h, k, s)
    
    ab <- x_k %>%
      left_join(df, by = character()) %>%
      mutate(wh=wh,
             a=a[, df$s.i], 
             aprime=x * a,
             b=exp(delta),
             asum=asums,
             bprime = -(wh * x * a / {asum^2}),
             gprime_h = x * (a * bprime + b * aprime))
    # beta_s %>%
    #   left_join(df, by = character())
    ab
  }
  
  df <- jlong_stub %>% filter(i==1, j==1)
  get_pd(df)
  
  
  tmp <- jlong_stub %>%
    filter(i %in% c(1, 4), j %in% 1) %>%
    group_by(i, j) %>%
    group_modify(~get_pd(.x))
  tmp
  tmp %>%
    group_by(i, j) %>%
    summarise(gprime=sum(gprime_h), .groups="drop")
  

  diffs <- diff_vec(betavec, wh, xmat, targets)
  etargets <- etargs_vec(betavec, wh, xmat, s)
  
  # start making data frames
  diffsdf <- skstub %>% 
    mutate(target=as.vector(targets),
           etarget=etargets,
           diff=diffs)
  
  whdf <- hstub %>% mutate(wh=wh)
  
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


# djb - come back here! ----
p <- make_problem(h=2, k=1, s=2)
p <- make_problem(h=2, k=2, s=2)
p <- make_problem(h=4, k=2, s=3) # use this
p <- make_problem(h=20, k=4, s=8)

sval <- rep(0, p$s * p$k)
# sval <- runif(p$s*p$k)
betavec <- rep(0, p$s * p$k)
targets <- p$targets
xmat <- p$xmat
wh <- p$wh
beta <- vtom(betavec, p$s)
delta <- get_delta(wh, beta, xmat)
whs <- get_weights(beta, delta, xmat)

# make targets that are all hit, except 1
etargets <- t(whs) %*% xmat
targets <- etargets
row <- 1; col <- 1
targets[row, col] <- etargets[row, col] + 1
diff_vec(betavec, wh, xmat, targets)

h <- nrow(xmat)
k <- ncol(xmat)
s <- nrow(targets)

# xvec <- as.vector(xmat) # this is in the proper order 
# xdf <- hkstub %>% mutate(x=as.vector(xmat))
xmat
diffs <- diff_vec(betavec, wh, xmat, targets)

irows_diff <- expand_grid(s.i=1:s, k.i=1:k) %>% 
  arrange(k.i, s.i) %>%
  mutate(i=row_number(), diff=diffs, beta.diff=betavec) %>%
  select(i, s.i, k.i, diff, beta.diff)

jcols_beta <- expand_grid(s.j=1:s, k.j=1:k) %>% 
  arrange(k.j, s.j) %>%
  mutate(j=row_number(), beta.pd=betavec) %>% # we'll need the beta from this iteration
  select(j, s.j, k.j, beta.pd)

p
# let's get J.d1b1
# now J.d1b3
# now J.d3b3
irows_diff # i=d1  # for d, we're still in row 1 but now the j col will be different -- diff b
# so idiff has s=1, k=1
i.s <- 1; i.k <- 1
i.s <- 1; i.k <- 2

jcols_beta
# jbeta has s=1, k=1
# so for our 2 people we have
j.s <- 1; j.k <- 2
j.s <- 1; j.k <- 2
xmat
# A <- exp(b.s1k1*xhk1 + b.s1k2*xhk2)
A_exponent <- xmat %*% beta[i.s, ] # (h x k) x (1 x k) # we need the beta for this difference -- its state
A <- exp(A_exponent) # this is their weight without delta
A
Aprime <- xmat[, j.k] * A # element by element, 1 per person we need j.k because it is wrt b.sjkj
Aprime

# now B, which is exp(delta), where delta=ln(wh / sum[s]exp(beta[s]X))
# B <- wh / sum[s]exp(beta[s]x)
# we need each exp(beta[s]X) -- 2 households, 2 states
B_exponent <- beta %*% t(xmat) # (s x k) x (k x h) = s x h
B_exponent # s x h
Bmat <- exp(B_exponent) # s x h
Bmat

B <- wh / colSums(Bmat) # 1 element per hh
B <- t(t(B))
B
Bprime <- - wh * xmat[, j.k] * Bmat[j.s, ] / colSums(Bmat)^2
Bprime

gprimeh <- xmat[, i.k] * (A * Bprime + B * Aprime)
gprimeh
colSums(gprimeh)

jacobian(diff_vec, x=sval, wh=p$wh, xmat=p$xmat, targets=targets)


ijsk <- expand_grid(s=1:s, k=1:k) %>% 
  arrange(k, s) %>%
  mutate(ij=row_number()) %>%
  select(ij, s, k) %>%
  as.matrix

idiff <- 2; jbeta <- 1
idiff <- 1; jbeta <- 1

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
  
  # we need the partial derivative, gprime, wrt a particular beta
  # gprime(beta), still for a single target -- product rule gives:
  #     gprime(beta)=sum[h] of X * (A * Bprime + B * Aprime)
  
  # calulate A and Aprime
  # A <- exp(b.s1k1*xhk1 + b.s1k2*xhk2)
  # first get the exponent and then compute the exponentiation
  # we want it for the state of the target in question i.e., for which we have the difference -- i.s
  A_exponent <- xmat %*% beta[i.s, ] # (h x k) x (k x 1) # we need the beta for this difference -- its state
  # A_exponent <- xmat %*% beta[j.s, ] # test (h x k) x (k x 1) # we need the beta for this difference -- its state
  # A_exponent has 1 row per household, and 1 column, with the exponent
  A <- exp(A_exponent) # this is their weight without delta
  # A # A also has 1 row per household, and 1 column, with the exponent
  
  # we get the x values that correspond to the characteristic for the beta, j.k, because
  # we are differentiating wrt that beta
  Aprime <- xmat[, j.k] * A # element by element multiplication
  # Aprime <- xmat[, i.k] * A # test
  # Aprime # 1 row per household, 1 column
  
  # calculate B and Bprime
  # B is exp(delta), where delta=ln(wh / sum[s]exp(beta[s]X))
  # B has 1 element per household
  #   this simplifies to:
  #   B <- wh / sum[s]exp(beta[s]x)
  # we need each exp(beta[s]X) -- h households, s states
  
  # create B_exponent: a matrix with 1 row per state and 1 column per household
  # for each column it has the exponent for a given state in the expression above
  # that is, it has beta for that state times the k characteristics for the household
  B_exponent <- beta %*% t(xmat) # (s x k) x (k x h) = s x h -- 1 row per state
  # B_exponent # s x h
  Bmat <- exp(B_exponent) # s x h
  # Bmat
  
  B <- wh / colSums(Bmat) # 1 element per hh
  B <- t(t(B))
  # B
  Bprime <- - wh * xmat[, j.k] * Bmat[j.s, ] / colSums(Bmat)^2
  # Bprime <- - wh * xmat[, j.k] * Bmat[i.s, ] / colSums(Bmat)^2 # test
  # Bprime
  
  # now gprime
  gprimeh <- - xmat[, i.k] * (A * Bprime + B * Aprime)
  # gprimeh
  element <- colSums(gprimeh)
  element
}



jacobian(diff_vec, x=bvec, wh=p$wh, xmat=p$xmat, targets=p$targets)

i <- 1; j <- 1
i <- 1; j <- 2
i <- 2; j <- 1
i <- 3; j <- 1
i <- 3; j <- 2
i <- 3; j <- 4
i <- 2; j <- 2

bvec <- rep(0, p$s * p$k)
idiff <- 2; jbeta <- 1; wh <- p$wh; xmat <- p$xmat; beta <- vtom(bvec, nrows=nrow(p$targets))
pd(i, j, ijsk, p$wh, p$xmat, beta=vtom(bvec, nrows=nrow(p$targets)))
pd(1, 1, ijsk, p$wh, p$xmat, beta=vtom(bvec, nrows=nrow(p$targets)))
pd(2, 1, ijsk, p$wh, p$xmat, beta=vtom(bvec, nrows=nrow(p$targets)))

pd2 <- function(idiff, jbeta, ijsk, wh, xmat, beta){
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
             exponents_sum(xmat[h.df,], i.s), # i.s we need this for the state that's in the target
           A=exp(bx_sum),
           # Aprime is the derivative wrt beta-j, which is the x value for that beta
           Aprime=A * xmat[h.df, j.k]) %>%
    ungroup
  
  # for each state, for this person, we need: B_exponent <- beta %*% t(xmat)
  # we also need the sums across states
  # I think? this is the same for all i, j for a given beta so could move out of here and pass it in
  Bsx <- hsstub %>%
    rowwise() %>%
    mutate(bsx=exponents_sum(xmat[h.df, ], s.df), # s.df here we want beta exponent-sum for each given state
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
           # state_betax =mBsx[h.df, j.s],
           state_betax =mBsx[h.df, i.s],
           betax_sum=mBsx_sum[h.df],
           Bprime= -wh * xh_kj * state_betax / betax_sum^2)
  
  gprimeh <- Adf %>% 
    left_join(Bdf, by = "h.df") %>%
    mutate(gprimeh=- xh_ki * (A * Bprime + B * Aprime))
  # - xmat[, i.k] * (A * Bprime + B * Aprime)
  # gprimeh
  element <- sum(gprimeh$gprimeh)
  print(gprimeh)
  element
  
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
  # print(gprimeh)
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


