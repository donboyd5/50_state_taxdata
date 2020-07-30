
library(nleqslv) # must use global dbldog is good; newton better than broyden unless too big I think

get_dweights <- function(prob, goal=100){
  # difference weights - a weight to be applied to each target in the difference function so that
  # it hits its goal
  dw <- prob$targets / goal
  dw <- ifelse(prob$targets!=0, 100 / prob$targets, 1)
  dw
}



p <- make_problem(h=4, s=3, k=2)
p <- make_problem(h=30, s=5, k=4)
p <- make_problem(h=50, s=10, k=6)
p <- make_problem(h=100, s=20, k=8) # 4 steps jac, 6 steps tpc * 100
p <- make_problem(h=1000, s=20, k=8) # 5 steps tpc * 1000
p <- make_problem(h=2000, s=25, k=10) # findiff 23 secs, serial 42, par8 42; 
p <- make_problem(h=4000, s=30, k=10) # findiff 62, serial 101, par8 80
p <- make_problem(h=6000, s=50, k=20) # par8 821 secs
p <- make_problem(h=10000, s=50, k=30)

prob
badbeta <- betavec

p


dw * p$targets


betavec <- runif(length(p$targets), -100, 100)
betavec <- rnorm(length(p$targets))

dw <- get_dweights(p$targets)
betavec <- rep(0, length(p$targets))

sse_fn(betavec, wh=p$wh, xmat=p$x, targets=p$targets, dweights = dw)

res <- tpc(betavec, p$wh, p$x, p$targets, dweights=dw, ssfactor=1 / nrow(p$targets), maxiter=10)
res <- tpc(betavec, p$wh, p$x, p$targets, dweights=dw, ssfactor=3, maxiter=10)
res$delta
res

t(res$weights) %*% p$xmat
p$targets

reldiff(betavec, p) %>% round(2)
reldiff(as.vector(res$beta), p) %>% round(2)

betavec <- rep(0, length(p$targets)); wh <- p$wh; xmat <- p$xmat; targets <- p$targets; dweights <- as.vector(dw); testiter <- 5
length(p$targets)

start_search <- function(betavec=NULL, wh, xmat, targets, dweights=rep(1, length(targets)), testiter){
  # binary search for the best scale_factor
  # repeatedly run tpc looking for best scale factor and starting beta
  minscale <- 1e-3
  maxscale <- length(targets) / 2
  
  vlen <- 20
  sse_vec <- rep(NA, vlen)
  scale_vec <- rep(NA, vlen)
  beta_list <- vector(mode = "list", length = vlen)
  
  # define the first 5 searches, before beginning binary search
  scale_vec[1:5] <- c(1, minscale, maxscale, (1 + minscale) / 2, (1 + maxscale) / 2)
  for(i in 1:5){
    print(sprintf("testing step scale %.2f", scale_vec[i]))
    resx <- tpc(betavec, wh, xmat, targets, dweights, ssfactor=scale_vec[i], maxiter=testiter, quiet=TRUE)
    sse_vec[i] <- resx$sse
    beta_list[[i]] <- resx$beta
  }
  
  i_best <- which.min(sse_vec)
  beta_best <- beta_list[[i_best]]
  # print(i_best); print(i_nextbest) #; print(scale_vec); print(sse_vec)
  # i <- 5
  
  for(i in 6:vlen){
    # i <- i + 1
    # print(sprintf("iter %i", i))
    best <- scale_vec[i_best]
    nextbest_sse <- min(sse_vec[-i_best], na.rm=TRUE)
    i_nextbest <- min(which(sse_vec == nextbest_sse))
    nextbest <- scale_vec[i_nextbest]
    # print(sprintf("best scale %.2f nextbest scale %.2f", best, nextbest))
    # print(sprintf("best sse %.2f nextbest sse %.2f actual %.2f", sse_vec[i_best], sse_vec[i_nextbest], nextbest_sse))
    
    ssfactor <- (best + nextbest) / 2
    scale_vec[i] <- ssfactor
    print(sprintf("testing step scale %.2f", scale_vec[i]))
    
    res <- tpc(betavec, wh, xmat, targets, dweights, ssfactor, maxiter=testiter, quiet=TRUE)
    sse_vec[i] <- res$sse
    i_best <- which.min(sse_vec) # get the FIRST best
    if(i==i_best) beta_best <- res$beta
    
    i_allbest <- which(sse_vec == min(sse_vec, na.rm = TRUE)) # get multiple minimums if there are multiples
    i_nextbest <- which.min(sse_vec[-i_allbest]) # make sure we remove any multiple minimums before getting next-best
    
    # if(scale_vec[i]==scale_vec[i - 1]) break
  }
  
  result <- list()
  result$vlen <- vlen
  result$i_best <- i_best
  result$scale_best <- scale_vec[i_best]
  result$sse_best <- sse_vec[i_best]
  result$scale_vec <- scale_vec
  result$sse_vec <- sse_vec
  result$beta_best <- as.vector(beta_best)
  result
}

p <- pacs
names(p)
p$targets


# approach #1: replace 0 targets with something else ----
p$targets <- ifelse(p$targets==0, 10e3, p$targets)
pacs$targets - p$targets

# alternatively, drop any column that has a zero and redefine p accordingly ----
zcols <- apply(p$targets, 2, function(x) any(x==0))
sum(zcols)
which(zcols)
p$targets <- p$targets[, -which(zcols)]
p$xmat <- p$xmat[, -which(zcols)]
p$k <- p$k - sum(zcols)
p

# continue on ---
dw <- get_dweights(p)
betavec <- rep(0, length(p$targets))

ss <- start_search(betavec, p$wh, p$xmat, p$targets, dw, testiter=5)
ss
cbind(ss$scale_vec, ss$sse_vec)

tmp <- tpc(ss$beta_best, p$wh, p$xmat, p$targets, dw, ssfactor=ss$scale_best, maxiter=30)
tmp <- tpc(betavec, p$wh, p$xmat, p$targets, dw, ssfactor=ss$scale_best, maxiter=30)
tmp <- tpc(betavec, p$wh, p$xmat, p$targets, dw, ssfactor=1, maxiter=30)
tmp <- tpc(betavec, p$wh, p$xmat, p$targets, dw, ssfactor=1 / nrow(p$xmat), maxiter=30)

tmp2 <- tpc(ss$beta_best, p$wh, p$xmat, p$targets, dw, ssfactor=10, maxiter=30)

tmp3 <- tpc(betavec, p$wh, p$xmat, p$targets, dw, ssfactor=8, maxiter=5)

out1 <- nls.lm(par=as.vector(tmp$beta), fn=diff_vec, 
              control=nls.lm.control(maxiter=50, nprint=1, factor=100),
              wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)
names(out1)
out1$message
out1$niter
out1$rsstrace
names(p)

p$h; p$k; p$s
t1 <- p$targets
t1 %>% round(1)
t2 <- etargs_mat(betavec, p$wh, p$xmat, p$s)
(t2 - t1) %>% round(1)
((t2 - t1) / t1 * 100) %>% round(1)


t2a <- etargs_mat(ss$beta_best, p$wh, p$xmat, p$s)
(t2a - t1) %>% round(1)
((t2a - t1) / t1 * 100) %>% round(1)


t2b <- etargs_mat(as.vector(tmp$beta), p$wh, p$xmat, p$s)
(t2b - t1) %>% round(1)
((t2b - t1) / t1 * 100) %>% round(1)


t3 <- etargs_mat(out3$x, p$wh, p$xmat, p$s)
(t3 - t1) %>% round(1)
((t3 - t1) / t1 * 100) %>% round(1)

t4 <- etargs_mat(out2$par, p$wh, p$xmat, p$s)
(t4 - t1) %>% round(1)
((t4 - t1) / t1 * 100) %>% round(1)

t5 <- etargs_mat(out4b$x, p$wh, p$xmat, p$s)
(t5 - t1) %>% round(1)
((t5 - t1) / t1 * 100) %>% round(1)

t6 <- etargs_mat(dfs$par, p$wh, p$xmat, p$s)
(t6 - t1) %>% round(1)
((t6 - t1) / t1 * 100) %>% round(1)


a <- proc.time()
out2 <- nls.lm(par=betavec, fn=diff_vec, # betavec
              control=nls.lm.control(maxiter=50, nprint=1, factor=100),
              wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)
b <- proc.time()
b - a


library(rootSolve)
a <- proc.time()
mr <- multiroot(f=diff_vec, start=betavec, maxiter=50, verbose=TRUE, useFortran=TRUE,
                wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=as.vector(dw))
b <- proc.time()
b - a
# names(mr)
mr$root
mr$f.root
mr$iter


library(BB)
# sane and dfsane do not seem able to solve this
bb1 <- proc.time()
# dfs <- sane(par=betavec, fn=diff_vec, method=1, control=list(M=15, NM=FALSE, maxit=400, trace=TRUE),
#               wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=as.vector(dw))
dfs <- BBsolve(par=betavec, fn=diff_vec, 
               wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=as.vector(dw))
bb2 <- proc.time()
bb2 - bb1
dfs$message
dfs$cpar

diff_vec(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=as.vector(dw))

names(out2)
out2$message
out2$rsstrace
out2$fvec

out2a <- nls.lm(par=betavec, fn=diff_vec,
               control=nls.lm.control(maxiter=50, nprint=1, factor=100),
               wh=p$wh, xmat=p$xmat, targets=p$targets)


out3 <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
                method = c("Newton"),
                global = c("dbldog"),
                xscalm = c("auto"),
                jacobian=FALSE,
                control=list(maxit=50, trace=1, allowSingular=TRUE, ftol=1e-1)) 
# ftol is % error as I formulate this with dweights, so 1e-2 means when all are < 0.1% we can stop

out3a <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
                method = c("Newton"),
                global = c("none"),
                xscalm = c("gline"),
                jacobian=FALSE,
                control=list(maxit=20, trace=1, allowSingular=TRUE))

out4b <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                 wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
                 method = c("Broyden"),
                 global = c("dbldog"),
                 xscalm = c("auto"),
                 jacobian=FALSE,
                 control=list(maxit=2000, trace=1, allowSingular=TRUE))

out4 <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=as.vector(dw),
                method = c("Newton"),
                global = c("dbldog"),
                xscalm = c("auto"),
                jacobian=FALSE,
                control=list(maxit=30, trace=1, allowSingular=TRUE))

out4a <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
                method = c("Newton"),
                global = c("hook"),
                xscalm = c("auto"),
                jacobian=FALSE,
                control=list(maxit=30, trace=1, allowSingular=TRUE)) #
out4a$fvec
out4a$message

out4a_xdw <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                 wh=p$wh, xmat=p$xmat, targets=p$targets,
                 method = c("Newton"),
                 global = c("hook"),
                 xscalm = c("auto"),
                 jacobian=FALSE,
                 control=list(maxit=30, trace=1, allowSingular=TRUE)) #

out4b <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                 wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
                 method = c("Newton"),
                 global = c("cline"),
                 xscalm = c("auto"),
                 jacobian=FALSE,
                 control=list(maxit=30, trace=1, allowSingular=TRUE)) # best
names(out4b)
out4b$fvec

out4b <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                 wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
                 method = c("Newton"),
                 global = c("qline"),
                 xscalm = c("auto"),
                 jacobian=FALSE,
                 control=list(maxit=5, trace=1, allowSingular=TRUE))

out4c <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                 wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
                 method = c("Newton"),
                 global = c("gline"),
                 xscalm = c("auto"),
                 jacobian=FALSE,
                 control=list(maxit=5, trace=1, allowSingular=TRUE))

out4d <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                 wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
                 method = c("Newton"),
                 global = c("none"),
                 xscalm = c("auto"),
                 jacobian=FALSE,
                 control=list(maxit=5, trace=1, allowSingular=TRUE))


out4b <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                 wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
                 method = c("Broyden"),
                 global = c("cline"),
                 xscalm = c("auto"),
                 jacobian=FALSE,
                 control=list(maxit=30, trace=1, allowSingular=TRUE))






out5 <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                wh=p$wh, xmat=p$x, targets=p$targets, dweights=dw,
                method = c("Newton"),
                global = c("dbldog"),
                xscalm = c("auto"),
                jacobian=FALSE,
                control=list(maxit=50, trace=1, allowSingular=TRUE))




out6 <- nleqslv(x=tmp$beta, fn=diff_vec, jac=NULL,
                wh=p$wh, xmat=p$x, targets=p$targets, dweights=dw,
                method = c("Broyden"),
                global = c("dbldog"),
                xscalm = c("auto"),
                jacobian=FALSE,
                control=list(maxit=1000, trace=1, allowSingular=TRUE))

out5 <- nleqslv(x=betavec, fn=diff_vec, jac=NULL,
                wh=p$wh, xmat=p$x, targets=p$targets, dweights=dw,
                method = c("Broyden"),
                global = c("dbldog"),
                xscalm = c("auto"),
                jacobian=FALSE,
                control=list(maxit=1000, trace=1, allowSingular=TRUE))


tpc <- function(betavec=NULL, wh, xmat, targets, dweights=rep(1, length(targets)), ssfactor, maxiter, quiet=FALSE){
  
  # items to store
  sse_vec_tpc <- rep(NA, maxiter)
  
  if(is.null(betavec)) betavec <- rep(0, length(targets))
  ebeta <- vtom(betavec, nrow(targets))
  
  best_ebeta <- ebeta
  best_sse <- Inf
  
  xpx <- t(xmat) %*% xmat
  invxpx <- solve(xpx) # TODO: add error check and exit if not invertible
  
  step_scale <- nrow(xmat) * ssfactor
  # step_scale <- 1
  # iter <- 0
  
  for(iter in 1:maxiter){
    # iter <- iter + 1
    edelta <- get_delta(wh, ebeta, xmat)
    ewhs <- get_weights(ebeta, edelta, xmat)
    ews <- colSums(ewhs)
    ewh <- rowSums(ewhs)
    
    # etargets <- t(ewhs) %*% xmat # not needed
    betavec <- as.vector(ebeta)
    
    d <- diff_vec(betavec, wh, xmat, targets, dweights)
    
    worst_diff <- max(abs(d), na.rm=TRUE)
    sse <- sum(d^2)
    sse_vec_tpc[iter] <- sse
    if(is.na(sse)) break # bad result, end it now, we have already saved the prior best result
    
    if(!quiet) print(sprintf("iter: %i, sse: %.4e, worst diff: %.4e", iter, sse, worst_diff))
    
    rel_err <- ifelse(targets==0, NA, abs(d / targets))
    max_rel_err <- max(rel_err, na.rm=TRUE)
    
    if(iter >= 2) {sse_rel_change <- sse_vec_tpc / lag(sse_vec_tpc) - 1} else sse_rel_change <- NA_real_
    # iter <- 5
    # test2 <- ifelse(iter >= 5, !any(sse_rel_change[iter - 0:2] < -.01), FALSE)
    # test2
    # any(sse_rel_change[c(4, 5, 6)] < -.01)
    
    best_sse <- min(sse_vec_tpc, na.rm=TRUE)
    if(sse==best_sse) best_ebeta <- ebeta
    if(sse < 1e-8) break
    
    prior_sse <- sse
    
    step_tpc <- -(1 / ews) * vtom(d, nrow(targets)) %*% invxpx
    
    ebeta <- ebeta - step_tpc * step_scale
  }
  result <- list()
  result$sse <- best_sse
  result$beta <- best_ebeta
  result$delta <- get_delta(wh, result$beta, xmat)
  result$weights <- get_weights(result$beta, result$delta, xmat)
  return(result)
}


