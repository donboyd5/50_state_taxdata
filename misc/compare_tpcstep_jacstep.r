

p <- make_problem(h=4, s=3, k=2)
p <- make_problem(h=30, s=5, k=4)
p <- make_problem(h=50, s=10, k=6)
p <- make_problem(h=100, s=20, k=8) # 4 steps jac, 6 steps tpc * 100
p <- make_problem(h=1000, s=20, k=8) # 5 steps tpc * 1000
p <- make_problem(h=2000, s=25, k=10) # findiff 23 secs, serial 42, par8 42; 
p <- make_problem(h=4000, s=30, k=10) # findiff 62, serial 101, par8 80
p <- make_problem(h=6000, s=50, k=20) # par8 821 secs
p <- make_problem(h=10000, s=50, k=30)

p <- make_problem(h=30e3, s=50, k=40)

p <- make_problem(h=30e3, s=50, k=30)

p

p <- pacs
p <- scale_problem_mdn(pacs, scale_goal=100)
p <- scale_problem(pacs, scale_goal=1000)
p <- scale_problem_mdn(pacs, scale_goal=10e3)
p <- scale_problem(pacs, scale_goal=100e3)

set.seed(2345); rbetavec <- runif(p$s * p$k)
rbeta <- vtom(rbetavec, nrow(p$targets))

# unbundle the problem list and create additional variables needed
targets <- p$targets
wh <- p$wh
xmat <- p$xmat

xpx <- t(xmat) %*% xmat
invxpx <- solve(xpx) # TODO: add error check and exit if not invertible

beta0 <- matrix(0, nrow=nrow(targets), ncol=ncol(targets))
delta0 <- get_delta(wh, beta0, xmat) # tpc uses initial delta based on initial beta 

# scale notes ----
# h=100, s=20, k=8  100  4 steps jac, 6 steps tpc * 100
# h=1000, s=20, k=8  5 steps tpc * 1000
# h=2000, s=25, k=10 scale 2000  5 iter # findiff 23 secs, serial 42, par8 42; 
# h=4000, s=30, k=10 scale 4000
# h=6000, s=50, k=20

# start here ----
ebeta <- beta0 # tpc uses 0 as beta starting point
edelta <- delta0 # tpc uses initial delta based on initial beta 

maxiter <- 27
sse_vec <- rep(NA_real_, maxiter)
# iter <- 0
steps <- array(dim=c(2, length(ebeta), maxiter))
sscale <- nrow(xmat)
sscale <- nrow(xmat) * .4

# ACS hks 1000 4 5
# 1 1.7606e-01
# .9 1.5264e-04 1.1 6.3403e+03
# .95 1.9640e-04 .85 4.4871e-03

# ACS hks 10k 4 20
# 1 NAN at 7
# 1 NAN 6 scale 1000
# scale 1000 .75 sscale is good 2.0895e-06
# sscale .85 5.5032e-09
# sscale .90 6.5301e-05
# .95 9.0354e-08


for(iter in 1:maxiter){
  # iter <- iter + 1
  edelta <- get_delta(wh, ebeta, xmat)
  ewhs <- get_weights(ebeta, edelta, xmat)
  ews <- colSums(ewhs)
  ewh <- rowSums(ewhs)
  
  etargets <- t(ewhs) %*% xmat
  d <- targets - etargets
  
  sse <- sum(d^2)
  print(sprintf("iter: %i, sse: %.4e", iter, sse))
  sse_vec[iter] <- sse

  # ad hoc step
  step_tpc <- -(1 / ews) * d %*% invxpx
  step_tpc <- step_tpc * sscale
  
  # gradient step
  # step_gr1 <- numDeriv::grad(sse_fn, as.vector(ebeta), wh=wh, xmat=xmat, targets=targets)
  # step_grad <- vtom(step_gr1, nrow(ebeta)) * .01
  
  # jac <- jacobian(diff_vec, x=as.vector(ebeta), wh=wh, xmat=xmat, targets=targets)
  # ijac <- solve(jac)
  # step_jac <- vtom(as.vector(d) %*% ijac, nrow(ebeta))
  step_jac <- step_tpc

  steps[1, , iter] <- step_tpc
  steps[2, , iter] <- step_jac
  
  ebeta <- ebeta - step_tpc
}

# end of loop ----

sse_vec

(d / targets * 100) %>% round(2)

steps
df <- adply(steps, 1:3) %>%
  setNames(c("type", "ij", "iter", "value")) %>%
  mutate(type=factor(type, levels=1:2, labels=c("tpc", "jac")),
         iter=as.integer(iter),
         ij=as.integer(ij)) %>%
  as_tibble
df
summary(df)

df %>%
  filter(ij==1)

df %>%
  filter(iter <= 10) %>%
  filter(ij <= 8) %>%
  ggplot(aes(iter, value, colour=type)) +
  geom_point() +
  scale_x_continuous(breaks=0:100) +
  facet_wrap(~ij, ncol = 2, scales = "free") +
  ggtitle("Facets are step sizes for a given coefficient, selected iterations")

df %>%
  filter(iter <= 9) %>%
  filter(ij <= 25) %>%
  ggplot(aes(ij, value, colour=type)) +
  geom_point() +
  scale_x_continuous(breaks=0:100) +
  facet_wrap(~iter, ncol = 3, scales = "free") +
  ggtitle("Facets are step sizes for a given iteration, selected coefficients")

df %>%
  filter(iter <= 9) %>%
  filter(ij <= 12) %>%
  pivot_wider(names_from = type) %>%
  ggplot(aes(jac, tpc)) +
  geom_point()+
  facet_wrap(~ij, ncol = 3, scales = "free") +
  geom_abline(slope=1, intercept=0) 

summary(df)
df %>%
  filter(iter <= 9) %>%
  filter(ij <= 20) %>%
  pivot_wider(names_from = type) %>%
  mutate(tjratio=tpc / jac) %>%
  ggplot(aes(iter, tjratio)) +
  geom_point()+
  scale_x_continuous(breaks=0:20) +
  facet_wrap(~ij, ncol = 4, scales = "free") +
  geom_hline(yintercept = 1)




df %>%
  pivot_wider(names_from = type) %>%
  mutate(ratio=jac / tpc) %>%
  arrange(ij, iter)

df %>%
  pivot_wider(names_from = type) %>%
  mutate(ratio=jac / tpc) %>%
  arrange(iter, ij)


as_tibble(steps)

lsteps <- list()
lsteps[[1]] <- c(1, 2, 3)
lsteps[[1]] <- 1
lsteps[1][[1]] <- c(1, 2, 3)
lsteps[1][["tpc"]] <- c(1, 2, 3)

lsteps[1] <- "don"
lsteps[2] <- "jim"
lsteps[[1]][["tpc"]] <- c(1, 2, 3)

best_edelta <- get_delta(ewh, best_ebeta, xmat)
ewhs <- get_weights(best_ebeta, best_edelta, xmat)
ewh <- rowSums(ewhs)
if(scale==TRUE) etargets <- sweep(etargets, 2, problem$scale_factor, "/")
final_step_scale <- step_scale