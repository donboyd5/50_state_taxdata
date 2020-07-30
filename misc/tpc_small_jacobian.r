
# This program is based largely on the methodology in:
#   Khitatrakun, Surachai, Gordon B T Mermin, and Norton Francis. “Incorporating State Analysis into the 
#   Tax Policy Center’s Microsimulation Model: Documentation and Methodology.” Working Paper, March 2016.
#   https://www.taxpolicycenter.org/sites/default/files/alfresco/publication-pdfs/2000697-Incorporating-State-Analysis-into-the-TPCs-Microsimulation-Model.pdf.


source(here::here("include", "libraries.r"))
library(mvtnorm)
library(numDeriv)

#.. setup ----
# xscale <- 10 and step_scale <- 10 works well with this
h <- 4 # number households
k <- 2 # number characteristics
s <- 3 # number states

xmatrix <- matrix(c(5, 6,
              7, 4,
              9, 2,
              3, 1), nrow=h, byrow=TRUE)
xpx <- t(xmatrix) %*% xmatrix
invxpx <- solve(xpx) # TODO: add error check and exit if not invertible

xxp <- xmatrix %*% t(xmatrix)

whs <- matrix(c(4, 3, 4,
                7, 12, 1,
                20, 6, 4,
                6, 4, 5), nrow=4, byrow=TRUE)

# bigger ----
h <- 20e3
k <- 30
s <- 50

xmatrix <- matrix(runif(h*k), nrow=h, byrow=TRUE)
xpx <- t(xmatrix) %*% xmatrix
invxpx <- solve(xpx) # TODO: add error check and exit if not invertible

whs <- matrix(runif(h*s, 10, 20), nrow=h, byrow=TRUE)

# (xxp <- x %*% t(x))
# (xpx <- t(x) %*% x)
(wh=rowSums(whs))
(ws=colSums(whs))

(targets <- t(whs) %*% xmatrix) # s x k

#.. end setup ----


calc_delta <- function(hh_weights, beta, xmatrix){
  beta_x <- exp(beta %*% t(xmatrix))
  log(hh_weights / colSums(beta_x))
}


calc_weights <- function(beta, delta, xmatrix){
  # calculate all weights
  beta_x <- beta %*% t(xmatrix)
  # add delta to every row of beta_x and transpose
  beta_xd <- apply(beta_x, 1, function(mat) mat + delta) 
  exp(beta_xd)
}



f <- function(betavec, wh, xmatrix, targets){
  # return the distances from targets as a vector
  beta <- matrix(betavec, nrow=nrow(targets), byrow=FALSE)
  delta <- calc_delta(wh, beta, xmatrix)
  whs <- calc_weights(beta, delta, xmatrix)
  etargets <- t(whs) %*% xmatrix
  d <- targets - etargets
  as.vector(d)
}


get_result <- function(){
  result <- list()
  result$iter <- i
  result$sse <- sse
  result$d <- d
  result$ebeta <- ebeta
  result$edelta <- edelta
  result$whs <- calc_weights(ebeta, edelta, xmatrix)
  result
}

jac <- function(ewhs, xmatrix){
  x2 <- xmatrix * xmatrix
  ddiag <- - t(ewhs) %*% x2 # note the minus sign in front
  diag(as.vector(ddiag)) 
}

# jac(ewhs, xmatrix)


ebeta <- matrix(0, nrow=s, ncol=k) # tpc uses 0 as beta starting point
edelta <- calc_delta(wh, ebeta, xmatrix) # tpc uses initial delta based on initial beta 

maxiter <- 1000
for(i in 1:maxiter){
  ewhs <- calc_weights(ebeta, edelta, xmatrix)
  ews <- colSums(ewhs)
  ewh <- rowSums(ewhs)
  
  etargets <- t(ewhs) %*% xmatrix
  d <- targets - etargets
  sse <- sum(d^2)
  if(i <= 20 | i %% 50==0) print(sse)
  if(sse < 1e-6) {
    # exit if good
    result <- get_result()
    break
  }
  
  if(FALSE) {
    jval <- jacobian(f, x=as.vector(ebeta), wh=wh, xmatrix=xmatrix, targets=targets)
    jval <- jac(ewhs, xmatrix)
    step <- solve(jval) %*% as.vector(d)
    step <- matrix(step, nrow=nrow(d), byrow=FALSE)
  } else step <- -(1 / ews) * d %*% invxpx * 10000
  
  ebeta <- ebeta - step
  edelta <- calc_delta(ewh, ebeta, xmatrix)
  if(i==maxiter) {result <- get_result(); break}
}

result

whs
ewhs %>% round(1)
ws; ews %>% round(1)
wh; ewh

wh %*% xxp
ewhs[, 1] %*% xxp
ewhs[, 1] %*% xmatrix[, 1] %*% xmatrix[, 1]
sum(ewhs[, 1]) %*% xmatrix[, 1] %*% xmatrix[, 1]
sum(ewhs[, 1]) %*% xmatrix[, 1] %>% sum

ewhs[, 1] * (xmatrix[, 1] * xmatrix[, 1])

sum(ewhs[, 1] * xmatrix * xmatrix)
sum(ewhs[, 1]) * (t(xmatrix) %*% xmatrix)
sum(ewhs[, 1]) * (xmatrix %*% t(xmatrix))

colSums(ewhs) * (t(xmatrix) %*% xmatrix)
rowSums(ewhs)


(jnum <- jacobian(f, x=as.vector(ebeta), wh=wh, xmatrix=xmatrix, targets=targets))

ewhs[, 1] * xxp %>% rowSums
ewhs[, 1] * xxp %>% colSums
diag(ewhs[, 1] * xxp) %>% sum
ewhs[, 1] * xpx

ewhs[, 1] %*% xmatrix[, 1] %*% xmatrix[, 1] %>% sum

sweep(xxp, 2, ewhs[, 1], "*")

