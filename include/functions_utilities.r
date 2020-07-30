

make_problem <- function(h, k, s){
  #  create a random problem of chosen size
  
  # h: # of households
  # k: # of characteristics per household
  # s: # of states
  
  # returns a list with items created below
  
  # example call: make_problem(8, 2, 3)
  
  set.seed(1234)
  xmat <- matrix(runif(h * k), nrow=h, byrow=TRUE)
  
  set.seed(1234)
  whs <- matrix(runif(h * s, 10, 20), nrow=h, byrow=TRUE)
  
  wh=rowSums(whs)
  ws=colSums(whs)
  
  targets <- t(whs) %*% xmat # s x k
  
  scale_factor <- 1
  
  keepnames <- c("h", "k", "s", "xmat", "whs", "wh", "ws", "targets", "scale_factor")
  problem <- list()
  for(var in keepnames) problem[[var]] <- get(var)
  problem
}

