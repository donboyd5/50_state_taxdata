names(inputs2)
inputs2$cc_addup_sparse

Amat <- matrix(0, nrow=nrow(inputs2$weight_totals), ncol=inputs2$n_variables)
dim(Amat)
ij <- inputs2$cc_addup_sparse$j
ii <- inputs2$cc_addup_sparse$ipid
max(ij)
max(ii)
ipairs <- cbind(ii, ij)
str(ipairs)
max(ipairs[, 1])
max(ipairs[, 2])
Amat[ipairs] <- inputs2$cc_addup_sparse$iweight_state
Amat[1:12, 1:12]

calc <- Amat %*% inputs2$x0
b <- inputs2$weight_totals$weight_total
sum(calc - b)
meq <- 1:length(b)

t1 <- proc.time()
test <- spg(par=inputs2$x0, fn=eval_f_wfs, gr=eval_grad_f_wfs, project="projectLinear", 
            projectArgs=list(A=Amat, b=b, meq=meq), control=list(maxit=10), inputs=inputs2)
t2 <- proc.time()
t2 - t1


fn <- function(x) (x[1] - 3/2)^2 + (x[2] - 1/8)^4

gr <- function(x) c(2 * (x[1] - 3/2) , 4 * (x[2] - 1/8)^3)

# This is the set of inequalities
# x[1] - x[2] >= -1
# x[1] + x[2] >= -1
# x[1] - x[2] <= 1
# x[1] + x[2] <= 1

# The inequalities are written in R such that:  Amat %*% x  >= b 
Amat <- matrix(c(1, -1, 1, 1, -1, 1, -1, -1), 4, 2, byrow=TRUE)
b <- c(-1, -1, -1, -1)
meq <- 0  # all 4 conditions are inequalities

p0 <- rnorm(2)
spg(par=p0, fn=fn, gr=gr, project="projectLinear", 
    projectArgs=list(A=Amat, b=b, meq=meq))

meq <- 1  # first condition is now an equality
spg(par=p0, fn=fn, gr=gr, project="projectLinear", 
    projectArgs=list(A=Amat, b=b, meq=meq))

projectLinear(par=p0, A=Amat, b=b, meq=meq)

spg(par=p0, fn=fn, gr=gr, project="projectLinear", 
    projectArgs=list(A=Amat, b=b, meq=meq))


# box-constraints can be incorporated as follows:
# x[1] >= 0
# x[2] >= 0
# x[1] <= 0.5
# x[2] <= 0.5

Amat <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), 4, 2, byrow=TRUE)
b <- c(0, 0, -0.5, -0.5)

meq <- 0
spg(par=p0, fn=fn, gr=gr, project="projectLinear", 
    projectArgs=list(A=Amat, b=b, meq=meq))

# Note that the above is the same as the following:
spg(par=p0, fn=fn, gr=gr, lower=0, upper=0.5)


# An example showing how to impose other constraints in spg()

fr <- function(x) { ## Rosenbrock Banana function
  x1 <- x[1] 
  x2 <- x[2] 
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2 
} 

# Impose a constraint that sum(x) = 1

proj <- function(x){ x / sum(x) }

spg(par=runif(2), fn=fr, project="proj") 

# Illustration of the importance of `projecting' the constraints, rather 
#   than simply finding a feasible point:

fr <- function(x) { ## Rosenbrock Banana function 
  x1 <- x[1] 
  x2 <- x[2] 
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2 
} 
# Impose a constraint that sum(x) = 1 

proj <- function(x){ 
  # Although this function does give a feasible point it is 
  #  not a "projection" in the sense of the nearest feasible point to `x'
  x / sum(x) 
} 

p0 <- c(0.93, 0.94)  

# Note, the starting value is infeasible so the next 
#   result is "Maximum function evals exceeded"

spg(par=p0, fn=fr, project="proj") 

# Correct approach to doing the projection using the `projectLinear' function

spg(par=p0, fn=fr, project="projectLinear", projectArgs=list(A=matrix(1, 1, 2), b=1, meq=1)) 

# Impose additional box constraint on first parameter

p0 <- c(0.4, 0.94)    # need feasible starting point

spg(par=p0, fn=fr,  lower=c(-0.5, -Inf), upper=c(0.5, Inf),
    project="projectLinear", projectArgs=list(A=matrix(1, 1, 2), b=1, meq=1)) 


