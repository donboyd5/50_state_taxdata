# https://community.rstudio.com/t/newton-method-with-r/23132/3
# https://math.stackexchange.com/questions/925838/constructing-a-while-loop-in-r-for-newtons-method


g <- function(x) {
  x^3 + 4*x^2 - 10
}

gPrime <- function(x) {
  3*x^2 + 8*x
}

g <- function(x) {x^2 - 1}
gprime <- function(x) {2 * x}

g <- function(xv) {xv[1]^2 + 3*xv[2]^2}
gprime <- function(xv) {c(2*xv[1], 6*xv[2])}
guess <- c(3, 7)

guess <- 1.5

tolerance <- .00001

root <- function(g, gPrime, guess, tolerance) {
  x <- guess
  while (abs(g(x)) > tolerance) {
    x = x - g(x) / gPrime(x)
  }
  x
}


f <- function(xv) {
  xv[1]^2 + 3*xv[2]^2 + 10
  }

g <- function(xv) {c(2*xv[1], 6*xv[2])}

hm <- matrix(c(2, 0,
              0, 6), ncol=2, byrow=TRUE)

h <- function(xv) hm
x0 <- c(3, 7)
x0 <- c(-100, 300)

fmin <- function(x0, f, g, h) {
  x <- x0
  while (sum(abs(g(x))) > .01) {
    x <- x - t(solve(h(x)) %*% g(x))
  }
  x
}

fmin(x0, f, g, h)

fmin2(x0, f)

fmin2 <- function(x0, f) {
  x <- x0
  while (sum(abs(gx)) > .01) {
    gx <- grad(f, x)
    x <- x - t(solve(numDeriv::hessian(func=f, x)) %*% gx)
  }
  x
}

sc2.f <- function(x){
  n <- length(x)
  sum((1:n) * (exp(x) - x)) / n
}

sc2.g <- function(x){
  n <- length(x)
  (1:n) * (exp(x) - 1) / n
}

x0 <- rnorm(5)
hess <- numDeriv::hessian(func=sc2.f, x=x0)
hessc <- hessian(func=sc2.f, x=x0, "complex")
all.equal(hess, hessc, tolerance = .Machine$double.eps)

#  Hessian = Jacobian of the gradient
jac  <- jacobian(func=sc2.g, x=x0)
jacc <- jacobian(func=sc2.g, x=x0, "complex")
all.equal(hess, jac, tolerance = .Machine$double.eps)
all.equal(hessc, jacc, tolerance = .Machine$double.eps)




root(g, gPrime, guess, tolerance)
optim(guess, fn=g, gr=gPrime, method="BFGS")

tibble(x=seq(-2, 2, .01), y=g(x)) %>%
  ggplot(aes(x, y)) +
  geom_line() +
  scale_x_continuous(breaks=seq(-5, 5, 1)) +
  scale_y_continuous(breaks=seq(-10, 50, 1)) +
  geom_hline(yintercept = 0)

library(maxLik)
f <- function(x) -g(x)
maxNR(f, start=1000, print.level=2)

sse <- function(xvec, targs){
  evals <- exp(xvec)
  dist <- (targs - evals)^2
  sum(dist)
}

ssemax <- function(xvec, targs) -sse(xvec, targs)

targs <- c(2, 5, 7)
x <- c(1, 2, 3)
opt <- maxNR(ssemax, start=x, print.level=1, targs=targs)
names(opt)
ssemax(opt$estimate, targs)


# https://en.wikipedia.org/wiki/Newton%27s_method
fx <- function(x) x^2 - 1
fpx <- function(x) 2 * x


tol <- 1e-6
eps <- 1e-14
maxit <- 20

sol <- FALSE
x0 <- .5
for(i in 1:maxit){
  print("")
  print(i)
  print(sprintf("x: %.5e", x0))
  y <- fx(x0)
  print(sprintf("y: %.5e", y))
  yp <- fpx(x0)
  print(sprintf("yp: %.5e", yp))
  print(abs(yp) < eps)
  if(abs(yp) < eps) break
  step <- - y / yp
  print(sprintf("step: %.5e", step))
  x1 <- x0 - y / yp
  print(sprintf("x1: %.5e", x1))
  print(sprintf("abs(x1 - x0) %.5e", abs(x1 - x0)))
  print(abs(x1 - x0) < tol)
  if(abs(x1 - x0) < tol) {sol <- TRUE; break}
  x0 <- x1
}
if(sol==TRUE) sprintf("Solution %.5e", x1) else sprintf("BAD")

optim(200, fn=fx, gr=fpx, method="BFGS")



x0 = 1  # The initial guess
f(x) = x^2 - 2  # The function whose root we are trying to find
fprime(x) = 2 * x  # The derivative of the function
tolerance = 10^(-7)  # 7 digit accuracy is desired
epsilon = 10^(-14)  # Do not divide by a number smaller than this
maxIterations = 20  # Do not allow the iterations to continue indefinitely
solutionFound = false  # Have not converged to a solution yet

for i = 1:maxIterations
y = f(x0)
yprime = fprime(x0)

if abs(yprime) < epsilon  # Stop if the denominator is too small
break
end

global x1 = x0 - y/yprime  # Do Newton's computation

if abs(x1 - x0) <= tolerance  # Stop when the result is within the desired tolerance
global solutionFound = true
break
end

global x0 = x1  # Update x0 to start the process again
end

if solutionFound
println("Solution: ", x1)  # x1 is a solution within tolerance and maximum number of iterations
else
  println("Did not converge")  # Newton's method did not converge
end


f <- function(xv){
  sum(xv[1]^2 + xv[2]^2)
}

g <- function(xv){
  c(2*xv[1], 2*xv[2])
}

h <- function(xv){
  matrix(c(2, 0, 0, 2), ncol=2, byrow=TRUE)
}

xv <- c(300, -90)
for(i in 1:10){
  print(xv <- xv - solve(h(xv)) %*% g(xv) )
}
  
}


root <- function(g, gPrime, guess, tolerance) {
  x <- guess
  while (abs(g(x)) > tolerance) {
    x = x - g(x) / gPrime(x)
  }
  x
}

