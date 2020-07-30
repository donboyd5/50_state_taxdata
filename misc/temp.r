# https://bookdown.org/rdpeng/advstatcomp/steepest-descent.html#example-multivariate-normal
# https://bookdown.org/rdpeng/advstatcomp/the-newton-direction.html#newtons-method-in-r

#Log-likelihood, 1st & 2nd derivative
ln <- function(p, Y, R) {
  m <- sum(R==1)
  ln <- m*log(p) - p*sum(Y)
  attr(ln,"gradient") <- m/p-sum(Y)
  attr(ln,"hessian") <- -m/p^2
  ln
}

#Newton-Raphson method
newmle <- function(p, ln, ...) {
  l <- ln(p, ...)
  pnew <- p - attr(l, "gradient") / attr(l,"hessian")
  pnew
}

#Simulate censored data~Exp(1/5)
Y <- rexp(10,1/5)
R <- ifelse(Y>10,0,1)
Y[R==0] = 10

#Plot first derivative of the log-likelihood
x <- seq(0.05, 0.6, 0.01)
plot(x, attr(ln(x, Y, R), "gradient"), type="l",
     xlab=expression(theta), ylab="Score function")
abline(0,0)

#Apply Newton-Raphson iteration 3 times
#Starting value p=0.3
p <- 0.3
p <- 0.4
p <- newmle(p, ln, Y=Y, R=R)
p
p <- newmle(p,ln,Y=Y,R=R)
p
p <- newmle(p,ln,Y=Y,R=R)
p




