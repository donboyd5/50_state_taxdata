

dx2x <- deriv(~ x^2, "x") ; dx2x


trig.exp <- expression(sin(cos(x + y^2)))
( D.sc <- D(trig.exp, "x") )
all.equal(D(trig.exp[[1]], "x"), D.sc)

deriv((y ~ sin(cos(x) * y)), c("x","y"), func = TRUE)

djb <- expression(t - a * exp(a*x + b*x + c*y + d))
djb <- expression(t - (a * exp(a*x + c*y + d) - (b * exp(f* x + g*z + h))))
D(djb, "x")

djb <- expression(t - (a * exp(a*beta + c*y + d) - (b * exp(f* x + g*z + h))))
D(djb, "beta")

# weight hh 1, state 1
w11 <- expression(exp(b11*x11 + b12*x12 + c1))
D(w11, "b11")

# weight hh 1, state 2
w12 <- expression(exp(b21*x11 + b22*x12 + c1))
D(w12, "b11")

ewhs
xmatrix
x2 <- xmatrix * xmatrix
x2 %*% ewhs[, 1]
ewhsl <- as.vector(ewhs)
db <- x2[, 1] * ewhs[, 1]

t(ewhs[, 1]) %*% x2[, 1] # this is fopd s=1, k=1
t(ewhs) %*% x2[, 1] # these are the fopds for s=1, 2, 3; k=1
ddiag <- t(ewhs) %*% x2
ddiag
diag(as.vector(ddiag))


sum(db)

a <- expression(x + y)
b <- expression(a^2)
D(b, "x")

c <- expression((x+y)^2)
D(c, "x")

xmatrix[, 1] * xmatrix[, 1] %*% t(ewhs[, 1])
xmatrix[, 1] * xmatrix[, 1] %*% t(ewhs[, 1])

# t: s, k
# w: h, s
# x: h, k
# targets for state 1, characteristics 1 and 2
d <- expression(
  (t11 - (w11 * x11 +
            w12 * x21)) +
    (t12 - (w11 * x12 +
              w21 * x22)))
D(d, "w11")

db <- expression(
  (t11 - (exp(b11*x11 + b12*x12 + c1) * x11 +
          exp(b11*x21 + b12*x22 + c2) * x21)) +
  (t12 - (exp(b11*x11 + b12*x12 + c1) * x12 +
          exp(b11*x21 + b12*x22 + c2) * x22)))
D(db, "b11")

db2 <- expression(
  ((t11 - (exp(b11*x11 + b12*x12 + c1) * x11 +
            exp(b11*x21 + b12*x22 + c2) * x21)) +
    (t12 - (exp(b11*x11 + b12*x12 + c1) * x12 +
              exp(b11*x21 + b12*x22 + c2) * x22)))^2)
D(db2, "b11")



tmp <- poisson_weights(hweights, targets, x)
str(tmp)

a <- proc.time()

b <- proc.time()


ba <- b - a

print(ba)

ba
as.numeric(ba)[3]

bas <- as.character(print(ba))
bas 


library(cOde)
jacobianSymb(c(A="A*B", B="A+B"))
jacobianSymb(c(A="A*B", B="A+B"))
jacobianSymb(c(x="A*B", y="A+B"))

jacobianSymb(c(x="A*B", y="A+B"), c("A", "B"))

jacobianSymb(c(A="A*B"))

jacobianSymb(c(d="x^2"), variables=c("5"))

jacobianSymb(c(d="A, B, C", y="A+B"), c("A", "B", "C"))


jacobianSymb(c(t="targets", v="values"), c(1, 2, 3))


beta_x <- exp(beta %*% t(x))
log(wh / colSums(beta_x))

# hh 1, state 2 weight
w_h1s2 <- expression(exp(b21*x11 + b22*x12 + log(wh1 / (b21 * x11 + b22 * x12))))
D(w_h1s2, "b21") # w_h1s2 * 


# bsk s 1:3, k 1:2
w_h1s2 <- expression(exp(b21*x11 + b22*x12 + c1))
D(w_h1s2, "b21") # w_h1s2 * x11
D(w_h1s2, "b22") # w_h1s2 * x12

w_h1s3 <- expression(exp(b31*x11 + b32*x12 + c1))
D(w_h1s3, "b31") # w_h1s3 * x11
D(w_h1s3, "b32") # w_h1s3 * x12

w_h4s3 <- expression(exp(b31*x41 + b32*x42 + c1))
D(w_h4s3, "b31") # w_h1s3 * x11



b <- expression((t1 - (x11*exp(b1*x11) + x21*exp(b1*x21))^2) +
                (t2 - (x11*exp(b2*x12) + x21*exp(b2*x22))^2))

b <- expression((t1 - (x11*exp(b1*x11) + x21*exp(b1*x21))^2) +
                  (t2 - (x11*exp(b2*x12) + x21*exp(b2*x22))^2))

D(b, "b1")
D(b, "b2")

deriv(b, namevec="b1")

me <- expression(b1 * (x %*% x))
deriv(me, namevec="b1")

library(Deriv)
Deriv("sin(x^2) * y", "x") 
Deriv(me, "b1")

Deriv(~solve(matrix(c(1, x, x**2, x**3), nrow=2, ncol=2)))

bx <- expression(x * exp(b1*x + c))
Deriv(bx, "b1")

bf <- function(x, b1, c){
  exp(b1*x + c)
}

df <- function(x, y) (x + y)
Deriv(df, "x")

df <- Deriv(bf, "b1")
df
df(1, 8, 1)
1*exp(1*8 + 1)


df <- function(x) x^2
Deriv(df)

df <- function(x) x^2
dg <- Deriv(df)

df(7)
dg(7)


df <- function(x, n){
  s <- sum(x)
  x
}
dg <- Deriv(df)
dg


a <- expression(whs*xh)
D()


bx <- expression(xh * exp(b1*xh + c))
Deriv(bx, "b1")

e1 <- expression(whs * xh1 + whs * xh2)
Deriv(e1, "whs")


f1b <- function(x) exp(x)
f1b <- expression(xh1k1*exp(xh1k1*b1 + xh1k2*b2))
f1b <- expression(xh1k1*exp(xh1k1*b1 + xh1k2*b2 + c))
# expression(xh1k1^2 * exp(b1 * xh1k1 + b2 * xh1k2 + c))
Deriv(f1b, "b1")



x <- c(5, 10)
w <- 3
xh <- c(11, 20)
w %*% xh
w * xh

t <- c(200, 100)
w <- c(1, 2, 3)
xh <- matrix(c(11, 20, 
               30, 15,
               40, 10), ncol=2, byrow=TRUE)
xh
et <- w %*% xh
et
d <- t - et
d

# sk
k <- 2
s <- 3
h <- 4
ks <- k * s

# s x k
targs <- matrix(c(1, 2,
                  3, 4,
                  5, 6), ncol=k, byrow=TRUE)
targs

# h x s
whs <- matrix(c(1, 2, 3,
                4, 5, 6,
                6, 7, 8,
                9, 10, 11), ncol=s, byrow=TRUE)
whs

# 1 x h
(wh <- rowSums(whs))
# 1 x s
(ws <- colSums(whs))

xhs <- matrix(c(1, 2, 
               3, 4,
               5, 6,
               7, 8), ncol=2, byrow=TRUE)
xhs

etargs <- t(whs) %*% xhs
etargs


dm <- targs - etargs
dm


