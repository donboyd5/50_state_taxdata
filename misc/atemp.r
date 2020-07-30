library(abind)

vector1 <- c(5,9,3)
vector2 <- c(10,11,12,13,14,15)

# Take these vectors as input to the array.
new.array <- array(c(vector1,vector2),dim = c(3,3,2))
print(new.array)

# Use apply to calculate the sum of the rows across all the matrices.
result <- apply(new.array, c(1), sum)
print(result)

apply(new.array, 1, sum)
apply(new.array, 2, sum)
apply(new.array, 3, sum)

apply(new.array, 1:2, sum)


a <- array(1:24, dim=c(3, 2, 4))
apply(a, 1, sum)
apply(a, 2, sum)
apply(a, 3, sum)

a
apply(a, 1:2, sum)

library(multiApply)

A <- array(1:20, c(5, 2, 2))
B <- array(1:20, c(5, 2, 2))

A
B

D <- Apply(data = list(A, B), 
           target_dims = c(2, 3), 
           fun = "%*%")$output1

D2 <- Apply(data = list(A, B),
           target_dims = c(2, 3),
           fun = "%*%",
           ncores = 4)$output1

Apply(data = list(A, B), 
      target_dims = c(2, 3), 
      fun = "%*%")$output1

dims <- c(store = 4, item = 6, day = 8)
sales_amount <- array(rnorm(prod(dims)), dims)

dims <- c(store = 4, item = 6)
sales_price <- array(rnorm(prod(dims)), dims)
dim(sales_price)

income_function <- function(x, y) {
  # Expected inputs:
  #  x: array with dimensions (item, day) -- ie a single day
  #  y: price point vector with dimension (item)
  sum(rowSums(x) * y) # multiply sum of amount over days, by price, for a single store
}

income_function2 <- function(sa, sp) {
  # store = 4, item = 6, day = 8
  # Expected inputs:
  #  sa: array with dimensions (item, day) -- ie sales amounts for a single store, by item and day
  #  sp: price point vector with dimension (item)
  rowSums(sa) * sp # multiply sum of amount over days, by price, for a single store
  # sa is 6 items x 8 days so rowsums is 6 rows, sum of amount, over days for each item, for a single store
  # and sp is 1 x 6, price for a single item
  # returns matrix with (items x store) each cell is sum of items x price
}

income_function2 <- function(sa, sp) {
  print("sa"); print(sa)
  # store = 4, item = 6, day = 8
  # Expected inputs:
  #  sa: array with dimensions (item, day) -- ie sales amounts for a single store, by item and day
  #  sp: price point vector with dimension (item)
  rowSums(sa) * sp # multiply sum of amount over days, by price, for a single store
  # sa is 6 items x 8 days so rowsums is 6 rows, sum of amount, over days for each item, for a single store
  # and sp is 1 x 6, price for a single item
  # returns matrix with (items x store) each cell is sum of items x price
  return(NULL)
}


# income
Apply(data = list(sales_amount, sales_price),
                target_dims = list(c('item', 'day'), 'item'), # sa has a store with item amount and day, sp has price - what varies is store
                income_function2)

atemp <- Apply(data = list(nwhx, Bxs_div_BXsum2),
             target_dims="h",
             fun = "*",
             ncores = ncpar)$output1


bpm1 <- Apply(data = list(Bxs_div_BXsum2, nwhx),
               target_dims="h",
               fun = "*",
               ncores = ncpar)$output1
bpm <- Apply(bpm1, margins=1, fun=c,
             ncores = ncpar)$output1 %>% t

bpm
Bprime_mat




dim(income$output1)
# store


xm <- xmat
am <- Amat
dim(xm) <- c(h=nrow(xm), k=ncol(xm))
xm
dim(xm)

dim(am) <- c(h=nrow(am), s=ncol(am))
dim(am)

# multiply every 
# lAprime
Amat * xm[, 1]

Apply(data=list(am, xm), 
      target_dims = list(c("h", "s"), c("h", "k")),
      function(x, y) x * y)


am
xm
Apply(data = list(am, xm), 
      margins=2,
      fun = "*")$output1

aAprime <- Apply(data = list(am, xm), 
      margins=list("s", "k"), # multiply the "s" dimension (columns) in Amat by the "k" dimention (columns) in xmat
      fun = "*")$output1




Apply(data = list(am, xm), 
      target_dims=list(c("h", "s")),
      fun = "*")$output1

Apply(data = list(am, xm), 
      target_dims = c(2, 3), 
      fun = "*")$output1


#Change in the rate of exceedance for two arrays, with different
#dimensions, for some matrix of exceedances.
data <- list(array(rnorm(1000), c(5, 10, 20)),
             array(rnorm(500), c(5, 10, 10)),
             array(rnorm(50), c(5, 10)))
test_fun <- function(x, y, z) {
  ((sum(x > z) / (length(x))) /
     (sum(y > z) / (length(y)))) * 100
}
test <- Apply(data, target = list(3, 3, NULL), test_fun)


A <- array(c(rep(1,20), rep(2,20), rep(3,20)),dim = c(10,2,3))
B <- matrix(c(1:10), nrow = 2)
# multiply each A[,,i]%*%B
A
B


C <- array(NA, dim=c(nrow(A), ncol(B), 3))
C[] <- apply(A, 3, function(x) x%*%B)


A <- array(c(rep(1,20), rep(2,20), rep(3,20)),dim = c(10,2,3))
B <- matrix(c(1:10), nrow = 2)
# multiply each A[,,i]%*%B
A
B


C <- array(NA, dim=c(nrow(A), ncol(B), 3))
C[] <- apply(A, 3, function(x) x%*%B)


C <- array(NA, dim=c(nrow(A), ncol(B), 3))
C[] <- apply(A, 3, function(x) x * B)


dim(aABprime)
f <- function(m1, m2){
  #print("new"); print(m1); print(m2)
  # print(m1 %*% m2)
  return(m1 %*% m2)
}
djb <- Apply(data = list(xmat, aABprime),
      target_dims=list(c("k", "h"), c("h", "js")),
      fun = function(m1, m2) m1 %*% m2,
      ncores = ncpar)$output1
dim(djb)
d2 <- dim(djb)
names(d2) <- c("ik", "is", "js", "jk")
dim(djb) <- d2
ik <- 1; is <- 2; js <- 3; jk <- 2
ik <- 2; is <- 3; js <- 2; jk <- 1
djb[ik, is, js, jk]



lxABprime_sums

