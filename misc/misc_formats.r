sprintf("%f", pi)
sprintf("%.3f", pi)
sprintf("%1.0f", pi)
sprintf("%5.1f", pi)
sprintf("%05.1f", pi)
sprintf("%+f", pi)
sprintf("% f", pi)
sprintf("%-10f", pi) # left justified
sprintf("%e", pi)
sprintf("%E", pi)
sprintf("%g", pi)
sprintf("%g",   1e6 * pi) # -> exponential
sprintf("%.9g", 1e6 * pi) # -> "fixed"
sprintf("%G", 1e-6 * pi)


xx <- pi * 10^(-5:4)
cbind(format(xx, digits = 4), formatC(xx))
cbind(formatC(xx, width = 9, flag = "-"))
cbind(formatC(xx, digits = 5, width = 8, format = "f", flag = "0"))
cbind(format(xx, digits = 4), formatC(xx, digits = 4, format = "fg"))
r <- c("76491283764.97430", "29.12345678901", "-7.1234", "-100.1","1123")
## American:
prettyNum(r, big.mark = ",")
prettyNum(r, big.mark = ",")

cbind(ff <- format(1.2345 + 10^(0:5), width = 11, big.mark = "'"))
## all with same width (one more than the specified minimum)

cbind(ff <- format(1.23456789 + 10^(0:5), width = 6, big.mark = ",", mode="integer", digits=0))

(x <- 1.23456789*10^c(-1, 0, 1:8))
format(x, digits = 2, nsmall = 1)
prettyNum(x, big.mark = ",")
formatC(x, big.mark = ",")
formatC(x, big.mark = ",", digits=0)
formatC(x, big.mark = ",", width=8, digits=3)
formatC(as.integer(x), big.mark = ",", width=8, digits=3)
formatC(as.integer(x), big.mark = ",")
formatC(as.integer(x), big.mark = ",", width=6)
formatC(as.integer(x), big.mark = ",", width=0, digits=6)
formatC(as.integer(x), big.mark = ",", digits=0, width=6, format="g", flag="#")

labvec <- function(x){
  xmin <- min(x)
  xmax <- max(x)
  
}


