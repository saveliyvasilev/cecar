library("parallel")

a <- function(s) { return (2*s) }
mclapply(c(1:10^6), a, mc.cores=9)
