# Cargamos numero de array_id del SLURM
args <- commandArgs(trailing = TRUE)
slurm_array_id <- as.integer(args[1])

runs <- 1000000
#runif samples from a uniform distribution
xs <- runif(runs,min=-0.5,max=0.5)
ys <- runif(runs,min=-0.5,max=0.5)
in.circle <- xs^2 + ys^2 <= 0.5^2
mc.pi <- (sum(in.circle)/runs)*4

dfrm <- data.frame(array_id=slurm_array_id, apx_pi=mc.pi)
write.table(dfrm, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE, file="out")
