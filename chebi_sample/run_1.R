source("main.R")


ns = c(50,150,300)
ps = c(5,40,80)
niter = 10
start.index = 1
alpha = 2.7

for(i in 1: length(ns)){
  for(j in 1:length(ps)){
    simul.many(method = "GLMNET", n = ns[i], p = ps[j], niter = niter, start.index = start.index, alpha = 1)
    # simul.many(method = "HARDWBYPEN", n = ns[i], p = ps[j], niter = niter, start.index = start.index)
    # simul.many(method = "HARDWBYMCP", n = ns[i], p = ps[j], niter = niter, start.index = start.index, alpha = 2.7)
    # simul.many(method = "HARDWBYPEN", n = ns[i], p = ps[j], niter = niter, start.index = start.index, eps = 0.05, out.sigma.code = 20)
    # simul.many(method = "HARDWBYMCP", n = ns[i], p = ps[j], niter = niter, start.index = start.index, alpha = 2.7, eps = 0.05, out.sigma.code = 20)
    # simul.many(method = "HARDWBYPEN", n = ns[i], p = ps[j], niter = niter, start.index = start.index, eps = 0.1, out.sigma.code = 20)
    # simul.many(method = "HARDWBYMCP", n = ns[i], p = ps[j], niter = niter, start.index = start.index, alpha = 2.7, eps = 0.1, out.sigma.code = 20)
  }
}

