

# Univariate function to minimize
# t is the value of the selected coordinate
# X, Y are the data covariate and responses
# lambda is the penalization parameter
# to.min.index is the coordinate to minimize
# beta, beta0 are current values of the whole vector (all coordinates but to.min.index remain constant with the values of beta)
univ.function.to.min = function(t,loss.function, weights = NULL, pen.function,X,Y,lambda,to.min.index,beta, beta0, pen.aditive){
  beta[to.min.index] = t
  # print(loss.function)
  if(is.null(weights)){
      if(!pen.aditive){
      return(loss.function(X,Y,beta,beta0) + pen.function(beta = beta, lambda = lambda))
    } else{
      return(loss.function(X,Y,beta,beta0) + pen.function(beta = beta, j = to.min.index, lambda = lambda))
    }
  } else{
    if(!pen.aditive){
      # print("df")
      # print(loss.function(X,Y,beta,beta0,weights))
      # print(pen.function(beta = beta, lambda = lambda))
      return(loss.function(X,Y,beta,beta0,weights) + pen.function(beta = beta, lambda = lambda))
    } else{
      return(loss.function(X,Y,beta,beta0,weights) + pen.function(beta = beta, j = to.min.index, lambda = lambda))
    }
  }
}

# Function defined to minimize the scale of the whole beta vector. 
# c is the corresponding scale
scale.function.to.min = function(c,loss.function, weights = NULL, pen.function,X,Y,lambda,beta, beta0){
  beta = c*beta
  # print(beta)
  if(is.null(weights)){
    return(loss.function(X,Y,beta,beta0) + pen.function(beta = beta, lambda = lambda))
  } else{
    return(loss.function(X,Y,beta,beta0, weights) + pen.function(beta = beta, lambda = lambda))
  }
}

# The same as scale.funcion.to.min, without the penalization term
scale.function.to.min.nopen = function(c,loss.function, weights = NULL, X,Y,beta, beta0){
  beta = c*beta
  if(is.null(weights)){
    return(loss.function(X,Y,beta,beta0))
  } else{
    return(loss.function(X,Y,beta,beta0,weights))
  }
}

# Evaluates the function to minimize in a given beta and intercept
evaluate = function(loss.function, weights = NULL, pen.function,X,Y,lambda,beta, beta0){
  # if(is.null(beta0)){
  #   beta0 = 0
  # }
  if(is.null(weights)){
    return(scale.function.to.min(c = 1,loss.function = loss.function, pen.function = pen.function, X = X, Y = Y, lambda = lambda, beta = beta, beta0 = beta0))
  } else{
    return(scale.function.to.min(c = 1,loss.function = loss.function, weights = weights, pen.function = pen.function, X = X, Y = Y, lambda = lambda, beta = beta, beta0 = beta0))
  }
}



# Performs optimization on a single coordinate
update.one.coord = function(loss.function, weights = NULL, pen.function, X, Y, lambda, current.beta, beta0, index.to.update, pen.aditive){
  # print(current.beta)
  # print(index.to.update)
  # print(univ.function.to.min)
  # print(beta0)
  # print(lambda)
  # print(loss.function(X,Y,current.beta,beta0))
  # print(pen.function(beta = current.beta, j = index.to.update, lambda = lambda))
  # print(scad(beta = current.beta, j = index.to.update, lambda = lambda, alpha = 2.7))
  # print(univ.function.to.min(t = 2, loss.function = loss.function, pen.function = pen.function, X = X, Y = Y, lambda = lambda, to.min.index = index.to.update, beta = current.beta, beta0 = beta0, pen.aditive = pen.aditive))
  if(is.null(weights)){
    opt = optim(fn = univ.function.to.min, par = current.beta[index.to.update],loss.function = loss.function, pen.function = pen.function, X = X, Y = Y, lambda = lambda, to.min.index = index.to.update, beta = current.beta, beta0 = beta0, pen.aditive = pen.aditive)
    opt.value = opt$value
    opt.par = opt$par
    zero.value = univ.function.to.min(t = 0, loss.function = loss.function, pen.function = pen.function, X = X, Y = Y, lambda = lambda, beta = current.beta, beta0 = beta0, to.min.index = index.to.update, pen.aditive = pen.aditive)
    
    if(zero.value > opt.value){
      current.beta[index.to.update] = opt.par
    } else{
      current.beta[index.to.update] = 0
    }
    return(current.beta)
  } else{
    opt = optim(fn = univ.function.to.min, par = current.beta[index.to.update],loss.function = loss.function, weights = weights, pen.function = pen.function, X = X, Y = Y, lambda = lambda, to.min.index = index.to.update, beta = current.beta, beta0 = beta0, pen.aditive = pen.aditive)
    opt.value = opt$value
    opt.par = opt$par
    zero.value = univ.function.to.min(t = 0, loss.function = loss.function, weights = weights, pen.function = pen.function, X = X, Y = Y, lambda = lambda, beta = current.beta, beta0 = beta0, to.min.index = index.to.update, pen.aditive = pen.aditive)
    
    if(zero.value > opt.value){
      current.beta[index.to.update] = opt.par
    } else{
      current.beta[index.to.update] = 0
    }
    return(current.beta)
  }
}


# Optimizes the intercept
optim.intercept = function(loss.function, weights = NULL, X,Y,beta){
  if(is.null(weights)){
    fun = function(beta0){
      return(loss.function(X,Y,beta,beta0))
    }
    return(optim(par = 0, fn = fun)$par)
  } else{
    fun = function(beta0){
      return(loss.function(X,Y,beta,beta0,weights))
    }
    return(optim(par = 0, fn = fun)$par)
  }
}


# Performs onle cycle of minimizations. Includes:
# - Minimizing (in random order) over all covariates
# - Minimizing intercept
# - Minimizing scale of beta
one.cycle.descent = function(loss.function, weights = NULL, der.loss.function, pen.function,X,Y,lambda,beta.init,beta0.init, ssr.set,
                             intercept, pen.scale.invariant, pen.aditive){
  p = ncol(X)
  current.beta = beta.init
  beta0 = 0
  
  if(!is.null(beta0.init)){
    beta0 = beta0.init
  }
  
  order = sample(1:p)
  coord.old.value = evaluate(loss.function = loss.function, weights = weights, pen.function = pen.function, X = X, Y = Y, lambda = lambda, beta = current.beta, beta0 = beta0)
  for(j in 1:p){
    coord.to.update = order[j]
    if(coord.to.update %in% ssr.set || is.null(ssr.set)){
      current.beta = update.one.coord(loss.function = loss.function, weights = weights, pen.function = pen.function, X = X, Y = Y, lambda = lambda,
                                      current.beta = current.beta,beta0 = beta0, index.to.update = coord.to.update, pen.aditive)
    }
  }
  
  if(intercept){
    beta0 = optim.intercept(loss.function = loss.function, weights = weights, X = X, Y = Y, beta = current.beta)
  }
  opt.scale = optim.scale(loss.function = loss.function, weights = weights, der.loss.function = der.loss.function, pen.function = pen.function, lambda = lambda, X = X, Y = Y, beta = current.beta, beta0 = beta0, pen.scale.invariant = pen.scale.invariant)
  return(list(beta = opt.scale * current.beta, beta0 = beta0))
}


# Complete minimization method. Stops ciclying when improvement ratio is less than tol.
cyclical.descent = function(loss.function, der.loss.function, weights = NULL, pen.function, X,Y,lambda,beta.init,ssr.set = NULL, max.cycles = 100,
                            tol = 1e-3, beta0.init, intercept = T, pen.scale.invariant, pen.aditive){
  cycle.index = 1
  current.beta = beta.init
  if(is.null(ssr.set)){
    ssr.set = 1:ncol(X)
  }
  
  
  old.opt.value = evaluate(loss.function = loss.function, weights = weights, pen.function = pen.function, beta = current.beta, Y = Y, X = X, lambda = lambda, beta0 = beta0.init)
  stop = F
  while(!stop && cycle.index <= max.cycles){
    cycle.index = cycle.index + 1
    cycle.optim = one.cycle.descent(loss.function = loss.function, weights = weights, der.loss.function = der.loss.function, pen.function = pen.function,
                                    X = X, Y = Y, lambda = lambda, beta.init = current.beta, beta0.init = beta0.init, ssr.set = ssr.set,
                                    intercept = intercept, pen.scale.invariant = pen.scale.invariant, pen.aditive = pen.aditive)
    current.beta = cycle.optim$beta
    current.beta0 = cycle.optim$beta0
    
    opt.value = evaluate(loss.function = loss.function, weights = weights, pen.function = pen.function,beta = current.beta, beta0 = current.beta0, Y = Y, X = X, lambda = lambda)
    
    if((old.opt.value - opt.value)/ opt.value < tol){
      stop = T
    }
    old.opt.value = opt.value
  }
  
  return(list(beta = current.beta, beta0 = current.beta0, value = opt.value))
}

# Returns a complete list with fitted models for all lambdas and chooses one model with Cross-Validation.
cv = function(loss.function, der.loss.function, weights = NULL, pen.function, lambda.max.function, init.values.function, ssr.function, X,Y,nfolds = 10, lambda.len = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv = T, pen.scale.invariant, pen.aditive){
  
  # print(pen.function(rep(1,10),lambda = 0.01))
  if(is.null(lambdas)){
    lambda.max = lambda.max.function(X,Y, weights = weights)
    lambdas = 2**seq(log2(lambda.max), log2(lambda.max * lambdas.eps), length.out = lambda.len)
    # print(lambdas)
  }
  
  cv.matrix = NULL
  if(do.cv){
    n = nrow(X)
    p = ncol(X)
    folds = separate.folds(n, nfolds = nfolds)
    cv.matrix = matrix(NA, nfolds, length(lambdas))
    
    for(iSet in 1:nfolds){
      X.train = X[-folds[[iSet]],]
      Y.train = Y[-folds[[iSet]]]
      X.test = X[folds[[iSet]],]
      Y.test = Y[folds[[iSet]]]
      
      if(!is.null(weights)){
        weights.train = weights[-folds[[iSet]]] 
        weights.test = weights[folds[[iSet]]]
      }
      
      cv.init = init.values.function(X.train, Y.train)
      # print(cv.init)
      beta.init = cv.init$beta.init
      beta0.init = NULL
      
      if(intercept){
        beta0.init = cv.init$beta0.init
      }
      
      ssr.set = 1:p
      
      for(i in  1:length(lambdas)){
        
        if(is.null(weights)){
                    
          lambda.results = cyclical.descent(loss.function = loss.function, der.loss.function = der.loss.function, pen.function = pen.function,
                                            X = X.train, Y = Y.train, lambda = lambdas[i], ssr.set = ssr.set, beta.init = beta.init, beta0.init = beta0.init,
                                            intercept = intercept, pen.scale.invariant = pen.scale.invariant, pen.aditive = pen.aditive)
          lambda.beta = lambda.results$beta
          lambda.beta0 = lambda.results$beta0
          cv.matrix[iSet,i] = loss.function(X = X.test, Y = Y.test, beta = lambda.beta, beta0 = lambda.beta0)
          
          ssr.set = ssr.function(X = X.train, Y = Y.train, beta = lambda.beta, beta0 = lambda.beta0, lambda = lambdas[i])
          # print(c("Beta : ", lambda.beta))
        } else{
          lambda.results = cyclical.descent(loss.function = loss.function, der.loss.function = der.loss.function, weights = weights.train, pen.function = pen.function,
                                            X = X.train, Y = Y.train, lambda = lambdas[i], ssr.set = ssr.set, beta.init = beta.init, beta0.init = beta0.init,
                                            intercept = intercept, pen.scale.invariant = pen.scale.invariant, pen.aditive = pen.aditive)
          lambda.beta = lambda.results$beta
          lambda.beta0 = lambda.results$beta0
          cv.matrix[iSet,i] = loss.function(X = X.test, Y = Y.test, beta = lambda.beta, beta0 = lambda.beta0, weights = weights.test)
          ssr.set = ssr.function(X = X.train, Y = Y.train, beta = lambda.beta, beta0 = lambda.beta0, lambda = lambdas[i], weights = weights.train)
        }

        # print(c("SSR set : ", ssr.set))
      }
    }
  }
  
  
  return(get.fit.summary(loss.function, weights = weights, der.loss.function, pen.function, init.values.function, ssr.function, cv.matrix = cv.matrix,
                         X = X, Y = Y, lambdas = lambdas, intercept = intercept, pen.scale.invariant = pen.scale.invariant, pen.aditive = pen.aditive))
}



# Given a CV matrix, returns the chosen beta
get.fit.summary = function(loss.function,weights = NULL, der.loss.function, pen.function, init.values.function, ssr.function,
                           cv.matrix, X,Y,lambdas, intercept, pen.scale.invariant, pen.aditive){
  p = ncol(X)
  llambda = length(lambdas)
  min.index = 1
  if(!is.null(cv.matrix)){
      means = apply(cv.matrix,2,mean)
      min.index = min(which(means == min(means)))
  }
  

  initial = init.values.function(X = X, Y = Y)
  current.beta = initial$beta.init
  current.beta0 = NULL
  
  if(intercept){
    current.beta0 = initial$beta0.init
  }
  
  
  all.betas = matrix(NA, llambda, p)
  all.betas0 = rep(NA, llambda)
  all.loss.eval = rep(NA,llambda)
  all.df = rep(NA,llambda)
  
  
  ssr.set = 1:p
  for(i in 1: length(lambdas)){
    
    opt = cyclical.descent(loss.function = loss.function, weights = weights, der.loss.function = der.loss.function, pen.function = pen.function,
                           X = X, Y = Y, lambda = lambdas[i], ssr.set = ssr.set, beta.init = current.beta, beta0.init = current.beta0,
                           intercept = intercept, pen.scale.invariant = pen.scale.invariant, pen.aditive = pen.aditive)
    opt.beta = opt$beta
    opt.beta0 = opt$beta0
    
    if(is.null(weights)){
      ssr.set = ssr.function(X = X, Y = Y, beta = opt.beta, lambda = lambdas[i], beta0 = opt.beta0)
      
      all.betas[i,] = opt.beta
      all.betas0[i] = opt.beta0
      all.loss.eval[i] = loss.function(X = X, Y = Y, beta = opt.beta, beta0 = opt.beta0)
      all.df[i] = sum(opt.beta != 0) 
    } else{
      ssr.set = ssr.function(X = X, Y = Y, beta = opt.beta, lambda = lambdas[i], beta0 = opt.beta0, weights = weights)
      
      all.betas[i,] = opt.beta
      all.betas0[i] = opt.beta0
      all.loss.eval[i] = loss.function(X = X, Y = Y, beta = opt.beta, beta0 = opt.beta0, weights = weights)
      all.df[i] = sum(opt.beta != 0) 
    }
  }
  return(list(all.lambdas = lambdas, all.betas = all.betas, all.betas0 = all.betas0, all.df = all.df, all.loss.eval = all.loss.eval, beta = all.betas[min.index,], beta0 = all.betas0[min.index], lambda = lambdas[min.index]))
}


# Returns the scale for beta that minimizes the objetive function
optim.scale = function(loss.function, weights = NULL, der.loss.function, pen.function,lambda,X,Y,beta,beta0, pen.scale.invariant){
  if(is.null(weights)){
    if(pen.scale.invariant){
      scale.derivative = function(c){
        scores = c * X %*% beta + beta0
        der.phis = mapply(FUN = der.loss.function, score = scores, y = Y)
        return(mean(der.phis * scores))
      }
    return(optim(fn = scale.function.to.min.nopen, par = 1, gr = scale.derivative, loss.function = loss.function, X = X, Y = Y, beta = beta, beta0 = beta0)$par)
    } else{
      return(optim(fn = scale.function.to.min, par = 1, loss.function = loss.function, X = X, Y = Y, beta = beta, beta0 = beta0, pen.function = pen.function, lambda = lambda)$par)
    }
  } else{
    if(pen.scale.invariant){
      scale.derivative = function(c){
        scores = c * X %*% beta + beta0
        der.phis = mapply(FUN = der.loss.function, score = scores, y = Y, weights = weights)
        return(mean(der.phis * scores))
      }
      return(optim(fn = scale.function.to.min.nopen, par = 1, gr = scale.derivative, loss.function = loss.function, weights = weights, X = X, Y = Y, beta = beta, beta0 = beta0)$par)
    } else{
      return(optim(fn = scale.function.to.min, par = 1, loss.function = loss.function, weights = weights, X = X, Y = Y, beta = beta, beta0 = beta0, pen.function = pen.function, lambda = lambda)$par)
    }
  }
  
}



