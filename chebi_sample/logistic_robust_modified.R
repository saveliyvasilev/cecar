source("util.R")
library(foreach)
# library(doParallel)
# library(doSNOW)


univ.function.to.min.mod = function(t,X,Y,lambda,to.min.index,beta, beta0){
  beta[to.min.index] = t
  probs = predict.probs(beta, X,beta0 = beta0)
  return(mean((Y - probs)^2) + lambda * norm.mod(beta))
}

scale.function.to.min.mod = function(c,X,Y,lambda,to.min.index,beta, beta0){
  beta = c*beta
  probs = predict.probs(beta,X,beta0 = beta0)
  return(mean((Y - probs)^2) + lambda * norm.mod(beta))
}

univ.function.to.min.laso = function(t,X,Y,lambda,to.min.index,beta, beta0){
  beta[to.min.index] = t
  probs = predict.probs(beta,X, beta0 = beta0)
  return(mean((Y - probs)^2) + lambda * norm.quotient(beta))
}

scale.function.to.min.laso = function(c,X,Y,lambda,to.min.index,beta, beta0){
  
  beta = c * beta
  probs = predict.probs(beta,X, beta0 = beta0)
  return(mean((Y - probs)^2) + lambda * norm.quotient(beta))
}

update.one.coord = function(cyclical.function, X, Y, lambda, current.beta, beta0, index.to.update, tol = 1e-6){
  # print(current.beta)
  # print(index.to.update)
  # opt = optim(fn = cyclical.function, par = current.beta[index.to.update], X = X, Y = Y, lambda = lambda, to.min.index = index.to.update, beta = current.beta)
  # if(is.null(beta0)){
  #   beta0 = 0
  # }
  opt = try(nlm(f = cyclical.function, p = current.beta[index.to.update], X = X, Y = Y, lambda = lambda, to.min.index = index.to.update, beta = current.beta, beta0 = beta0))
  # print(opt)
  # opt.value = opt$value
  # opt.par = opt$par
  
  if(class(opt) == "try-error"){
    return(current.beta)
  }
  
  opt.value = opt$minimum
  opt.par = opt$estimate
  zero.value = cyclical.function(t = 0, X = X, Y = Y, lambda = lambda, beta = current.beta, beta0 = beta0, to.min.index = index.to.update)
  # print(zero.value)
  # print(opt.value)
  if((zero.value - opt.value) / opt.value > tol){
    # print(opt.par)
    current.beta[index.to.update] = opt.par
  } else{
    # print(c("Borro coordenada", index.to.update))
    current.beta[index.to.update] = 0
  }
  return(current.beta)
}

optim.intercept.cslaso = function(X,Y,beta, max.tries = 6, max = 20){
  
  # fun = function(beta0){
  #   scores = X %*% beta + beta0
  #   probs = inv.logit(scores)
  #   return(sum((Y - probs) * probs * (1 - probs)))
  # }
  
  fun = function(beta0){
    scores = X %*% beta  + beta0
    probs = inv.logit(scores)
    return(mean((Y - probs)**2))
  }
  
  return(optim(par = 0, fn = fun)$par)
  # print(fun(-max))
  # print(fun(max))
  # print(fun(1))
  
  # current.max = max
  # try = 1
  # beta0 = NA
  # while(try < max.tries && is.na(beta0)){
  #   beta0 = try(uniroot(f = fun, interval = c(-current.max,current.max))$root)
  #   if(class(beta0) == "try-error"){
  #     current.max = current.max / 2
  #     beta0 = NA
  #     try = try + 1
  #   }
  # }
  # return(beta0)
}

one.cycle.descent = function(method, X,Y,lambda,beta.init,beta0.init, ssr.set){
  p = ncol(X)
  current.beta = beta.init
  beta0 = 0
  if(!is.null(beta0.init)){
    beta0 = beta0.init
  }
  order = sample(1:p)
  
  # print(c("ssr set", ssr.set))
  if(identical(method, "CSMOD")){
    cyclical.function = univ.function.to.min.mod
    scale.function = scale.function.to.min.mod
  }
  if(identical(method, "CSLASO")){
    cyclical.function = univ.function.to.min.laso
    scale.function = scale.function.to.min.laso
  }
  coord.old.value = evaluate.cs(method = method, X = X, Y = Y, lambda = lambda, beta = current.beta, beta0 = beta0)
  for(j in 1:p){
    
    coord.to.update = order[j]
    # print(c("Coord to update ", coord.to.update))
    if(coord.to.update %in% ssr.set || is.null(ssr.set)){
      # print("entre!!")
      current.beta = update.one.coord(cyclical.function = cyclical.function, X = X, Y = Y, lambda = lambda, current.beta = current.beta,beta0 = beta0, index.to.update = coord.to.update)
    }
  }
  if(identical(method,"CSMOD")){
    opt.scale = optim(scale.function, par = 1, X = X, Y = Y, lambda = lambda, beta = current.beta, beta0 = beta0)$par
  }
  if(identical(method, "CSLASO")){
    opt.scale = try(optim.scale.border.check(X = X, Y = Y, direction = current.beta, beta0 = beta0))
    if(class(opt.scale) == "try-error"){
      opt.scale = 1
    }
    beta0 = try(optim.intercept.cslaso(X = X, Y = Y, beta = opt.scale * current.beta))
    if(is.null(beta0.init) || class(beta0) == "try-error"){
      beta0 = 0
    }
    # print(c("direction pre scale opt ",current.beta))
    # print(c("opt scale ",opt.scale))
  }
  return(list(beta = opt.scale * current.beta, beta0 = beta0))
}

evaluate.cs = function(method, X,Y,lambda,beta, beta0){
  if(is.null(beta0)){
    beta0 = 0
  }
  
  if(identical(method, "CSMOD")){
    return(scale.function.to.min.mod(c = 1, X = X, Y = Y, lambda = lambda, beta = beta, beta0 = beta0))
  }
  if(identical(method, "CSLASO")){
    return(scale.function.to.min.laso(c = 1, X = X, Y = Y, lambda = lambda, beta = beta, beta0 = beta0))
  }
}


cyclical.descent = function(method, X,Y,lambda,beta.init,ssr.set = NULL, max.cycles = 100, tol = 1e-3, beta0.init = NULL){
  cycle.index = 1
  current.beta = beta.init
  old.opt.value = evaluate.cs(method = method, beta = current.beta, Y = Y, X = X, lambda = lambda, beta0 = beta0.init)
  stop = F
  while(!stop && cycle.index <= max.cycles){
    # print(c("Cycle number : ", cycle.index))
    cycle.index = cycle.index + 1
    cycle.optim = one.cycle.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = current.beta, beta0.init = beta0.init, ssr.set = ssr.set)
    current.beta = cycle.optim$beta
    current.beta0 = cycle.optim$beta0
    # print(c("Current beta", current.beta))
    # print(c("Current beta0", current.beta0))
    opt.value = evaluate.cs(method = method, beta = current.beta, beta0 = current.beta0, Y = Y, X = X, lambda = lambda)
    
    if((old.opt.value - opt.value)/ opt.value < tol){
      stop = T
    }
    old.opt.value = opt.value
  }
  # print(c("Number of cycles : ", cycle.index))
  return(list(beta = current.beta, beta0 = current.beta0, value = opt.value))
}


#TODO : add method intercept
cv.cs.variations = function(method,X,Y,nfolds = 10 , lambda.len = 15, try = 1, max.tries = 5, lambdas = NULL, cl.size = NULL, lambdas.eps = 0.3, intercept = T){
  # TODO : adapt the lambda grid for this case
  # init = get.initial.values.cs(X,Y)
  # beta.init = init$beta.init
  # beta0.init = init$beta0.init
  # lambdas = try(get.bisection.grid.border.check(method = method,min = 0.001, max = 1, X =  X,Y =Y, beta.init = beta.init, len = lambda.len))
  # if(class(lambdas) == "try-error"){
  #   lambdas = get.lambda.grid.cs(X = X, Y = Y, alpha = 1, len = 15)
  # }
  lambda.max = get.lambda.max.2(X,Y)
  lambdas = 2**seq(log2(2*lambda.max), log2(lambda.max * lambdas.eps), length.out = lambda.len)
  n = nrow(X)
  p = ncol(X)
  folds = separate.folds(n, nfolds = nfolds)
  # print(folds)
  
  if(is.null(cl.size)){
    cv.matrix = matrix(NA, nfolds, length(lambdas))
  
    for(iSet in 1:nfolds){
      # print(c("fold : ", iSet))
      X.train = X[-folds[[iSet]],]
      Y.train = Y[-folds[[iSet]]]
      X.test = X[folds[[iSet]],]
      Y.test = Y[folds[[iSet]]]
      cv.init = get.initial.values.cs(X.train, Y.train)
      beta.init = cv.init$beta.init
      beta0.init = NULL
      
      if(intercept){
        beta0.init = cv.init$beta0.init
      }
      
      ssr.set = 1:p
      
      for(i in  1:lambda.len){
        # print(c("lambda ", lambdas[i]))
        # print(c("ssr set" , ssr.set))
        lambda.results = cyclical.descent(method = method, X = X.train, Y = Y.train, lambda = lambdas[i], ssr.set = ssr.set, beta.init = beta.init, beta0.init = beta0.init)
        lambda.beta = lambda.results$beta
        lambda.beta0 = lambda.results$beta0
        # print(c("ssr set" , ssr.set))
        # print(c("beta", lambda.beta))
        # print(c("beta 0",lambda.beta0))
        
        # active.set = which(cv.par.1 != 0)
        # cv.par.2 = cyclical.descent(X = X.train[,active.set], Y = Y.train, lambda = 0, beta.init = cv.par.1[active.set])$par
        cv.probs = predict.probs(beta = lambda.beta, beta0 = lambda.beta0, X = X.test)
        cv.matrix[iSet,i] = mean((Y.test - cv.probs)^2)
        # print(c("cv.par", cv.par))
        # print(c("active.set", active.set))
        # print(c("cv.par.2", cv.par.2))
        # print(c("Y.test",Y.test))
        # print(c("cv.probs", cv.probs))
        #cv.init = cv.par
        ssr.set = get.ssr.active.set(X = X.train, Y = Y.train, beta = lambda.beta, beta0 = lambda.beta0, lambda = lambdas[i])
        # print(c("ssr set" , ssr.set))
      }
    }
    
    
  } else {
    cl = makeCluster(cl.size)
    registerDoParallel(cl)
    
    tryCatch({
      #Cluster
      # cl = makeCluster(cl.size)
      # registerDoParallel(cl)
      # functions = c(lsf.str())
      # print(functions)
      
      
      
      # cv.matrix = matrix(NA, nfolds, length(lambdas))
      
      cv.matrix = foreach(iSet = 1:nfolds, .combine = 'rbind', .inorder = F, .verbose = T, .export = functions)%dopar%{
        
        ret.row = rep(NA, length(lambdas))
        
        print(c("fold : ", iSet))
        X.train = X[-folds[[iSet]],]
        Y.train = Y[-folds[[iSet]]]
        X.test = X[folds[[iSet]],]
        Y.test = Y[folds[[iSet]]]
        cv.init = get.initial.values.cs(X.train, Y.train)$beta.init
        
        for(i in  1:lambda.len){
          cv.par = cyclical.descent(method = method, X = X.train, Y = Y.train, lambda = lambdas[i], beta.init = cv.init)$par
          # active.set = which(cv.par.1 != 0)
          # cv.par.2 = cyclical.descent(X = X.train[,active.set], Y = Y.train, lambda = 0, beta.init = cv.par.1[active.set])$par
          cv.probs = predict.probs(beta = cv.par, X = X.test)
          ret.row[i] = mean((Y.test - cv.probs)^2)
          # print(c("lambda ", lambdas[i]))
          # print(c("cv.par", cv.par))
          # print(c("active.set", active.set))
          # print(c("cv.par.2", cv.par.2))
          # print(c("Y.test",Y.test))
          # print(c("cv.probs", cv.probs))
          #cv.init = cv.par
        }
        ret.row
      }
      stopCluster(cl)
      return(get.optimal.lambda.cs.mod(method = method, cv.matrix = cv.matrix, X = X, Y = Y, lambdas = lambdas))
    }, error = function(cond){
      stopCluster(cl)
      cl = makeCluster(cl.size)
      registerDoParallel(cl)
      print(c("Error: try number ", try))
      print(cond)
      if(try != max.tries){
        return(cv.cs.variations(method,X,Y,nfolds = 10 , lambda.len = 30, try = try + 1, max.tries = 5, lambdas = lambdas, cl.size = cl.size))
      }
      else{
        return(NA)
      }
    })
    
  }
  
  return(get.optimal.beta.cs.mod(method = method, cv.matrix = cv.matrix, X = X, Y = Y, lambdas = lambdas, intercept = intercept))
}



get.ssr.active.set = function(X,Y,beta,beta0,lambda,tol = 2){
  if(is.null(beta0)){
    beta0 = 0
  }
  
  p = ncol(X)
  norm2 = sqrt(sum(beta**2))
  norm1 = sum(abs(beta))
  importance = rep(NA, p)
  importance[beta!= 0] = Inf
  for(j in 1:p){
    if(beta[j] == 0){
      loss.der = get.loss.der(X = X, Y = Y, beta = beta, beta0 = beta0, coord.index = j)
      importance[j] = loss.der * norm2
    }
  }
  # print(lambda)
  # print(importance)
  return(which(abs(importance) > (lambda / tol)))
}



get.optimal.beta.cs.mod = function(method, cv.matrix, X,Y,lambdas, intercept){
  # print(cv.matrix)
  p = ncol(X)
  means = apply(cv.matrix,2,mean)
  # print(c("means : ", means))  
  min.index = min(which(means == min(means)))
  # print(c("min index : ", min.index))
  # print(c("chosen lambda : ", lambdas[min.index]))
  initial = get.initial.values.cs(X = X, Y = Y)
  lambdas.df = rep(NA,length(lambdas))
  lambdas.2norm = rep(NA,length(lambdas))
  norm.quotients = rep(NA, length(lambdas))
  current.beta = initial$beta.init
  current.beta0 = NULL
  if(intercept){
    current.beta0 = initial$beta0.init
  }
  
  ssr.set = 1:p
  for(i in 1: length(lambdas)){
    opt = cyclical.descent(method = method, X = X, Y = Y, lambda = lambdas[i], ssr.set = ssr.set, beta.init = current.beta, beta0.init = current.beta0)
    opt.beta = opt$beta
    opt.beta0 = opt$beta0
    # print(c("Lambda index: ", i))
    # print(c("Lambda value: ", lambdas[i]))
    # print(c("Optimal parameters: ", opt.beta))
    if(i == min.index){
      return(list(beta = opt.beta, beta0 = opt.beta0, lambda = lambdas[i]))
    }
    lambdas.df[i] = sum(opt.beta != 0) 
    lambdas.2norm[i] = sqrt(sum(opt.beta**2))
    norm.quotients[i] = norm.quotient(opt.beta)
    ssr.set = get.ssr.active.set(X = X, Y = Y, beta = opt.beta, lambda = lambdas[i], beta0 = opt.beta0)
    #current.beta = opt.beta
  }
  
}



# class.beta = function(beta, fixed.df.bound = 20,proportional.df.bound = 0.5, norm1.bound = 1000){
#   df = sum(beta != 0)
#   if(df <= 1){
#     return(2)
#   }
#   if((df > fixed.df.bound && df > proportional.df.bound * length(beta)) || sum(abs(beta))>norm1.bound){
#     return(0)
#   }
#   else{
#     return(1)
#   }
# }
# 
# get.class.for.lambda = function(method,X,Y,lambda, beta.init){
#   beta = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = beta.init)$par
#   return(class.beta(beta))
# }
# 
# get.bisection.grid.border.check = function(method, min = 0.001, max = 1, X, Y, beta.init, len = 30, iter = 1, max.iter = 5){
#   if(iter == max.iter){
#     return(2**seq(log2(0.2), log2(0.002), length.out = len))
#   }
#   min.class = get.class.for.lambda(method = method, X = X, Y = Y, lambda = min, beta.init = beta.init)
#   max.class = get.class.for.lambda(method = method, X = X, Y = Y, lambda = max, beta.init = beta.init)
#   print(c("Min : ", min))
#   print(c("Max : ", max))
#   print(c("Min class : ", min.class))
#   print(c("Max class : ", max.class))
#   if(min.class != 0){
#     return(get.bisection.grid.border.check(method, min = min/2, max = max, X = X, Y = Y, beta.init = beta.init, len = len, iter = iter + 1))
#   }
#   if(max.class != 2){
#     return(get.bisection.grid.border.check(method, min = min, max = max*2, X = X, Y = Y, beta.init = beta.init, len = len, iter = iter))
#   }
#   return(get.bisection.grid(method = method, min = min, max = max, X = X, Y = Y, beta.init = beta.init, len = len))
# }
# 
# get.bisection.grid = function(method, min , max , X, Y, beta.init, len, iter = 1, max.iter = 20){
#   if(iter == max.iter){
#     return(2**seq(log2(max), log2(min), length.out = len))
#   }
#   mid = 2**((log2(min) + log2(max))/2)
#   mid.class = get.class.for.lambda(method = method, X = X, Y = Y, lambda = mid, beta.init = beta.init)
#   print(c("Min : ", min))
#   print(c("Max : ", max))
#   print(c("Mid class : ", mid.class))
#   if(mid.class == 0){
#     return(get.bisection.grid(method = method, min = mid , max = max, X = X, Y = Y, beta.init = beta.init, len = len, iter = iter + 1))
#   }
#   if(mid.class == 2){
#     return(get.bisection.grid(method = method, min = min , max = mid, X = X, Y = Y, beta.init = beta.init, len = len, iter = iter + 1))
#   }
#   lambda.max = get.lambda.max(method = method, min = mid, max = max, X= X, Y = Y, beta.init = beta.init)
#   lambda.min = get.lambda.min(method = method, min = min, max = mid, X= X, Y = Y, beta.init = beta.init)
#   print(c("Lambda Min : ", lambda.min))
#   print(c("Lambda Max : ", lambda.max))
#   
#   return(2**seq(log2(lambda.max), log2(lambda.min), length.out = len))
# }
# 
# get.lambda.max = function(method, min,max,X,Y,beta.init, iter = 1, max.iter = 6){
#   mid = 2**((log2(min) + log2(max))/2)
#   
#   mid.class = get.class.for.lambda(method = method, X = X, Y = Y, lambda = mid, beta.init = beta.init)
#   if(iter == max.iter){
#     if(mid.class == 1){
#       return(max)
#     } else{
#       return(mid)
#     }
#   }
#   print(c("Min : ", min))
#   print(c("Max : ", max))
#   print(c("For get lambda max, mid class in iter", iter, " is ", mid.class))
#   if(mid.class <= 1){
#     return(get.lambda.max(method = method, min = mid ,max = max ,X = X,Y = Y,beta.init = beta.init, iter = iter + 1))
#   }
#   if(mid.class == 2){
#     return(get.lambda.max(method = method, min = min ,max = mid ,X = X,Y = Y,beta.init = beta.init, iter = iter + 1))
#   }
# }
# 
# get.lambda.min = function(method, min,max,X,Y,beta.init, iter = 1, max.iter = 6){
#   mid = 2**((log2(min) + log2(max))/2)
#   
#   mid.class = get.class.for.lambda(method = method, X = X, Y = Y, lambda = mid, beta.init = beta.init)
#   print(c("Min : ", min))
#   print(c("Max : ", max))
#   print(c("For get lambda min, mid class in iter", iter, " is ", mid.class))
#   if(iter == max.iter){
#     if(mid.class == 1){
#       return(mid)
#     } else{
#       return(max)
#     }
#   }
#   
#   if(mid.class == 1){
#     return(get.lambda.min(method = method, min = min ,max = mid ,X = X,Y = Y,beta.init = beta.init, iter = iter + 1))
#   }
#   if(mid.class == 0){
#     return(get.lambda.min(method = method, min = mid ,max = max ,X = X,Y = Y,beta.init = beta.init, iter = iter + 1))
#   }
# }

optim.scale.border.check = function(X, Y, direction, beta0, min= 0, max = 5, iter = 1, max.iter = 5){
  # print(c("Min : ", min))
  # print(c("Max : ", max))
  if(iter == max.iter){
    stop("Could not pass optim scale border check")
  }
  
  # #Debug
  # x = seq(0,5,0.01)
  # x.der = rep(NA, length(x))
  # for(k in 1:length(x)){
  #   x.der[k] = loss.derivative.for.scale(scale = x[k], X = X, Y = Y, direction = direction)
  # }
  # 
  # print(x.der)
  
  
  min.sign = sign(loss.derivative.for.scale(scale = min, X = X, Y = Y, direction = direction, beta0 = beta0))
  max.sign = sign(loss.derivative.for.scale(scale = max, X = X, Y = Y, direction = direction, beta0 = beta0))
  # print(c("Min sign: ", min.sign))
  # print(c("Max sign: ", max.sign))
  if(min.sign * max.sign == 1){
    return(optim.scale.border.check(X,Y,direction, min = min /2, max = max * 2, beta0 = beta0, iter = iter + 1))
  }
  return(optim.scale(X,Y,direction,beta0,min,max, min.sign, max.sign))
}

optim.scale = function(X,Y,direction,beta0,min,max, min.sign, max.sign, iter = 1, max.iter = 15){
  mid = (min + max)/2
  # print(c("Mid : ", mid))
  if(iter == max.iter){
    return(mid)
  }
  if(min.sign == 0){
    return(min)
  }
  if(max.sign == 0){
    return(max)
  }
  mid.sign = sign(loss.derivative.for.scale(scale = mid, X = X, Y = Y, direction = direction, beta0 = beta0))
  # print(c("Mid sign : ", mid.sign))
  if(mid.sign == 0){
    return(mid)
  }
  if(mid.sign == min.sign){
    return(optim.scale(X = X,Y = Y,direction = direction, beta0 = beta0, min = mid, max = max, min.sign = min.sign, max.sign = max.sign, iter = iter + 1))
  }
  if(mid.sign == max.sign){
    return(optim.scale(X = X,Y = Y,direction = direction, beta0 = beta0, min = min, max = mid, min.sign = min.sign, max.sign = max.sign, iter = iter + 1))
  }
}

loss.derivative.for.scale = function(scale,X,Y, direction, beta0){
  direction.scores = X %*% direction
  scores = scale* direction.scores + beta0
  probs = inv.logit(scores)
  return((-2/nrow(X))*sum( (Y - probs) * probs * (1-probs) * direction.scores))
}

kkt.check = function(X,Y,beta = NULL,lambda, seed = 1234){
  set.seed(seed)
  p = ncol(X)
  n = nrow(X)
  if(is.null(beta)){
    beta = get.beta.for.data(X = X, Y = Y, lambda = lambda)
  }
  # print(class(beta))
  # print(class(X))
  # print(beta)
  norm2 = sqrt(sum(beta**2))
  norm1 = sum(abs(beta))
  for(j in 1:p){
    loss.der = get.loss.der(X = X, Y = Y, beta = beta, coord.index = j)
    pen.der = (sign(beta[j]) * (norm2**2) - beta[j] * norm1) / (norm2**3)
    
    # print(c("Coord ", j, " loss.der ", loss.der))
    # print(c("Coord ", j, " pen.der ", pen.der))
    # print(c("Coord ", j, " loss.der/pen.der ", loss.der/pen.der))
    # 
    if(beta[j] != 0){
      print(c("Coordinate ", j, " not equal to zero"))
      print("The following two quantities should be similar")
      print(c("Loss derivative / pen der :" , loss.der / pen.der))
      print(c("Lambda :" , lambda))
      print(c("loss.der ", loss.der))
      print(c("pen.der ", pen.der))
    } else{
      print(c("Coordinate ", j, " equal to zero"))
      print("The first quantity should be smallet than the second")
      print(c("(abs) Norm2 * Loss derivative :" , abs(loss.der * norm2)))
      print(c("Lambda", lambda))
    }
  }
}

get.loss.der = function(X,Y,beta,beta0, coord.index){
  probs = predict.probs(beta = beta, X = X, beta0 = beta0)
  variances = probs * (1 - probs)
  norm2 = sqrt(sum(beta**2))
  norm1 = sum(abs(beta))
  return((2/nrow(X)) * sum(X[,coord.index] * (Y - probs) * variances))
}

get.lambda.max.2 = function(X,Y, max = 100, max.tries = 6){
  n = nrow(X)
  p = ncol(X)
  current.lambda = 0
  beta0 = log(mean(Y) / (1-mean(Y)))
  for(j in 1:p){
      
    fun = function(x){
    beta = rep(0,p)
    beta[j] = x
    return(get.loss.der(X,Y,beta,beta0 = beta0,j))
    }
    
    current.max = max
    try = 1
    beta.coord = NA
    
    while(try < max.tries && is.na(beta.coord)){
      beta.coord = try(uniroot(f = fun, interval = c(-current.max,current.max))$root)
      if(class(beta.coord) == "try-error"){
        current.max = current.max / 2
        beta.coord = NA
        try = try + 1
      }
    }
    

    
    if(!is.na(beta.coord)){
      beta = rep(0,p)
      beta[j] = beta.coord
      
      for(k in 1:p){
        if(k != j){
          lambda.candidate = abs(beta.coord) * get.loss.der(X = X, Y = Y, beta = beta, beta0 = beta0, coord.index = k) 
          current.lambda = max(current.lambda, lambda.candidate)
          print(c("j",j))
          print(c("k",k))
          print(c("current lambda :" , current.lambda))
        }
      }
    } else{
      warning(c("Could not find root for coordinate ", j))
    }
    
    
    
    
      
      
  }
  return(current.lambda)
}

ff = function(x){
  ret = rep(NA, length(x))
  for(i in 1:length(x)){
    beta = rep(0,80)
    beta[18] = x[i]
    ret[i] = get.loss.der(data$X, data$Y, beta, coord.index = 18)
  }
  return(ret)
}

