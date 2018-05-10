

# Returns the loss function for given data and beta
# The chosen loss function is same as in Croux and Haesbroeck
eval.loss.function.bypen = function(X,Y,beta,beta0){
  # ret = rep(NA,nrow(X))
  # for(i in 1:nrow(X)){
  #   score = beta0 + X[i,] %*% beta
  #   ret[i] = phi(score, Y[i])
  # }
  return(eval_loss_function_bypen_cpp(X,Y,beta,beta0))
}

eval.loss.function.wbypen = function(X,Y,beta,beta0, weights){
  # ret = rep(NA,nrow(X))
  # for(i in 1:nrow(X)){
  #   score = beta0 + X[i,] %*% beta
  #   ret[i] = phi(score, Y[i])
  # }
  return(eval_loss_function_wbypen_cpp(X,Y,beta,beta0, weights))
}

# Evaluates the loss derivative function, respect to one coordinate.
get.loss.der.for.coord.bypen = function(X,Y,beta,beta0,coord, weights = NULL){
  scores = X %*% beta + beta0
  if(is.null(weights)){
    return(mean(X[,coord] * mapply(der_phi_bypen_cpp, score = scores, y = Y)))
  }else{
    return(mean(X[,coord] * weights * mapply(der_phi_bypen_cpp, score = scores, y = Y)))
  }
}

# Gets the gradient of the loss function
get.loss.der.bypen = function(X,Y,beta,beta0){
  return(eval_der_loss_function_bypen_cpp(X,Y,beta,beta0))
}

# Gets the gradient of the loss function
get.loss.der.wbypen = function(X,Y,beta,beta0, weights){
  return(eval_der_loss_function_wbypen_cpp(X,Y,beta,beta0, weights))
}


eval.loss.function.pregpen = function(X,Y,beta,beta0){
  # ret = rep(NA,nrow(X))
  # for(i in 1:nrow(X)){
  #   score = beta0 + X[i,] %*% beta
  #   ret[i] = phi(score, Y[i])
  # }
  return(eval_loss_function_pregpen_cpp(X,Y,beta,beta0))
}



# Evaluates the loss derivative function, respect to one coordinate.
get.loss.der.for.coord.pregpen = function(X,Y,beta,beta0,coord){
  scores = X %*% beta + beta0
  return(mean(X[,coord] * mapply(der_phi_pregpen_cpp, score = scores, y = Y)))
}

# Gets the gradient of the loss function
get.loss.der.pregpen = function(X,Y,beta,beta0){
  return(eval_der_loss_function_pregpen_cpp(X,Y,beta,beta0))
}


# Perform an heuristic version of the Strong Safe Rules
get.ssr.active.set.bypen = function(X,Y,beta,beta0,lambda,tol = 1.2, weights = NULL){
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
      scores = X %*% beta + beta0
      if(is.null(weights)){
        loss.der = mean(X[,j] * mapply(der_phi_bypen_cpp, score = scores, y = Y))
      } else{
        loss.der = mean(X[,j] *  mapply(der_phi_bypen_cpp, score = scores, y = Y))
      }
      importance[j] = loss.der * norm2
    }
  }
  # print(lambda)
  # print(importance)
  return(which(abs(importance) > (lambda / tol)))
}

# Gets a lambda value that makes only one covariate active
get.lambda.max.bypen = function(X,Y, max = 100, max.tries = 6, intercept = T){
  n = nrow(X)
  p = ncol(X)
  current.lambda = 0
  beta0 = log(mean(Y) / (1-mean(Y)))
  if(!intercept){
    beta0 = 0
  }
  for(j in 1:p){
    fun = function(x){
      beta = rep(0,p)
      beta[j] = x
      return(get.loss.der.for.coord.bypen(X,Y,beta,beta0,j))
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
    
    # print(beta.coord)
    
    if(!is.na(beta.coord)){
      beta = rep(0,p)
      beta[j] = beta.coord
      
      for(k in 1:p){
        if(k != j){
          lambda.candidate = abs(beta.coord) * abs(get.loss.der.for.coord.bypen(X = X, Y = Y, beta = beta, beta0 = beta0, coord = k))
          # print(lambda.candidate)
          current.lambda = max(current.lambda, lambda.candidate)
          # print(c("j",j))
          # print(c("k",k))
          # print(c("current lambda :" , current.lambda))
        }
      }
    } else{
      # warning(c("Could not find root for coordinate ", j))
    }
    
  }
  print(2*current.lambda)
  return(2*current.lambda)
  # warning("Corrida excepcional de get.lambda.max!")
  # return(current.lambda)
}


get.lambda.max.bypen.2 = function(X,Y, max = 100, max.tries = 6, intercept = T, max.multip = 6, factor = 1.5){
  current.lambda = get.lambda.max.bypen(X,Y,max,max.tries,intercept)
  dfs = rep(NA,max.multip)
  init = get.robust.initial.2(X,Y)
  df.init = sum(init$beta.init != 0)
  for(i in 1:max.multip){
      df = get.df.for.lambda.bypen(X,Y, current.lambda, intercept = intercept, init = init)
      print(c("DF: ", df))
      if(df <= df.init){
        return(current.lambda)
      }
      dfs[i] = df
      current.lambda = factor * current.lambda
  }
  index = min(which(dfs == min(dfs)))
  return(lambda * (factor ** (index - 1)))
}

get.df.for.lambda.bypen = function(X,Y,lambda, intercept = T, init = NULL){
  
  if(is.null(init)){
    init = get.robust.initial.2(X = X, Y = Y)
  }
  
  fit = cyclical.descent(loss.function = eval.loss.function.bypen, der.loss.function = get.loss.der.bypen, pen.function = norm.quotient,
                   X = X, Y = Y, lambda = lambda, beta.init = init$beta.init,
                   beta0.init = init$beta0.init,ssr.set = NULL, pen.scale.invariant = T,
                   pen.aditive = F, intercept = intercept)
  return(sum(fit$beta!= 0))
}


# Performs cross-validation
cv.bypen = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv){
  return(cv(loss.function = eval.loss.function.bypen, der.loss.function = get.loss.der.bypen,
            pen.function = norm.quotient, lambda.max.function = get.lambda.max.lasso,
            init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.bypen,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps,
            intercept = intercept, do.cv = do.cv, pen.scale.invariant = T, pen.aditive = F))
}

# Performs cross-validation
cv.hubwbypen = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv){
  weights = get.weights(X, method = "HUBER", constant = 1.345)
  return(cv(loss.function = eval.loss.function.wbypen, der.loss.function = get.loss.der.wbypen, weights = weights,
            pen.function = norm.quotient, lambda.max.function = get.lambda.max.lasso,
            init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.bypen,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps,
            intercept = intercept, do.cv = do.cv, pen.scale.invariant = T, pen.aditive = F))
}

# Performs cross-validation
cv.bisqwbypen = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv){
  weights = get.weights(X, method = "BISQUARE", constant = 4.685)
  
  
  return(cv(loss.function = eval.loss.function.wbypen, der.loss.function = get.loss.der.wbypen,  weights = weights,
            pen.function = norm.quotient, lambda.max.function = get.lambda.max.lasso,
            init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.bypen,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps,
            intercept = intercept, do.cv = do.cv, pen.scale.invariant = T, pen.aditive = F))
}

# Performs cross-validation
cv.hardwbypen = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv){
  # weights = get.weights(X, method = "HARD", constant = sqrt(qchisq(0.975, df = ncol(X))))
  weights = get.weights.box(X)
  
  return(cv(loss.function = eval.loss.function.wbypen, der.loss.function = get.loss.der.wbypen, weights = weights,
            pen.function = norm.quotient, lambda.max.function = get.lambda.max.lasso,
            init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.bypen,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps,
            intercept = intercept, do.cv = do.cv, pen.scale.invariant = T, pen.aditive = F))
}

cv.hubwbyscad = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  weights = get.weights(X, method = "HUBER", constant = 1.345)
  
  scad.fixed.alpha = function(beta, j = NULL, lambda){
    # print(alpha)
    return(scad(beta, alpha, j, lambda))
  }
  # print(scad.fixed.alpha(c(1,1,1,0,0,0), j = 2, lambda = 0.01))
  return(cv(loss.function = eval.loss.function.wbypen, der.loss.function = get.loss.der.wbypen, weights = weights, pen.function = scad.fixed.alpha,
            lambda.max.function = get.lambda.max.lasso, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.lasso,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}

# Performs cross-validation
cv.hubwbymcp = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  weights = get.weights(X, method = "HUBER", constant = 1.345)
  
  mcp.fixed.alpha = function(beta, j = NULL, lambda){
    return(mcp(beta, alpha, j, lambda))
  }
  return(cv(loss.function = eval.loss.function.wbypen, der.loss.function = get.loss.der.wbypen, weights = weights, pen.function = mcp.fixed.alpha,
            lambda.max.function = get.lambda.max.lasso, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.lasso,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}

cv.bisqwbyscad = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  
  weights = get.weights(X, method = "BISQUARE", constant = 4.685)
  
  scad.fixed.alpha = function(beta, j = NULL, lambda){
    # print(alpha)
    return(scad(beta, alpha, j, lambda))
  }
  # print(scad.fixed.alpha(c(1,1,1,0,0,0), j = 2, lambda = 0.01))
  return(cv(loss.function = eval.loss.function.wbypen, der.loss.function = get.loss.der.wbypen, weights = weights, pen.function = scad.fixed.alpha,
            lambda.max.function = get.lambda.max.lasso, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.lasso,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}

# Performs cross-validation
cv.bisqwbymcp = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  weights = get.weights(X, method = "BISQUARE", constant = 4.685)
  print(weights)
  mcp.fixed.alpha = function(beta, j = NULL, lambda){
    return(mcp(beta, alpha, j, lambda))
  }
  return(cv(loss.function = eval.loss.function.wbypen, der.loss.function = get.loss.der.wbypen, weights = weights, pen.function = mcp.fixed.alpha,
            lambda.max.function = get.lambda.max.lasso, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.lasso,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}

cv.hardwbyscad = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  weights = get.weights.box(X)
  
  scad.fixed.alpha = function(beta, j = NULL, lambda){
    # print(alpha)
    return(scad(beta, alpha, j, lambda))
  }
  # print(scad.fixed.alpha(c(1,1,1,0,0,0), j = 2, lambda = 0.01))
  return(cv(loss.function = eval.loss.function.wbypen, der.loss.function = get.loss.der.wbypen, weights = weights, pen.function = scad.fixed.alpha,
            lambda.max.function = get.lambda.max.lasso, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.lasso,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}

# Performs cross-validation
cv.hardwbymcp = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  # weights = get.weights(X, method = "HARD", constant = sqrt(qchisq(0.975, df = ncol(X))))
  
  weights = get.weights.box(X)
  
  mcp.fixed.alpha = function(beta, j = NULL, lambda){
    return(mcp(beta, alpha, j, lambda))
  }
  return(cv(loss.function = eval.loss.function.wbypen, der.loss.function = get.loss.der.wbypen, weights = weights, pen.function = mcp.fixed.alpha,
            lambda.max.function = get.lambda.max.lasso, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.lasso,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}

# Performs cross-validation
cv.pregpen = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv){
  return(cv(loss.function = eval.loss.function.pregpen, der.loss.function = get.loss.der.pregpen,
            pen.function = norm.quotient, lambda.max.function = get.lambda.max.lasso,
            init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.bypen,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps,
            intercept = intercept, do.cv = do.cv, pen.scale.invariant = T, pen.aditive = F))
}

# Performs cross-validation
cv.bypenmod = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, constant = 1){
  norm.quotient.mod.fixed.constant = function(beta, lambda){
    return(norm.quotient.mod(beta = beta, lambda = lambda, constant = constant))
  }
  return(cv(loss.function = eval.loss.function.bypen, der.loss.function = get.loss.der.bypen,
            pen.function = norm.quotient.mod.fixed.constant, lambda.max.function = get.lambda.max.lasso,
            init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.bypen,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps,
            intercept = intercept, do.cv = do.cv, pen.scale.invariant = F, pen.aditive = F))
}


#####################
# Clarke conditions
#####################




# Returns penalty derivative for non-zero coordinates. 
get.penalty.der.bypen = function(beta){
  signs = sign(beta)
  norm1 = sum(abs(beta))
  norm2 = sqrt(sum(beta^2))
  ret = (signs * norm2^2 - beta * norm1) / norm2^3
  # print(ret)
  ret[signs == 0] = NA
  return(ret)
}


# Verifies Clarke conditions for local minimum
verify.clarke = function(X,Y,beta.est,beta0.est, lambda){
  loss.der = get.loss.der.bypen(X,Y,beta.est,beta0.est)
  pen.der = get.penalty.der.bypen(beta.est)
  norm2 = sqrt(sum(beta.est^2))
  ####Must be zero for non-zero components of beta.est
  print(loss.der + lambda * pen.der)
  
  ####Zero components must be less than 1 in absolute value
  print(norm2 * loss.der / lambda) 
}


