# Performs cross-validation
cv.byscad = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  scad.fixed.alpha = function(beta, j = NULL, lambda){
    # print(alpha)
    return(scad(beta, alpha, j, lambda))
  }
  # print(scad.fixed.alpha(c(1,1,1,0,0,0), j = 2, lambda = 0.01))
  return(cv(loss.function = eval.loss.function.bypen, der.loss.function = get.loss.der.bypen, pen.function = scad.fixed.alpha,
            lambda.max.function = get.lambda.max.lasso, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.lasso,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}

# Performs cross-validation
cv.bymcp = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  mcp.fixed.alpha = function(beta, j = NULL, lambda){
    return(mcp(beta, alpha, j, lambda))
  }
  return(cv(loss.function = eval.loss.function.bypen, der.loss.function = get.loss.der.bypen, pen.function = mcp.fixed.alpha,
            lambda.max.function = get.lambda.max.lasso, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.lasso,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}

# Performs cross-validation
cv.pregscad = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  scad.fixed.alpha = function(beta, j = NULL, lambda){
    # print(alpha)
    return(scad(beta, alpha, j, lambda))
  }

  # print(scad.fixed.alpha(c(1,1,1,0,0,0), j = 2, lambda = 0.01))
  return(cv(loss.function = eval.loss.function.pregpen, der.loss.function = get.loss.der.pregpen, pen.function = scad.fixed.alpha,
            lambda.max.function = get.lambda.max.lasso, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.lasso,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}
# Performs cross-validation
cv.pregmcp = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  mcp.fixed.alpha = function(beta, j = NULL, lambda){
    return(mcp(beta, alpha, j, lambda))
  }
  return(cv(loss.function = eval.loss.function.pregpen, der.loss.function = get.loss.der.pregpen, pen.function = mcp.fixed.alpha,
            lambda.max.function = get.lambda.max.lasso, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.lasso,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}


get.lambda.max.lasso = function(X,Y, weights = NULL){
  ret = 0
  for(i in 1: ncol(X)){
    ret = max(ret, abs(get.loss.der.for.coord.bypen(X = X, Y = Y, beta = rep(0, ncol(X)), beta0 = 0, coord = i, weights = weights)))
  }
  print(ret)
  return(ret)
}

# Perform an heuristic version of the Strong Safe Rules
get.ssr.active.set.lasso = function(X,Y,beta,beta0,lambda,tol = 2, weights = NULL){
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
        loss.der = mean(X[,j] * weights * mapply(der_phi_bypen_cpp, score = scores, y = Y))
      }
      importance[j] = loss.der
    }
  }
  # print(c("Lambda: ", lambda))
  # print(c("Importance: ", importance))
  return(which(abs(importance) > (lambda / tol)))
}