# Performs cross-validation
cv.bynet = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  elastic.net.fixed.alpha = function(beta, j = NULL, lambda){
    return(elastic.net(beta, alpha, j, lambda))
  }
  get.lambda.max.bynet.fixed.alpha = function(X,Y){
    return(get.lambda.max.bynet(X,Y,alpha))
  }
  get.ssr.active.set.bynet.fixed.alpha = function(X,Y,beta,beta0,lambda,tol = 2){
    return(get.ssr.active.set.bynet(X,Y,beta,beta0,lambda,alpha,tol = 2))
  }
  return(cv(loss.function = eval.loss.function.bypen, der.loss.function = get.loss.der.bypen, pen.function = elastic.net.fixed.alpha,
            lambda.max.function = get.lambda.max.bynet.fixed.alpha, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.bynet.fixed.alpha,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}


get.lambda.max.bynet = function(X,Y, alpha){
  ret = 0
  for(i in 1: ncol(X)){
    ret = max(ret, abs(get.loss.der.for.coord.bypen(X = X, Y = Y, beta = rep(0, ncol(X)), beta0 = 0, coord = i))/alpha)
  }
  return(ret)
}

# Perform an heuristic version of the Strong Safe Rules
get.ssr.active.set.bynet = function(X,Y,beta,beta0,lambda,alpha,tol = 2){
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
      loss.der = mean(X[,j] * mapply(der_phi_bypen_cpp, score = scores, y = Y))
      importance[j] = loss.der / alpha
    }
  }
  # print(lambda)
  # print(importance)
  return(which(abs(importance) > (lambda / tol)))
}

