# Returns the loss function for given data and beta
# The chosen loss function is same as in Croux and Haesbroeck
eval.loss.function.dev = function(X,Y,beta,beta0){
  return(eval_loss_function_dev_cpp(X,Y,beta,beta0))
}

# Gets the gradient of the loss function
get.loss.der.dev = function(X,Y,beta,beta0){
  return(eval_der_loss_function_dev_cpp(X,Y,beta,beta0))
}


## TODO: cambiar esto!!
get.lambda.max.dev = function(X,Y,alpha = 1, len = 50){
  return(max(glmnet(X,Y, alpha = alpha)$lambda))
}

# Perform an heuristic version of the Strong Safe Rules
get.ssr.active.set.dev = function(X,Y,beta,beta0,lambda,tol = 1.2){
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
      loss.der = mean(X[,j] * mapply(der_dev_cpp, score = scores, y = Y))
      importance[j] = loss.der * norm2
    }
  }
  # print(lambda)
  # print(importance)
  return(which(abs(importance) > (lambda / tol)))
}

cv.devpen = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T){
  return(cv(loss.function = eval.loss.function.dev, der.loss.function = get.loss.der.dev, pen.function = norm.quotient,
            lambda.max.function = get.lambda.max.dev, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.dev,
            X = X, Y = Y, nfolds = nfolds, lambdas = lambdas, lambda.len = nlambda,
            lambdas.eps = lambdas.eps, intercept = intercept, do.cv = T, pen.scale.invariant = T, pen.aditive = F))
}


init.values.flags = function(X,Y){
  return(list(beta.init = res.bypen$beta, beta0.init = res.bypen$beta0))
}