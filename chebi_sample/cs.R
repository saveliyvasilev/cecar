# Returns the loss function for given data and beta
# The chosen loss function is same as in Croux and Haesbroeck
eval.loss.function.cs = function(X,Y,beta,beta0){
  return(eval_loss_function_cs_cpp(X,Y,beta,beta0))
}

# Gets the gradient of the loss function
get.loss.der.cs = function(X,Y,beta,beta0){
  return(eval_der_loss_function_cs_cpp(X,Y,beta,beta0))
}

get.lambda.max.cs = function(X,Y,alpha = 1, weights = NULL){
  y.mean = mean(Y)
  factor = max(abs(t(X) %*% Y))
  lambda.max = y.mean * (1 - y.mean) * factor *2 /(nrow(X) * alpha)
  return(lambda.max)
}

# Perform an heuristic version of the Strong Safe Rules
get.ssr.active.set.cs = function(X,Y,beta,beta0,lambda,tol = 1.2){
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
      loss.der = mean(X[,j] * mapply(der_phi_cs_cpp, score = scores, y = Y))
      importance[j] = loss.der * norm2
    }
  }
  # print(lambda)
  # print(importance)
  return(which(abs(importance) > (lambda / tol)))
}

cv.cs = function(X,Y,nfolds = 10, lambda.len = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T){
  return(cv(loss.function = eval.loss.function.cs, der.loss.function = get.loss.der.cs, pen.function = norm.quotient, lambda.max.function = get.lambda.max.cs, init.values.function = get.robust.initial.2 , ssr.function = get.ssr.active.set.cs,
            X = X, Y = Y, nfolds = nfolds, lambda.len = lambda.len, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept, do.cv = T, pen.scale.invariant = T, pen.aditive = F))
}

