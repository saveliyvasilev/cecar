

eval.loss.function.basu = function(X,Y,beta,beta0){
  # ret = rep(NA,nrow(X))
  # for(i in 1:nrow(X)){
  #   score = beta0 + X[i,] %*% beta
  #   ret[i] = phi(score, Y[i])
  # }
  return(eval_loss_function_basu_cpp(X,Y,beta,beta0))
}


# Evaluates the loss derivative function, respect to one coordinate.
get.loss.der.for.coord.basu = function(X,Y,beta,beta0,coord){
  scores = X %*% beta + beta0
  return(mean(X[,coord] * mapply(der_phi_basu_cpp, score = scores, y = Y)))
}

# Gets the gradient of the loss function
get.loss.der.basu = function(X,Y,beta,beta0){
  return(eval_der_loss_function_basu_cpp(X,Y,beta,beta0))
}


# Perform an heuristic version of the Strong Safe Rules
get.ssr.active.set.basu = function(X,Y,beta,beta0,lambda,tol = 1.2){
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
      loss.der = mean(X[,j] * mapply(der_phi_basu_cpp, score = scores, y = Y))
      importance[j] = loss.der * norm2
    }
  }
  # print(lambda)
  # print(importance)
  return(which(abs(importance) > (lambda / tol)))
}

# Gets a lambda value that makes only one covariate active
get.lambda.max.basusp = function(X,Y, max = 100, max.tries = 6, intercept = T){
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
      return(get.loss.der.for.coord.basu(X,Y,beta,beta0,j))
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
          lambda.candidate = abs(beta.coord) * abs(get.loss.der.for.coord.basu(X = X, Y = Y, beta = beta, beta0 = beta0, coord = k))
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


get.lambda.max.basusp.2 = function(X,Y, max = 100, max.tries = 6, intercept = T, max.multip = 6, factor = 1.5){
  init.lambda = get.lambda.max.basusp(X,Y,max,max.tries,intercept)
  current.lambda = init.lambda
  dfs = rep(NA,max.multip)
  init = get.robust.initial.2(X,Y)
  df.init = sum(init$beta.init != 0)
  for(i in 1:max.multip){
    df = get.df.for.lambda.basusp(X,Y, current.lambda, intercept = intercept, init = init)
    print(c("DF: ", df))
    if(df <= df.init){
      return(current.lambda)
    }
    dfs[i] = df
    current.lambda = factor * current.lambda
  }
  index = min(which(dfs == min(dfs)))
  return(init.lambda * (factor ** (index - 1)))
}

get.df.for.lambda.basusp = function(X,Y,lambda, intercept = T, init = NULL){

  if(is.null(init)){
    init = get.robust.initial.2(X = X, Y = Y)
  }

  fit = cyclical.descent(loss.function = eval.loss.function.basu, der.loss.function = get.loss.der.basu, pen.function = norm.quotient,
                         X = X, Y = Y, lambda = lambda, beta.init = init$beta.init,
                         beta0.init = init$beta0.init,ssr.set = NULL, pen.scale.invariant = T,
                         pen.aditive = F, intercept = intercept)
  return(sum(fit$beta!= 0))
}


# Performs cross-validation
cv.basusp = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv){
  return(cv(loss.function = eval.loss.function.basu, der.loss.function = get.loss.der.basu,
            pen.function = norm.quotient, lambda.max.function = get.lambda.max.lasso,
            init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.basu,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps,
            intercept = intercept, do.cv = do.cv, pen.scale.invariant = T, pen.aditive = F))
}

# Performs cross-validation
cv.basuscad = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  # print(alpha)
  scad.fixed.alpha = function(beta, j = NULL, lambda){
    return(scad(beta, alpha, j, lambda))
  }
  # print(scad.fixed.alpha(c(1,1,1,0,0,0), j = 2, lambda = 0.01))
  return(cv(loss.function = eval.loss.function.basu, der.loss.function = get.loss.der.basu, pen.function = scad.fixed.alpha,
            lambda.max.function = get.lambda.max.lasso, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.lasso,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}

# Performs cross-validation
cv.basumcp = function(X,Y,nfolds = 10, nlambda = 12, lambdas = NULL, lambdas.eps = 0.2, intercept = T, do.cv, alpha = 1){
  mcp.fixed.alpha = function(beta, j = NULL, lambda){
    return(mcp(beta, alpha, j, lambda))
  }
  return(cv(loss.function = eval.loss.function.basu, der.loss.function = get.loss.der.basu, pen.function = mcp.fixed.alpha,
            lambda.max.function = get.lambda.max.lasso, init.values.function = get.robust.initial.2, ssr.function = get.ssr.active.set.lasso,
            X = X, Y = Y, nfolds = nfolds, lambda.len = nlambda, lambdas = lambdas, lambdas.eps = lambdas.eps, intercept = intercept,
            do.cv = do.cv, pen.scale.invariant = F, pen.aditive = T))
}
