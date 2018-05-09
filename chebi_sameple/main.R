library(mvtnorm)
library(glmnet)
# library(inline)
library(Rcpp)
library(RcppArmadillo)
library(robustbase)
library(rrcov)
library(rrcovHD)
library(cvplogistic)
library(Gmedian)
library(glasso)
#library(caret)
source("croux_by.R")
# source("rLogistic.R") source("robust_CV_cs.R")
# source("logistic_robust_modified.R")
source("bypen.R")
source("bynet.R")
source("cs.R")
source("dev.R")
source("byscad_bymcp.R")
source("basu.R")
source("general_sparse.R")
sourceCpp("C++/util_bypen.cpp")
sourceCpp("C++/util_cs.cpp")
sourceCpp("C++/util_dev.cpp")
sourceCpp("C++/util_basu.cpp")
sourceCpp("C++/penalties.cpp")




# source("noverlap.R")


## Function to obtain the design matrix X. Outliers can be generated with one fixed value for X and Y
## or randomly with specified mean and covariance parameters.
##
## Update 3/4/17 : Outilers are added and not replaced in data
##
## n total sample size
## p dimension size
## mu distribution mean
## sigma distribution covariance matrix
## nonzero.beta beta value for nonzero components (the remaining components are assumed to be zero)
## eps outlier proportion
## intercept value for the intercept
## out.X X value for all outliers
## out.Y Y value for all outliers
## out.mu outiler distribution mean
## out.sigma outlier distribution covariance matrix
get.data.norm = function(n,p,mu = NULL,sigma = NULL,nonzero.beta,eps = 0, intercept = NULL, out.X = NULL, out.Y = NULL, out.mu = NULL, out.sigma = NULL){
  if(is.null(mu)){
    mu = rep(0,p)
  }
  if(is.null(sigma)){
    sigma = diag(p)
  }
  if(eps == 0){
    out.sigma = NULL
  }
  beta = c(nonzero.beta, rep(0,p - length(nonzero.beta)))
  out.numb = floor(eps * n)
  X = rmvnorm(n = n, mean = mu, sigma = sigma)
  if(!is.null(intercept)){
    X = cbind(rep(1,n),X)
    beta = c(intercept, beta)
  }
  X.lc = X %*% beta
  probs = inv.logit(X.lc)
  Y = rbinom(n = n, size = 1, prob = probs)
  
  X.out = NULL
  if(!is.null(out.X)){
    X.out = matrix(data = rep(out.X,out.numb), nrow = out.numb, ncol = p, byrow = T)
    if(!is.null(intercept)){
      X.out = cbind(rep(1,out.numb),X.out)
    }
    X = rbind(X, X.out)
    Y = c(Y, rep(out.Y, out.numb))
    
  }
  
  if(!is.null(out.mu) || !is.null(out.sigma)){
    if(is.null(out.mu)){
      out.mu = rep(0,p)
    }
    if(is.null(out.sigma)){
      out.sigma = diag(p)
    }
    X.out = rmvnorm(out.numb, out.mu, out.sigma)
    if(!is.null(intercept)){
      X.out = cbind(rep(1,out.numb),X.out)
    }
    X = rbind(X, X.out)
    Y.out = apply(X.out,1,get.out.Y.for.X, beta = beta)
    Y = c(Y, Y.out)
  }
  
  if(!is.null(X.out)){
    X.out.lc = X.out %*% beta
    probs = c(probs, inv.logit(X.out.lc))
  }
  
  if(!is.null(intercept)){
    X = X[,-1]
  }
  return(list(X = X, Y = Y, probs = probs))
}


## Gets the outlier corresponding to a given X data and beta
get.out.Y.for.X = function(x,beta){
  lc = sum(x*beta)
  return(1*(inv.logit(lc) < 0.5))
}

###################################
## GLM Method
###################################

# Runs the GLM Method
# It estimates the intercept only if intercept=T. It is set to zero otherwise. 
#
# Returns the corresponding beta and intercept
run.glm = function(X,Y, intercept = T){
  if(intercept){
    coefs = coef(glm(Y~X, family = "binomial"))
    return(list(beta = coefs[-1], beta0 = coefs[1], lambda = 0))
  } else{
    coefs = coef(glm(Y~X-1, family = "binomial"))
    return(list(beta = coefs, beta0 = 0, lambda = 0))
  }
}


###################################
## GLMNET Method
###################################

# Runs the GLMNET Method
# If a value for lambda is provided, the estimation is performed only for that value of lambda. It chooses the lambda value with CV otherwise.
# Alpha is the elastic net parameter. The value 1 corresponds to the LASSO penalty. The value 0 corresponds to Ridge penalty.
# It estimates the intercept only if intercept=T. It is set to zero otherwise. 
#
# Returns the corresponding beta and intercept
run.glmnet = function(X,Y,lambda = NULL, alpha = 1, intercept = T){
  if(is.null(lambda)){
    cv = cv.glmnet(x = X, y = Y,alpha = alpha, family = "binomial", intercept = intercept)
    lambda = cv$lambda.1se
  }
  coefs = coef(glmnet(x = X, y = Y, alpha = alpha, family = "binomial", intercept = intercept), s = lambda)
  # print(coefs)
  return(list(beta = coefs[-1], beta0 = coefs[1], lambda = lambda))
}

###################################
## Bianco Yohai Method
###################################

# Runs the robust Bianco and Yohai method.
#
# Returns the corresponding beta and intercept
run.by = function(X,Y){
  fit.by = BYlogreg(X,Y)
  return(list(beta = fit.by$coef[-1], beta0 = fit.by$coef[1], lambda = 0))
}

run.wby = function(X,Y){
  full.coef = glmrob(Y~X, family = "binomial", method = "WBY")$coefficients
  return(list(beta = full.coef[-1], beta0 = full.coef[1], lambda = 0))
}


###################################
## Chi Scott Method
###################################

# alpha is the elastic-net parameter (1 for LASSO, 0 for Ridge)
# 
# Returns a grid with the lambda candidates for given data and alpha
get.lambda.grid.cs = function(X,Y,alpha, len = 50){
  y.mean = mean(Y)
  factor = max(abs(t(X) %*% Y))
  lambda.max = y.mean * (1 - y.mean) * factor *2 /(nrow(X) * alpha)
  lambda.min = 0.02 * lambda.max
  log.seq = seq(from = log10(lambda.max), to = log10(lambda.min), length.out = len)
  return(10**(log.seq))
}

# Returns initial estimates for beta and the intercept
#
get.initial.values.cs = function(X,Y){
  y.mean = mean(Y)
  scores = abs(y.mean * (1- y.mean) * t(X) %*% (Y - rep(y.mean,length(Y))))
  # max.score = max(scores)
  beta.init = 1* (scores > quantile(scores,0.9))
  beta0.init = log(y.mean / (1-y.mean))
  return(list(beta.init = beta.init, beta0.init = beta0.init))
}

# Runs the Chi-Scott Method
# If a value for lambda is provided, the estimation is performed only for that value of lambda. It chooses the lambda value with CV otherwise.
# Alpha is the elastic net parameter. The value 1 corresponds to the LASSO penalty. The value 0 corresponds to Ridge penalty.
#
# Returns the corresponding beta and intercept
run.cs = function(X,Y,lambda = NULL, alpha = 1){
  initial.values = get.initial.values.cs(X = X, Y = Y)
  if(is.null(lambda)){
    lambda = robust.cv.2(X = X, Y = Y,  alpha = alpha, lambda.len = 50, tol = 1e-6)
  }
  fit = rLogistic(X = X, Y = Y, alpha = alpha, lambda = lambda, beta0 = initial.values$beta0.init, beta = initial.values$beta.init)
  return(list(beta = fit$beta, beta0 = fit$a0, lambda = lambda))    
}

###################################
## Chi Scott Modified and LASO Method
###################################


# Runs the variation of the Chi-Scott Method
# If a value for lambda is provided, the estimation is performed only for that value of lambda. It chooses the lambda value with CV otherwise.
# method equal to "CSMOD" corresponds to the penalty sum(abs(beta_j) / (1 + abs(beta_j))). method equal to "CSLASO" corresponds to the penalty norm_1(beta) / norm_2(beta).
# It estimates the intercept only if intercept=T. It is set to zero otherwise. 
#
# Returns the corresponding beta and intercept
run.cs.variations = function(method, X,Y,lambda = NULL, intercept = T){
  # if(is.null(lambda)){
  #   fit = cv.cs.variations(method = method, X = X, Y = Y, lambda.len = 15, cl.size = cl.size, intercept = intercept)
  # }
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.cs(X,Y,lambdas = lambda, intercept = intercept))
}

###################################
## BYpen Method
###################################


# Runs the Penalized Bianco-Yohai Method
# If a value for lambda is provided, the estimation is performed only for that value of lambda. It chooses the lambda value with CV otherwise.
# It estimates the intercept only if intercept=T. It is set to zero otherwise. 
#
# Returns the corresponding beta and intercept
run.bypen = function(X,Y,lambda = NULL, intercept = T, do.cv, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  
  return(cv.bypen(X = X, Y = Y, intercept = intercept, lambdas = lambda, do.cv = do.cv, nlambda = nlambda))
}

###################################
## DEVpen Method
###################################


# Runs the Penalized Deviance Method (penalizing with the quotient of the norms)
# If a value for lambda is provided, the estimation is performed only for that value of lambda. It chooses the lambda value with CV otherwise.
# It estimates the intercept only if intercept=T. It is set to zero otherwise. 
#
# Returns the corresponding beta and intercept
run.devpen = function(X,Y,lambda = NULL, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.devpen(X = X, Y = Y, intercept = intercept, lambdas = lambda, nlambda = nlambda))
}

###################################
## GLMMCP Method
###################################



# Returns the corresponding beta and intercept
run.glmmcp = function(X,Y,lambda = NULL, intercept = T){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  cv.res = cv.cvplogistic(y = Y, x = X, penalty = "mcp")
  return(list(beta = cv.res[[3]][-1], beta0 = cv.res[[3]][1], lambda = cv.res[[2]]))
}



###################################
## GLMSCAD Method
###################################



# Returns the corresponding beta and intercept
run.glmscad = function(X,Y,lambda = NULL, intercept = T){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  cv.res = cv.cvplogistic(y = Y, x = X, penalty = "scad" )
  return(list(beta = cv.res[[3]][-1], beta0 = cv.res[[3]][1], lambda = cv.res[[2]]))
}

###################################
## BYNET Method
###################################



run.bynet = function(X,Y,lambda = NULL, alpha = 1, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.bynet(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}

###################################
## BYSCAD Method
###################################



run.byscad = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.byscad(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}

###################################
## BYMCP Method
###################################



run.bymcp = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.bymcp(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}

###################################
## PREGPEN Method
###################################



run.pregpen = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.pregpen(X = X, Y = Y, intercept = intercept, lambdas = lambda, do.cv = T, nlambda = nlambda))
}

###################################
## PREGSCAD Method
###################################



run.pregscad = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.pregscad(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}

###################################
## BYMCP Method
###################################



run.pregmcp = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.pregmcp(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}

###################################
## BYPENMOD Method
###################################



run.bypenmod = function(X,Y,lambda = NULL, alpha = 1, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.bypenmod(X = X, Y = Y, intercept = intercept, lambdas = lambda, do.cv = T, nlambda = nlambda, constant = alpha))
}


###################################
## BASUSP Method
###################################


run.basusp = function(X,Y,lambda = NULL, intercept = T, do.cv, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  
  return(cv.basusp(X = X, Y = Y, intercept = intercept, lambdas = lambda, do.cv = do.cv, nlambda = nlambda))
}

###################################
## BASUSCAD Method
###################################



run.basuscad = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.basuscad(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}

###################################
## BYMCP Method
###################################



run.basumcp = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.basumcp(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}

run.hubwbypen = function(X,Y,lambda = NULL, intercept = T, do.cv, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  
  return(cv.hubwbypen(X = X, Y = Y, intercept = intercept, lambdas = lambda, do.cv = do.cv, nlambda = nlambda))
}
run.bisqwbypen = function(X,Y,lambda = NULL, intercept = T, do.cv, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  
  return(cv.bisqwbypen(X = X, Y = Y, intercept = intercept, lambdas = lambda, do.cv = do.cv, nlambda = nlambda))
}
run.hardwbypen = function(X,Y,lambda = NULL, intercept = T, do.cv, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  
  return(cv.hardwbypen(X = X, Y = Y, intercept = intercept, lambdas = lambda, do.cv = do.cv, nlambda = nlambda))
}
run.hubwbyscad = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.hubwbyscad(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}
run.hubwbymcp = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.hubwbymcp(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}
run.bisqwbyscad = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.bisqwbyscad(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}
run.bisqwbymcp = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.bisqwbymcp(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}
run.hardwbyscad = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.hardwbyscad(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}
run.hardwbymcp = function(X,Y,lambda = NULL, alpha = 2.7, intercept = T, nlambda){
  # print(c("Chosen beta : ", cv.best$beta))
  # print(c("Chosen beta0 : ", cv.best$beta0))
  # fit = cyclical.descent(method = method, X = X, Y = Y, lambda = lambda, beta.init = initial.values$beta.init)
  # print(c("Final fit parameters: ", fit$par))
  return(cv.hardwbymcp(X = X, Y = Y, intercept = intercept, alpha = alpha, lambdas = lambda, do.cv = T, nlambda = nlambda))
}

###################################
## General simulation methods
###################################


###############################################
## Run a method for given data
###############################################

#Auxiliary function to match method string with its corresponding method.
run.method = function(X,Y, method, alpha = 1, intercept = T, nlambda = 12, lambda = NULL, cl.size = NULL, do.cv = T){
  if(alpha == 1 && (identical(method, "BASUSCAD") ||identical(method, "BASUMCP") || identical(method, "BYSCAD") || identical(method, "BYMCP") || identical(method, "PREGSCAD") || identical(method, "PREGMCP"))){
    alpha = 2.7
  }
  if(identical(method, "GLM")){
    betas = run.glm(X = X, Y = Y, intercept = intercept)
  }
  if(identical(method, "GLMNET")){
    betas = run.glmnet(X = X, Y = Y, intercept = intercept,lambda = lambda, alpha = alpha)
  }
  
  ### Es el mismo que BY!!! (creo)
  if(identical(method, "WBY")){
    betas = run.wby(X = X, Y = Y)
  }
  if(identical(method, "BY")){
    betas = run.by(X = X, Y = Y)
  }
  if(identical(method, "CS")){
    betas = run.cs(X = X, Y = Y,lambda = lambda, alpha = alpha)
  }
  if(identical(method, "CSMOD")){
    betas = run.cs.variations(method = "CSMOD", X = X, Y = Y, lambda = lambda)
  }
  if(identical(method, "CSLASO")){
    betas = run.cs.variations(method = "CSLASO", X = X, Y = Y, lambda = lambda, intercept = intercept)
  }
  if(identical(method, "BYPEN")){
    # print("a3")
    betas = run.bypen(X = X, Y = Y, lambda = lambda, intercept = intercept, do.cv = do.cv, nlambda = nlambda)
  }
  if(identical(method, "DEVPEN")){
    # print("a3")
    betas = run.devpen(X = X, Y = Y, lambda = lambda, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "GLMMCP")){
    # print("a3")
    betas = run.glmmcp(X = X, Y = Y, lambda = lambda, intercept = intercept)
  }
  if(identical(method, "GLMSCAD")){
    # print("a3")
    betas = run.glmscad(X = X, Y = Y, lambda = lambda, intercept = intercept)
  }
  if(identical(method, "BYNET")){
    # print("a3")
    betas = run.bynet(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "BYSCAD")){
    # print("a3")
    betas = run.byscad(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "BYMCP")){
    # print("a3")
    betas = run.bymcp(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "PREGPEN")){
    betas = run.pregpen(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "PREGSCAD")){
    betas = run.pregscad(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "PREGMCP")){
    betas = run.pregmcp(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "BYPENMOD")){
    betas = run.bypenmod(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "BASUSP")){
    betas = run.basusp(X = X, Y = Y, lambda = lambda, intercept = intercept, do.cv = do.cv, nlambda = nlambda)
  }
  if(identical(method, "BASUSCAD")){
    betas = run.basuscad(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "BASUMCP")){
    betas = run.basumcp(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "HUBWBYPEN")){
    betas = run.hubwbypen(X = X, Y = Y, lambda = lambda, intercept = intercept, do.cv = do.cv, nlambda = nlambda)
  }
  if(identical(method, "HUBWBYSCAD")){
    if(alpha == 1){
      alpha = 3.7
    }
    betas = run.hubwbyscad(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "HUBWBYMCP")){
    if(alpha == 1){
      alpha = 2.7
    }
    betas = run.hubwbymcp(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "BISQWBYPEN")){
    betas = run.bisqwbypen(X = X, Y = Y, lambda = lambda, intercept = intercept, do.cv = do.cv, nlambda = nlambda)
  }
  if(identical(method, "BISQWBYSCAD")){
    if(alpha == 1){
      alpha = 3.7
    }
    betas = run.bisqwbyscad(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "BISQWBYMCP")){
    if(alpha == 1){
      alpha = 2.7
    }
    betas = run.bisqwbymcp(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "HARDWBYPEN")){
    betas = run.hardwbypen(X = X, Y = Y, lambda = lambda, intercept = intercept, do.cv = do.cv, nlambda = nlambda)
  }
  if(identical(method, "HARDWBYSCAD")){
    if(alpha == 1){
      alpha = 3.7
    }
    betas = run.hardwbyscad(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  if(identical(method, "HARDWBYMCP")){
    if(alpha == 1){
      alpha = 2.7
    }
    betas = run.hardwbymcp(X = X, Y = Y, lambda = lambda, alpha = alpha, intercept = intercept, nlambda = nlambda)
  }
  return(betas)
}





###############################################
## Generates data and simulates
###############################################


# Single replication method
# method corresponds to the method for estimating beta and intercept
# n is the number of individuals
# p is the number of covariates in each individual
# nonzero.beta is the true value of beta (only the nonzero part)
# model.intercept is the true value of the intercept (NULL value implies no intercept)
# mu is the mean value for the generation of each row in the X matrix.
# sigma is the covariance matrix for the generation of each row in the X matrix.
# eps is the proportion of outliers that are added (and not replaced) in the data
#
# Two types of contamination are available:
# First type: fixed contamination
# out.X is a given X matrix of outliers (will be attached to the X matrix generated without outliers)
# out.Y is a given Y vector of outliers (will be attached to the Y vector generated without outliers)
#
# Second type: random contamination
# out.mu is the mean value for the generation of each row in the outlier matrix
# out.sigma is the covariance matrix for the generation of each row in the outlier matrix
# The corresponding Y value is the opposite of the most likely value that each X row would generate
# n.test is the number of OOS individuals to estimate errors
# alpha is the elastic net parameter
# method.intercept indicates if the intercept should be estimated. If FALSE, the intercept is set to zero.
# if a value for lambda is provided, it estimates beta and the intercept with that penalization parameter. It chooses lambda with CV otherwise.

simul.method = function(method, n,p,nonzero.beta = rep(1,5) ,model.intercept = NULL,
                        mu = NULL, sigma = NULL, eps = 0, out.X = NULL, out.Y = NULL,
                        out.mu = NULL, out.sigma = NULL, n.test = NULL, alpha = 1, seed = 1234,
                        method.intercept = T, lambda = NULL, nlambda = 12, overlapping = NA, cl.size = NULL, do.cv = T){
  set.seed(seed)
  data.train = get.data.norm(n = n, p = p,mu = mu, sigma = sigma, nonzero.beta = nonzero.beta, eps = eps, intercept = model.intercept, out.X = out.X, out.Y = out.Y, out.mu = out.mu, out.sigma = out.sigma)
  # print("a1")
  overlapping = T
  
  # if(is.na(overlapping)){
  #   overlapping = noverlap.for(cbind(data.train$X,data.train$Y))$NOVERLAP != 0
  # }
  # 
  # if(!overlapping){
  #   return(rep("NO",5))
  # }
  
  if(is.null(n.test)){
    n.test = n
  }
  
  data.test = get.data.norm(n = n.test, p = p,mu = mu, sigma = sigma, nonzero.beta = nonzero.beta, eps = 0, intercept = model.intercept)
  # print("a2")
  
  fit = run.method(X = data.train$X, Y = data.train$Y, method = method,alpha = alpha, intercept = method.intercept, lambda = lambda, cl.size = cl.size, do.cv = do.cv, nlambda = nlambda)
  errors = get.errors(X = data.train$X ,Y = data.train$Y ,probs = data.train$probs, new.X = data.test$X, new.Y = data.test$Y, new.probs = data.test$probs, beta.est = fit$beta, beta0.est = fit$beta0, nonzero.beta.true = nonzero.beta)
  return(list(probs.mse.oos = errors$probs.mse.oos, correct.prop.oos = errors$correct.prop.oos, probs.mse.is = errors$probs.mse.is, correct.prop.is = errors$correct.prop.is, beta.se = errors$beta.se, tn = errors$tn, tp = errors$tp,
              lambda = fit$lambda, beta0 = fit$beta0, beta = fit$beta, 
              all.lambdas = fit$all.lambdas, all.betas = fit$all.betas, all.betas0 = fit$all.betas0, all.df = fit$all.df, all.loss.eval = fit$all.loss.eval ))
} 


###################################
## Run
###################################


###############################################
## Run several iterations for simulation
###############################################
# Runs several iterations for a given method
# niter is the number of iterations to perform
# other parameters are described in simul.method
# mu.code, sigma.code, oux.X.code, out.Y.code, out.mu.code, out.sigma.code are codes to generate the corresponding matrices or vectors.

simul.many = function(method, n, p, nonzero.beta = rep(1,5) ,niter = 1000,mu.code = NULL, sigma.code = NULL, eps = 0, model.intercept = NULL, out.X.code = NULL, out.Y.code = NULL, out.mu.code = NULL,
                      out.sigma.code = NULL, n.test = NULL, alpha = 1, lambda = NULL, method.intercept = T, overlapping = NULL, cl.size = NULL, start.index = 1){
  dir = "results"
  if(!dir.exists(dir)){
    dir.create(dir)
  }
  file = get.file.name(method = method, n = n, p =p, nonzero.beta = nonzero.beta, mu.code = mu.code, sigma.code = sigma.code, eps = eps, intercept = model.intercept, out.X.code = out.X.code, out.Y.code = out.Y.code, out.mu.code = out.mu.code, out.sigma.code = out.sigma.code,n.test = n.test, alpha = alpha, lambda = lambda)
  print(file)
  
  ##TODO : Create codes
  
  out.X = NULL
  mu = rep(mu.code,p)
  if(is.null(sigma.code)){
    sigma = NULL
  } else{
    sigma = sigma.code * diag(p)
  }
  if(!is.null(out.X.code)){
    out.X = c(rep(out.X.code * sqrt(p) ,length(nonzero.beta)), rep(0, p - length(nonzero.beta)))
  }
  out.Y = out.Y.code
  out.mu = rep(out.mu.code, p)
  if(is.null(out.sigma.code)){
    out.sigma = NULL
  } else{
    out.sigma = out.sigma.code * diag(p)
  }
  
  if(!is.null(sigma)){
    beta = c(nonzero.beta, rep(0,p-length(nonzero.beta)))
    var.cl = t(beta) %*% sigma %*% beta
  } else{
    var.cl = sum(nonzero.beta**2)
  }
  correct.prop.optim = get.best.prob(var.cl)
  
  if(is.null(overlapping)){
    overlapping = rep(NA, niter)
  }
  
  for(i in start.index:niter){
    print(c("File", file, "Iter n. ", i, sep = " "))
    res = try(simul.method(method, n,p,nonzero.beta,model.intercept =  model.intercept, mu = mu, sigma = sigma, eps = eps, out.X = out.X, out.Y = out.Y, out.mu = out.mu, out.sigma = out.sigma, n.test = n.test, alpha = alpha, seed = 1234 + i, method.intercept = method.intercept, lambda = lambda, overlapping = overlapping[i], cl.size = cl.size))
    if(class(res) != "try-error"){
      # File printing for error measures
      write(x = c(i,correct.prop.optim$value,unlist(res[1:7])), file = paste(dir,file,sep="/"), ncolumns = 9, append = T, sep = " " )
      # File printing for fit (lambda + beta0 + beta)
      write(x = c(i,unlist(res[8:10])), file = paste(dir, paste(file, "fit", sep = "_"), sep = "/"), ncolumns = p + 3, append = T, sep = " ")
      # File printing fot fit for all lambdas: (lambda + beta0 + beta)
      lambdas = res$all.lambdas
      if(!is.null(lambdas)){
        for(j in 1:length(lambdas)){
          lam = lambdas[j]
          beta0 = res$all.betas0[j]
          beta = res$all.betas[j,]
          df = res$all.df[j]
          loss.eval = res$all.loss.eval[j]
          write(x = c(i,lam,beta0,beta,df,loss.eval), file = paste(dir, paste(file, "fit", sep = "_"), sep = "/"), ncolumns = p + 5, append = T, sep = " ")
        }
      }
    }
  }
}


###############################################
## Return the file names given the parameters
###############################################
get.file.name = function(method, n, p, nonzero.beta,mu.code, sigma.code, eps, intercept, out.X.code, out.Y.code, out.mu.code, out.sigma.code, n.test, alpha, lambda){
  res = paste(method,"n",n,"p",p,"k",length(nonzero.beta), sep = "_")
  if(!is.null(mu.code)){
    res = paste(res,"mu", mu.code, sep = "_")
  }
  if(!is.null(sigma.code)){
    res = paste(res,"sigma", sigma.code, sep = "_")
  }
  if(!is.null(intercept)){
    res = paste(res,"int", intercept, sep = "_")
  }
  if(!is.null(eps)){
    res = paste(res,"eps", eps, sep = "_")
  }
  if(!is.null(out.X.code)){
    res = paste(res,"outx", out.X.code, sep = "_")
  }
  if(!is.null(out.Y.code)){
    res = paste(res,"outy", out.Y.code, sep = "_")
  }
  if(!is.null(out.mu.code)){
    res = paste(res,"outmu", out.mu.code, sep = "_")
  }
  if(!is.null(out.sigma.code)){
    res = paste(res,"outsigma", out.sigma.code, sep = "_")
  }
  if(!is.null(alpha)){
    res = paste(res, "alpha",alpha, sep = "_")
  }
  if(!is.null(lambda)){
    res = paste(res, "lambda",lambda, sep = "_")
  }
  if(!is.null(n.test)){
    res = paste(res, "ntest",n.test, sep = "_")
  }
  return(res)
}

################################
## Auxiliary functions
################################

# Computes error measures
# new.X, new.Y and new.probs corresponds to the new data, for testing out of sample errors
# beta.est and beta0.est are the estimates for beta and the intercept
# nonzero.beta.true is the true value of beta (only the non-zero part)
#
# probs.mse corresponds to the in-sample mean squared error for the probabilities
# correct.prop corresponds to the out for sample estimate of correct classification
# beta.se corresponds the 2-norm of the difference of the true beta and its estimate
# tp corresponds to the proportion of active variable that are estimated as active
# tn corresponds to the proportion of non-active variables that are estimated as non-active
#
get.errors = function(X,Y,probs,new.X,new.Y,new.probs,beta.est,beta0.est, nonzero.beta.true){
  scores.oos = cbind(1,new.X) %*% c(beta0.est,beta.est)
  probs.est.oos = inv.logit(scores.oos)
  class.est.oos = 1*(probs.est.oos > (1/2))
  probs.mse.oos = mean((new.probs - probs.est.oos)**2)
  correct.prop.oos = sum(class.est.oos == new.Y) / length(new.Y)
  
  scores.is = cbind(1,X) %*% c(beta0.est,beta.est)
  probs.est.is = inv.logit(scores.is)
  class.est.is = 1*(probs.est.is > (1/2))
  probs.mse.is = mean((probs - probs.est.is)**2)
  correct.prop.is = sum(class.est.is == Y) / length(Y)
  
  beta.error = get.beta.error(beta.est = beta.est, nonzero.beta = nonzero.beta.true)
  return(list(probs.mse.oos = probs.mse.oos, correct.prop.oos = correct.prop.oos, probs.mse.is = probs.mse.is, correct.prop.is = correct.prop.is, beta.se = beta.error$sq.error, tn = beta.error$tn, tp = beta.error$tp))
}

# Returns errors from comparing true beta with estimated beta
get.beta.error = function(beta.est, nonzero.beta){
  p = length(beta.est)
  k = length(nonzero.beta)
  beta = c(nonzero.beta, rep(0,p-k))
  sq.error = sum((beta.est - beta)^2)
  active.set = 1:k
  active.set.est = which(beta.est != 0)
  true.pos = length(which(is.element(active.set.est,active.set)))
  return(list(sq.error = sq.error, tp = true.pos / k, tn = (p-k - (length(active.set.est) - true.pos))/(p-k)))
}

# Returns the probability of correct classification in the case of knowing the true beta.
# TO DO: check. Maybe Monte-Carlo is better.

get.best.prob = function(sigma.sq){
  fun = function(x){
    return(2*(1/sqrt(2*pi*sigma.sq))*exp(-(x^2/(2*sigma.sq)) + x) / (1 + exp(x)))
  }
  return(integrate(fun, lower = 0, upper = Inf))
}

# Implementation of the function f(s) = (1+exp(-s))^(-1)
inv.logit = function(u){
  return(ifelse(test = (u > 16),yes = 1,no = ifelse(test = (u < -16), yes = 0, no = (1+exp(-u))**(-1))))
}

# Predicts probabilities for given X, beta and intercept
predict.probs = function(beta,X, beta0 = 0){
  scores = X %*% beta + beta0
  return(inv.logit(scores))
}

separate.folds = function(total, nfolds = 10){
  shuffled = sample(1:total,size = total, replace = F)
  return(split(shuffled,cut(1:total,breaks = nfolds)))
}


#############
# Penalties
#############

# Returns the quotient norm1(vec)/norm2(vec)
norm.quotient = function(beta, lambda){
  return(norm_quotient_cpp(beta = beta, lambda = lambda))
}

norm.quotient.mod = function(beta, lambda, constant = 1){
  return(norm_quotient_mod_cpp(beta = beta, lambda = lambda, constant = constant))
}

elastic.net = function(beta, alpha, j = NULL, lambda){
  if(is.null(j)){
    return(lambda * alpha * sum(abs(beta)) + lambda * ((1-alpha)/2) * sum(beta^2))
  } else{
    return(lambda * alpha * abs(beta[j]) + lambda * ((1-alpha)/2) * beta[j]^2)
  }
}

scad = function(beta, alpha, j = NULL, lambda){
  if(is.null(j)){
    return(sum(mapply(FUN = scad_cpp, t = beta, lambda = lambda , gamma = alpha)))
  } else{
    return(scad_cpp(t = beta[j], lambda = lambda, gamma = alpha))
  }
}

mcp = function(beta, alpha, j = NULL, lambda){
  if(is.null(j)){
    return(sum(mapply(FUN = mcp_cpp, t = beta, lambda = lambda , gamma = alpha)))
  } else{
    return(mcp_cpp(t = beta[j], lambda = lambda, gamma = alpha))
  }
}

# Returns the penalty for the CSMOD method
norm.mod = function(beta, lambda){
  return(lambda * sum(abs(beta) / (1 + abs(beta))))
}





# #### TRUE means overlapping
# #### A bit hardcoded
# 
# analyze.overlapping = function(data.list){
#   
#   ret = rep(NA,length(data.list))
#   for(i in 1:length(data.list)){
#     
#     if(nrow(data.list[[i]]$X) > 100){
#       ret[i] = TRUE
#     } else if(nrow(data.list[[i]]$X) > ncol(data.list[[i]]$X)){
#       ret[i] = noverlap.for(cbind(data.list[[i]]$X,data.list[[i]]$Y))$NOVERLAP != 0
#     } else{
#       ret[i] = FALSE
#     }
#   }
#   return(ret)
# }

# proportion: the proportion of covariates that will be set as active
# full.fit: if TRUE, the BY method is applied to the set of selected covariates. If FALSE, the initial estimates are the univariate fits for the active variables.
get.robust.initial = function(X,Y, proportion = 0.1, full.fit = T){
  n = nrow(X)
  p = ncol(X)
  signif.values = rep(NA,p)
  univ.fits = rep(NA,p)
  for(i in 1:p){
    univ.fit = glmrob(Y ~ X[,i] , family = "binomial", method = "BY")
    print(i)
    # print(summary(univ.fit)$coefficients)
    signif.values[i] = summary(univ.fit)$coefficients[2,4]
    if(!full.fit){
      univ.fits[i] = summary(univ.fit)$coefficients[2,1]
    }
  }
  
  # print(signif.values)
  n.active = ceiling(p * proportion)
  active.set = order(signif.values)[1:n.active]
  active.set = sort(active.set)
  beta = rep(0,p)
  if(!full.fit){
    for(i in 1:p){
      if(i %in% active.set){
        beta[i] = univ.fits[i]
      }
    }
    beta0 = glmrob(Y ~ 1, family = "binomial", method = "WBY")$coefficients[1]
    return(list(beta.init = beta, beta0.init = beta0))
  }
  else{
    full.fit.coef = glmrob(Y~X[,active.set], family = "binomial", method = "WBY")$coefficients
    beta0 = full.fit.coef[1]
    beta[active.set] = full.fit.coef[-1]
    return(list(beta.init = beta, beta0.init = beta0))
  }
}

get.robust.initial.2 = function(X,Y, prop = 0.1, trim = 0.15){
  n = nrow(X)
  p = ncol(X)
  scores = rep(NA,p)
  p.bar = mean(Y)
  for(j in 1:p){
    scores[j] = p.bar * (1-p.bar) * abs(mean(X[,j] * (Y - p.bar), trim = trim))
  }
  # print(scores)
  n.active = ceiling(p * prop)
  active = order(-scores)[1:n.active]
  ret = rep(0,p)
  ### Poner un Try-catch aca!!
  full.fit.coef = try(glmrob(Y~X[,active], family = "binomial", method = "WBY")$coefficients)
  if(class(full.fit.coef) == "try-error" || is.na(full.fit.coef[1])){
    full.fit.coef = c(0,rep(1,n.active))
  }
  ret[active] = full.fit.coef[-1]
  return(list(beta.init = ret, beta0.init = full.fit.coef[1]))
}


### Weights function

get.rob.hd.cov.estimate = function(X){
  p = ncol(X)
  kendall.cor = cor(X, method = "kendall")
  mads = apply(X = X, MARGIN = 2, FUN = mad)
  kendall.cov = matrix(NA, p,p)
  for(j in 1:p){
    for(k in 1:j){
      kendall.cov[j,k] = kendall.cov[k,j] = mads[k] * mads[j] * kendall.cor[j,k]
    }
  }
  return(kendall.cov)
}

get.weights.box = function(X){
  cov = get.rob.hd.cov.estimate(X)
  med = Gmedian(X)
  prec = glasso(cov, rho = 0.1)$w
  dist = rep(NA, nrow(X))
  for(i in 1: nrow(X)){
    dist[i] = sqrt((X[i,] - med) %*% prec %*% t((X[i,] - med)))
  }
  q3 = quantile(dist,0.75)
  q1 = quantile(dist,0.25)
  lbox = q3 - q1
  weights = rep(NA,nrow(X))
  for(i in 1: nrow(X)){
    weights[i] = 1*(dist[i] < q3 + 1.5*lbox)
  }
  return(weights)
}

get.weights = function(X, method, constant){
  # if(ncol(X) >= nrow(X)){
  #   return(rep(1,nrow(X)))
  # }
  
  ## Esto no funciona xq invierte la matriz de covarianza
  #return(getWeight(OutlierMahdist(X)))
  cov = get.rob.hd.cov.estimate(X)
  med = Gmedian(X)
  prec = glasso(cov, rho = 0.1)$w
  dist = rep(NA, nrow(X))
  # print(prec)
  # print(X[1,] - med)
  for(i in 1: nrow(X)){
    dist[i] = sqrt((X[i,] - med) %*% prec %*% t((X[i,] - med)))
  }
  print(dist)
  boxplot(dist)$out
  if(identical(method, "HARD")){
    return(1*(dist < constant))
  } else if(identical(method, "BISQUARE")){
    return(sapply(X = dist, FUN = bisquare.weight, constant = constant))
  } else if(identical(method, "HUBER")){
    return(sapply(X = dist, FUN = huber.weight, constant = constant))
  } else{
    stop("Enter a valid method")
  }
}


## Huber weights function

huber.weight = function(x, constant){
  
  if(abs(x) <= constant){
    return(1)
  } else{
    return(constant / abs(x))
  }
  
}



bisquare.weight = function(x, constant){
  if(abs(x) <= constant){
    return((1-(x/constant)**2)**2)
  } else{
    return(0)
  }
}


####################
## Functions for real data
###################


separate.folds = function(total, nfolds){
  shuffled = sample(1:total,size = total, replace = F)
  return(split(shuffled,cut(1:total,breaks = nfolds)))
}


split.data = function(X,Y, folds = 5){
  ret = list()
  X.ones = X[which(Y == 1),]
  X.zeros = X[which(Y == 0),]
  ones.folds.ind = separate.folds(nrow(X.ones), folds)
  zeros.folds.ind = separate.folds(nrow(X.zeros), folds)
  # print(ones.folds.ind)
  # print(zeros.folds.ind)
  for(i in 1: folds){
    ret[[i]] = list()
    ret[[i]]$X = rbind(X.ones[ones.folds.ind[[i]],], X.zeros[zeros.folds.ind[[i]],])
    ret[[i]]$Y = c(rep(1,length(ones.folds.ind[[i]])), rep(0,length(zeros.folds.ind[[i]])))
  }
  return(ret)
}

get.ccp = function(res, X.test, Y.test){
  return(sum((X.test %*% res$beta + res$beta0 > 0) == Y.test) / nrow(X.test))
}

standarize = function(vec){
  return((vec - mean(vec)) / sd(vec))
}

standarize.X = function(X){
  
  for(i in 1:ncol(X)){
    if(all(X[,i] == X[1,i])){
      X[,i] = rep(0,nrow(X))
    } else{
      X[,i] = standarize(X[,i])
    }
  }
  return(X)
}


cv.data = function(method,X,Y, alpha = 1){
  set.seed(1234)
  print(c("Method: ", method))
  betas = list()
  ccp = betas0 = c()
  split = split.data(X,Y)
  nfolds = length(split)
  for(i in 1:nfolds){
    print(c("Fold: ", i))
    X.train = data.frame()
    Y.train = c()
    for(j in 1:nfolds){
      if(j != i){
        X.train = rbind(X.train,split[[j]]$X)
        Y.train = c(Y.train, split[[j]]$Y)
      }
    }
    res.fold = run.method(X = X, Y = Y, method = method, alpha = alpha)
    betas[[i]] = res.fold$beta
    betas0[i] = res.fold$beta0
    ccp[i] = get.ccp(res.fold, split[[i]]$X, split[[i]]$Y)
    print(ccp[i])
  }
  return(list(betas = betas, ccp = ccp, betas0 = betas0))
}


get.table.name.from.method.code = function(method, lambda, alpha){
  print(method)
  print(lambda)
  print(alpha)
  if(identical(method, "GLM")){
    return("$\\wbbe_{GLM}$")
  }
  if(identical(method, "BY")){
    return("$\\wbbe_{BY}$")
  }
  if(identical(method, "CS") && !is.na(lambda) && lambda == 0){
    return("$\\wbbe_{CS,0}$")
  }
  if(identical(method, "GLMNET") && !is.na(alpha) && alpha == 0){
    return("$\\wbbe_{NET,0}$")
  }
  if(identical(method, "GLMNET") && !is.na(alpha)&& alpha == 0.5){
    return("$\\wbbe_{NET,0.5}$")
  }
  if(identical(method, "GLMNET") && !is.na(alpha)&& alpha == 0.75){
    return("$\\wbbe_{NET,0.75}$")
  }
  if(identical(method, "GLMNET") && !is.na(alpha)&& alpha == 1){
    return("$\\wbbe_{NET,1}$")
  }
  if(identical(method, "CS") && !is.na(alpha)&& alpha == 0.6){
    return("$\\wbbe_{CS,0.6}$")
  }
  if(identical(method, "CS") && !is.na(alpha) && alpha == 1){
    return("$\\wbbe_{CS,1}$")
  }
  if(identical(method, "BYPEN")){
    return("$\\wbbe_{BY,SP}$")
  }
  if(identical(method, "DEVPEN")){
    return("$\\wbbe_{GLM,SP}$")
  }
  if(identical(method, "CSLASO")){
    return("$\\wbbe_{CS,SP}$")
  }
  if(identical(method, "GLMSCAD")){
    return("$\\wbbe_{GLM,SCAD}$")
  }
  if(identical(method, "GLMMCP")){
    return("$\\wbbe_{GLM,MCP}$")
  }
  if(identical(method, "BYSCAD")){
    return("$\\wbbe_{BY,SCAD}$")
  }
  if(identical(method, "BYMCP")){
    return("$\\wbbe_{BY,MCP}$")
  }
  if(identical(method, "PREGPEN")){
    return("$\\wbbe_{PREG,SP}$")
  }
  if(identical(method, "PREGSCAD")){
    return("$\\wbbe_{PREG,SCAD}$")
  }
  if(identical(method, "PREGMCP")){
    return("$\\wbbe_{PREG,MCP}$")
  }
  if(identical(method, "BASUSP")){
    return("$\\wbbe_{BASU,SP}$")
  }
  if(identical(method, "BASUSCAD")){
    return("$\\wbbe_{BASU,SCAD}$")
  }
  if(identical(method, "BASUMCP")){
    return("$\\wbbe_{BASU,MCP}$")
  }
  if(identical(method, "HUBWBYPEN")){
    return("$\\wbbe_{BY,WHUB,SP}$")
  }
  if(identical(method, "HUBWBYMCP")){
    return("$\\wbbe_{BY,WHUB,MCP}$")
  }
  if(identical(method, "HARDWBYPEN")){
    return("$\\wbbe_{BY,WHARD,SP}$")
  }
  if(identical(method, "HARDWBYMCP")){
    return("$\\wbbe_{BY,WHARD,MCP}$")
  }
  if(identical(method, "BYPENMOD") && !is.na(alpha) && alpha == 0.1){
    return("$\\wbbe_{BY,SP,0.1}$")
  }
  if(identical(method, "BYPENMOD") && !is.na(alpha) && alpha == 0.5){
    return("$\\wbbe_{BY,SP,0.5}$")
  }
  if(identical(method, "BYPENMOD") && !is.na(alpha) && alpha == 1){
    return("$\\wbbe_{BY,SP,1}$")
  }
  if(identical(method, "BYPENMOD") && !is.na(alpha) && alpha == 3){
    return("$\\wbbe_{BY,SP,3}$")
  }
}
