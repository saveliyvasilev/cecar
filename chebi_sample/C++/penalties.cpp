#include <RcppArmadillo.h>
#include <stdio.h>


using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
float mcp_cpp(float t, float lambda, float gamma){
	if(fabs(t) < (lambda * gamma)){
    return (2*lambda * gamma * fabs(t) - pow(t,2)) / (2 * gamma);
	} else {
		return pow(lambda,2) * gamma / 2;
	}
}

// [[Rcpp::export]]
float scad_cpp(float t, float lambda, float gamma){
  if(fabs(t) < (lambda)){
    return lambda * fabs(t);
  } else{
    if(fabs(t) < gamma * lambda){
      return (gamma * lambda * fabs(t) - 0.5 * (pow(t,2) + pow(lambda,2))) / (gamma - 1);
    } else{
      return pow(lambda,2) * (pow(gamma,2) -1) / (2 * (gamma - 1));
    }
  }
}

// [[Rcpp::export]]
float norm_quotient_cpp(arma::vec beta, float lambda){
  float n1 = norm(beta,1);
  if(n1 == 0){
    return 1;
  } else{
    return lambda * n1 / norm(beta,2);
  }
}

// [[Rcpp::export]]
float norm_quotient_mod_cpp(arma::vec beta, float lambda, float constant){
  return lambda * norm(beta,1) / (norm(beta,2) + constant);
}
