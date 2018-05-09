#include <RcppArmadillo.h>
#include <stdio.h>


using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]




// [[Rcpp::export]]
float aux_fun_1(float s){
	if(s > 100){
		return s;
	} else{
		return log(1 + exp(s));
	}
}


// [[Rcpp::export]]
float inv_logit_cpp(float u){
	if(u > 16){
		return 1;
	}
	if(u < -16){
		return 0;
	}
	return 1/(1+exp(-u));
}

// [[Rcpp::export]]
float phi_cs_cpp(float score, int y){
  return pow(y - inv_logit_cpp(score),2);
}

// [[Rcpp::export]]
float der_phi_cs_cpp(float score, int y){
	float prob = inv_logit_cpp(score);
	return -2 * (y-prob) * prob * (1-prob);
}

// [[Rcpp::export]]
arma::vec get_scores(arma::mat X, arma::vec beta, float beta0){
	return(X*beta + beta0);
}

// [[Rcpp::export]]
float eval_loss_function_cs_cpp(arma::mat X, arma::vec Y, arma::vec beta, float beta0){
	int n = X.n_rows;
	arma::vec scores = get_scores(X,beta,beta0);
	float acum = 0;
	int i;
	for(i=0;i < n; i++){
		acum += phi_cs_cpp(scores[i], Y[i]);
	}
	return acum/n;
}

// [[Rcpp::export]]
arma::vec eval_der_loss_function_cs_cpp(arma::mat X, arma::vec Y, arma::vec beta, float beta0){
	int n = X.n_rows;
	int p = X.n_cols;
	arma::vec scores = get_scores(X,beta,beta0);
	arma::vec acum = zeros<vec>(p) ;
	int i;
	for(i=0;i < n; i++){
		acum += der_phi_cs_cpp(scores[i], Y[i]) * X.row(i).t();
	}
	return acum/n;
}



