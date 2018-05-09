#include <RcppArmadillo.h>
#include <stdio.h>


using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

const float PII = 3.14159265358979323846;


// [[Rcpp::export]]
float pnorm_cpp(float x)
{
    // constants
    float a1 =  0.254829592;
    float a2 = -0.284496736;
    float a3 =  1.421413741;
    float a4 = -1.453152027;
    float a5 =  1.061405429;
    float p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    float t = 1.0/(1.0 + p*x);
    float y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}


// [[Rcpp::export]]
float aux_fun_1(float s){
	if(s > 100){
		return s;
	} else{
		return log(1 + exp(s));
	}
}

// [[Rcpp::export]]
float rho_bypen_cpp(float t, float k = 0.5){
  if(t <= k){
    return t*exp(-sqrt(k));
  } else{
    return -2*exp(-sqrt(t))*(1+sqrt(t)) + exp(-sqrt(k)) * (2 * (1+sqrt(k)) + k);
  }
}

// [[Rcpp::export]]
float dev_cpp(float score, int y){
  float aux = aux_fun_1(score);
	return -y * (score - aux) + (1-y) * aux;
}

// [[Rcpp::export]]
float der_rho_bypen_cpp(float t, float k = 0.5){
  if(t <= k){
    return exp(-sqrt(k));
  } else{
    return exp(-sqrt(t));
  }
}

// [[Rcpp::export]]
float G_bypen_cpp(float t, float k = 0.5){
  if(t <= exp(-k)){
    return t*exp(-sqrt(-log(t))) + exp(0.25)*sqrt(PII) * (pnorm_cpp(sqrt(2)*(0.5 + sqrt(-log(t)))) - 1);
  } else{
    return exp(-sqrt(k))*t + exp(0.25)*sqrt(PII)*(pnorm_cpp(sqrt(2)*(0.5 + sqrt(k)))-1);
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
float phi_bypen_cpp(float score, int y){
  arma_rng::set_seed(1);
  float prob = inv_logit_cpp(score);
  return rho_bypen_cpp(dev_cpp(score,y)) + G_bypen_cpp(prob) + G_bypen_cpp(1-prob) - G_bypen_cpp(1);
}

// [[Rcpp::export]]
float phi_pregpen_cpp(float score, int y){
  arma_rng::set_seed(1);
  return rho_bypen_cpp(dev_cpp(score,y));
}


// [[Rcpp::export]]
float der_phi_0_bypen_cpp(float score){
  float prob = inv_logit_cpp(score);
  float aux = aux_fun_1(score);
  return der_rho_bypen_cpp(aux) * pow(prob,2) + der_rho_bypen_cpp(-score + aux) * prob * (1 - prob);
}

// [[Rcpp::export]]
float der_phi_0_pregpen_cpp(float score){
  float prob = inv_logit_cpp(score);
  float aux = aux_fun_1(score);
  return der_rho_bypen_cpp(aux) * prob;
}

// [[Rcpp::export]]
float der_phi_bypen_cpp(float score, int y){
  if(y == 0){
    return der_phi_0_bypen_cpp(score);
  } else{
    return -der_phi_0_bypen_cpp(-score);
  }
}

// [[Rcpp::export]]
float der_phi_pregpen_cpp(float score, int y){
  if(y == 0){
    return der_phi_0_pregpen_cpp(score);
  } else{
    return -der_phi_0_pregpen_cpp(-score);
  }
}

// [[Rcpp::export]]
arma::vec get_scores(arma::mat X, arma::vec beta, float beta0){
	return(X*beta + beta0);
}

// [[Rcpp::export]]
float eval_loss_function_bypen_cpp(arma::mat X, arma::vec Y, arma::vec beta, float beta0){
	int n = X.n_rows;
	arma::vec scores = get_scores(X,beta,beta0);
	float acum = 0;
	int i;
	for(i=0;i < n; i++){
		acum += phi_bypen_cpp(scores[i], Y[i]);
	}
	return acum/n;
}

// [[Rcpp::export]]
float eval_loss_function_wbypen_cpp(arma::mat X, arma::vec Y, arma::vec beta, float beta0, arma::vec weights){
  int n = X.n_rows;
  arma::vec scores = get_scores(X,beta,beta0);
  float acum = 0;
  int i;
  for(i=0;i < n; i++){
    acum += phi_bypen_cpp(scores[i], Y[i]) * weights[i];
  }
  return acum/n;
}

// [[Rcpp::export]]
float eval_loss_function_pregpen_cpp(arma::mat X, arma::vec Y, arma::vec beta, float beta0){
  int n = X.n_rows;
  arma::vec scores = get_scores(X,beta,beta0);
  float acum = 0;
  int i;
  for(i=0;i < n; i++){
    acum += phi_pregpen_cpp(scores[i], Y[i]);
  }
  return acum/n;
}

// [[Rcpp::export]]
arma::vec eval_der_loss_function_bypen_cpp(arma::mat X, arma::vec Y, arma::vec beta, float beta0){
	int n = X.n_rows;
	int p = X.n_cols;
	arma::vec scores = get_scores(X,beta,beta0);
	arma::vec acum = zeros<vec>(p) ;
	int i;
	for(i=0;i < n; i++){
		acum += der_phi_bypen_cpp(scores[i], Y[i]) * X.row(i).t();
	}
	return acum/n;
}

// [[Rcpp::export]]
arma::vec eval_der_loss_function_wbypen_cpp(arma::mat X, arma::vec Y, arma::vec beta, float beta0, arma::vec weights){
  int n = X.n_rows;
  int p = X.n_cols;
  arma::vec scores = get_scores(X,beta,beta0);
  arma::vec acum = zeros<vec>(p) ;
  int i;
  for(i=0;i < n; i++){
    acum += der_phi_bypen_cpp(scores[i], Y[i]) * weights[i] * X.row(i).t();
  }
  return acum/n;
}

// [[Rcpp::export]]
arma::vec eval_der_loss_function_pregpen_cpp(arma::mat X, arma::vec Y, arma::vec beta, float beta0){
  int n = X.n_rows;
  int p = X.n_cols;
  arma::vec scores = get_scores(X,beta,beta0);
  arma::vec acum = zeros<vec>(p) ;
  int i;
  for(i=0;i < n; i++){
    acum += der_phi_pregpen_cpp(scores[i], Y[i]) * X.row(i).t();
  }
  return acum/n;
}




