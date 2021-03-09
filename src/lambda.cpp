#include "lambda.h"

RcppExport SEXP cdraw_lambda_i(SEXP lambda, SEXP mean, SEXP kmax, SEXP thin) {
  arn gen;
  return Rcpp::wrap(draw_lambda_i(Rcpp::as<double>(lambda), 
                                  Rcpp::as<double>(mean),
                                  Rcpp::as<int>(kmax), 
                                  Rcpp::as<int>(thin), gen));
}

//draw lambda from its (infinite mixture) prior
double draw_lambda_prior(double *psii, int kmax, rn& gen)
{
  double lambda;
  int k;
  
  lambda = 0.0;
  for(k=0; k<=kmax; k++) {
    lambda += psii[k] * gen.exp();
  }
  
  return lambda;
}


//Metropolis-Hastings algorithm for drawing lambda from its full conditional -- uses proposals from the prior
double draw_lambda_i(double lambda_old, double xbeta,
                     int kmax, int thin, rn& gen)
{
  int t, k;
  double lambda, lp, lpold, m, s;
  double *psii;
  
  //calculate the probability og the previous lambda
  s = sqrt(lambda_old);
  m = xbeta;
  lpold = pnorm(0.0, m, s, 0, 1);
  
  //allocate psii
  psii = (double*) malloc(sizeof(double) * (kmax+1));
  for(k=0; k<=kmax; k++) psii[k] =  2.0/((1.0+k)*(1.0+k));
  
  //thinning is essential when kappa is large
  for(t=0; t<thin; t++) {
    
    //propose a new lambda from the prior
    lambda = draw_lambda_prior(psii, kmax, gen);
    
    //calculate the probability of the propsed lambda
    s = sqrt(lambda);
    m = xbeta;
    lp = pnorm(0.0, m, s, 0, 1);
    
    //MH accept or reject
    if(gen.uniform() < exp(lp - lpold)) {
      lambda_old = lambda;
      lpold = lp;
    }
  }
  
  //possibly clean up psii
  free(psii);
  
  return lambda_old;
}