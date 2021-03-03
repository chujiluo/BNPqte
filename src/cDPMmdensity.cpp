#include "dpm.h"

// [[Rcpp::export]]
Rcpp::List cDPMmdensity (
    const arma::uword n,
    const arma::uword d,
    const arma::mat & data,   // nxd matrix
    const arma::colvec & y,   // vector of length n
    const arma::mat & x,  // nx(d-1) matrix
    const bool diag,
    const bool pdf,
    const bool cdf,
    const arma::uword ngrid,
    const arma::colvec & grid,
    const arma::uword npred,
    const arma::mat & xpred,
    const bool updateAlpha,
    const bool useHyperpriors,
    double alpha,
    const double a0,
    const double b0,
    double lambda,
    const double gamma1,
    const double gamma2,
    const int nu0,
    const int nu,
    const arma::uword nclusters,
    const arma::uword nskip,
    const arma::uword ndpost,
    const arma::uword keepevery,
    const arma::rowvec & diri,
    const arma::colvec & probs,
    const arma::uword nprobs
) {
  //------------------------------------------------------------------
  // process args
  arma::colvec m0(d);
  arma::mat S0(d, d, arma::fill::zeros);
  arma::mat Psi0(d, d);
  arma::colvec m(d);
  arma::mat Psi(d, d);
  
  arma::mat Zeta(d, nclusters); // each column is of length d
  arma::cube Omega(d, d, nclusters); // each slice is dxd
  arma::cube cholOmega(d, d, nclusters); // cholOmega.slice(i) = arma::chol(Omega.slice(i))
  arma::cube icholOmega(d, d, nclusters); // icholOmega.slice(i) = arma::inv(cholOmega.slice(i))
  arma::colvec othersOmega(nclusters); // terms excluding data in the log pdf of Normal(Zeta.col(i), Omega.slice(i))
  arma::rowvec a_gd(nclusters-1, arma::fill::ones);
  arma::rowvec b_gd(nclusters-1, arma::fill::ones);
  arma::rowvec lw(nclusters);  // log(weight)
  arma::uvec kappa(n);  // support: 0 ~ nclusters-1
  
  double lmpp;  // log marginal partition posterior
  arma::mat evalyPDFs(ndpost, ngrid);
  arma::mat evalyCDFs(ndpost, ngrid);
  arma::mat quantiles(ndpost, nprobs);
  
  
  //------------------------------------------------------------------
  // initialize hyperparameters
  if(alpha < 0) {
    alpha = arma::randg<double>(arma::distr_param(a0, (1.0/b0)));
  }
  b_gd = b_gd * alpha;
  
  if(useHyperpriors) {
    m0 = arma::mean(data, 0).t();
    
    arma::rowvec tmp_vec = arma::range(data, 0);
    tmp_vec = tmp_vec % tmp_vec / 16.0;
    S0.diag() = tmp_vec.t();
    
    m = arma::mvnrnd(m0, S0);
    
    lambda = arma::randg<double>(arma::distr_param(gamma1, (1.0/gamma2)));
    
    Psi0 = S0 / nu0;
    Psi = arma::wishrnd(Psi0, nu0);
    
  } else {
    m = arma::mean(data, 0).t();
    
    arma::rowvec tmp_vec = arma::range(data, 0);
    tmp_vec = tmp_vec % tmp_vec / 16.0;
    
    Psi.diag() = tmp_vec.t();
  }
  
  
  //------------------------------------------------------------------
  // initialize parameters
  setparam(n, nclusters, m, lambda, nu, Psi, alpha, Omega, cholOmega, Zeta, lw, a_gd, b_gd, kappa);
  
  
  //------------------------------------------------------------------
  // return data structures
  Rcpp::NumericVector lmpps(ndpost);
  
  
  //------------------------------------------------------------------
  // start mcmc
  for(arma::uword i=0; i<(nskip+ndpost); i++){
    if(i<nskip){
      // update (hyper)parameters
      drawparam(n, d, nclusters, data, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, 
                alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, lw, a_gd, b_gd, kappa, diag, lmpp);
    } else {
      // update (hyper)parameters
      for(arma::uword j=0; j<keepevery; j++){
        drawparam(n, d, nclusters, data, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, 
                  alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, lw, a_gd, b_gd, kappa, diag, lmpp);
      }
      
      // prediction
      if(pdf || cdf) {
        arma::rowvec tmp_pdf(ngrid);
        arma::rowvec tmp_cdf(ngrid);
        arma::rowvec tmp_quantile(nprobs);
        
        predict_marginal(ngrid, npred, d, nclusters, nprobs, grid, xpred, probs, Zeta, Omega, 
                         lw, pdf, cdf, diri, tmp_pdf, tmp_cdf, tmp_quantile);
        
        if(pdf) {
          evalyPDFs.row(i-nskip) = tmp_pdf;
        }
        if(cdf) {
          evalyCDFs.row(i-nskip) = tmp_cdf;
          quantiles.row(i-nskip) = tmp_quantile;
        }
      }
      
      if(diag) 
        lmpps[i-nskip] = lmpp;
    }
  }
  
  
  //------------------------------------------------------------------
  // return
  Rcpp::List res;
  if(diag)
    res["logMPPs"] = lmpps;
  
  if(pdf) {
    res["predict.pdfs"] = Rcpp::wrap(evalyPDFs);
  }
  
  if(cdf) {
    res["predict.cdfs"] = Rcpp::wrap(evalyCDFs);
    res["predict.quantiles"] = Rcpp::wrap(quantiles);
  }
  
  return res;
}

















