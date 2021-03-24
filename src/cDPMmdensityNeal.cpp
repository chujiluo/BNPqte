#include "dpmNeal.h"

// [[Rcpp::export]]
Rcpp::List cDPMmdensityNeal(
    const arma::uword n,
    const arma::uword d,
    const arma::mat & data,   // nxd matrix
    const arma::colvec & y,   // vector of length n
    const arma::mat & x,  // nx(d-1) matrix
    const bool diag,
    const bool pdf,
    const bool cdf,
    const arma::uword ngrid,
    arma::colvec & grid,
    const arma::uword npred,
    arma::mat & xpred,
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
    arma::uword & nclusters,
    const arma::uword nskip,
    const arma::uword ndpost,
    const arma::uword keepevery,
    const arma::rowvec & diri,
    const arma::colvec & probs,
    const arma::uword nprobs
) {
  //------------------------------------------------------------------
  // initialize hyperparameters
  arma::colvec m0(d);
  arma::mat S0(d, d, arma::fill::zeros);
  arma::mat Psi0(d, d);
  arma::mat invS0(d, d);
  arma::colvec invS0m0(d);
  arma::mat invPsi0(d, d);
  arma::colvec m(d);
  arma::mat Psi(d, d);
  
  if(useHyperpriors) {
    m0 = arma::mean(data, 0).t();
    arma::rowvec tmp_vec = arma::range(data, 0);
    tmp_vec = tmp_vec % tmp_vec / 16.0;
    S0.diag() = tmp_vec.t();
    m = m0 + 100.0 * arma::randn<arma::colvec>(d);
    Psi0 = S0 / nu0;
    Psi = S0;
    
    invS0 = arma::inv_sympd(S0);
    invS0m0 = invS0 * m0;
    invPsi0 = arma::inv_sympd(Psi0);
  } else {
    m = arma::mean(data, 0).t();
    arma::rowvec tmp_vec = arma::range(data, 0);
    tmp_vec = tmp_vec % tmp_vec / 16.0;
    Psi.diag() = tmp_vec.t();
  }
  
  
  //------------------------------------------------------------------
  // initialize parameters
  arma::mat Zeta(d, nclusters); // each column is of length d
  arma::cube Omega(d, d, nclusters); // each slice is dxd
  arma::cube cholOmega(d, d, nclusters); // cholOmega.slice(i) = arma::chol(Omega.slice(i))
  arma::cube icholOmega(d, d, nclusters); // icholOmega.slice(i) = arma::inv(cholOmega.slice(i))
  arma::colvec othersOmega(nclusters); // terms excluding data in the log pdf of Normal(Zeta.col(i), Omega.slice(i))
  arma::uvec kappa(n);  // support: 0 ~ nclusters-1
  arma::urowvec clusterSize(n+2, arma::fill::zeros);
  
  double yloglik;
  
  setparamNeal(n, d, nclusters, m, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, kappa, clusterSize);
  
  
  //------------------------------------------------------------------
  // return data structures
  arma::mat evalyPDFs(ndpost, ngrid);
  arma::mat evalyCDFs(ndpost, ngrid);
  arma::mat quantiles(ndpost, nprobs);
  
  Rcpp::NumericVector ylogliks(ndpost);
  
  //------------------------------------------------------------------
  // start mcmc
  for(arma::uword i=0; i<(nskip+ndpost); i++){
    if(i<nskip){
      Rcpp::checkUserInterrupt();
      // update (hyper)parameters
      drawparamNeal(n, d, nclusters, data, updateAlpha, useHyperpriors, a0, b0, m0, S0, invS0, invS0m0, gamma1, gamma2, nu0, Psi0, invPsi0,
                    alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, kappa, clusterSize, diag, yloglik);
    } else {
      // update (hyper)parameters
      for(arma::uword j=0; j<keepevery; j++){
        Rcpp::checkUserInterrupt();
        drawparamNeal(n, d, nclusters, data, updateAlpha, useHyperpriors, a0, b0, m0, S0, invS0, invS0m0, gamma1, gamma2, nu0, Psi0, invPsi0,
                      alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, kappa, clusterSize, diag, yloglik);
      }
      
      // prediction
      if(pdf || cdf) {
        arma::rowvec tmp_pdf(ngrid);
        arma::rowvec tmp_cdf(ngrid);
        arma::rowvec tmp_quantile(nprobs);
        
        predict_marginal_Neal(ngrid, npred, n, d, nclusters, nprobs, grid, xpred, probs, Zeta, Omega, alpha, m, lambda, nu, Psi, 
                              clusterSize, pdf, cdf, diri, tmp_pdf, tmp_cdf, tmp_quantile);
        
        if(pdf) {
          evalyPDFs.row(i-nskip) = tmp_pdf;
        }
        if(cdf) {
          evalyCDFs.row(i-nskip) = tmp_cdf;
          quantiles.row(i-nskip) = tmp_quantile;
        }
      }
      
      if(diag)
        ylogliks[i-nskip] = yloglik;
    }
  }
  
  
  //------------------------------------------------------------------
  // return
  Rcpp::List res;
  
  if(diag)
    res["ylogliks"] = ylogliks;
  
  if(pdf) {
    res["predict.pdfs"] = Rcpp::wrap(evalyPDFs);
  }
  
  if(cdf) {
    res["predict.cdfs"] = Rcpp::wrap(evalyCDFs);
    res["predict.quantiles"] = Rcpp::wrap(quantiles);
  }
  
  return res;
}

















