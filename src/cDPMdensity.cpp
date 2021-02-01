#include "dpm.h"

// [[Rcpp::export]]
Rcpp::List cDPMdensity (
    const arma::uword n,
    const arma::uword d,
    const arma::mat & y,   //nxd matrix
    const bool prediction,
    const arma::uword ngrid,
    const arma::colvec & grid1,  //vector of length ngrid
    const arma::colvec & grid2,  //vector of length ngrid
    const bool updateAlpha,
    const bool useHyperpriors,
    const double a0,
    const double b0,
    double alpha,
    const arma::colvec & m0,  //vector of length d
    const arma::mat & S0, //dxd matrix
    arma::colvec & m, //vector of length d
    const double gamma1,
    const double gamma2,
    double lambda,
    const int nu0,
    const arma::mat & Psi0, //dxd matrix
    const int nu,
    arma::mat & Psi, //dxd matrix
    const arma::uword nclusters,
    const arma::uword nskip,
    const arma::uword ndpost,
    const arma::uword keepevery,
    Rcpp::Nullable<Rcpp::NumericMatrix> Zeta_ = R_NilValue,
    Rcpp::Nullable<Rcpp::List> Omega_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> a_gd_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> b_gd_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> lw_ = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> kappa_ = R_NilValue
) {
  
  // declare variables
  arma::mat Zeta(d, nclusters); // each column is of length d
  arma::cube Omega(d, d, nclusters); // each slice is dxd
  arma::rowvec a_gd(nclusters-1, arma::fill::ones);
  arma::rowvec b_gd(nclusters-1, arma::fill::ones);
  b_gd = b_gd * alpha;
  arma::rowvec lw(nclusters);  // log(weight)
  arma::uvec kappa(n);  // support: 0 ~ nclusters-1
  arma::mat yeval(ngrid*ngrid, d);
  if(prediction) {
    arma::uword tmp_idx = 0;
    for(arma::uword i=0; i<ngrid; i++){
      for(arma::uword j=0; j<ngrid; j++){
        yeval(tmp_idx, 0) = grid1(i);
        yeval(tmp_idx, 1) = grid2(j);
        tmp_idx = tmp_idx + 1;
      }
    }
  }
  arma::colvec evalDensity(ngrid*ngrid);
  
  // process optional arguments
  bool previous = false;
  if(Zeta_.isNotNull()) {
    previous = true;
    Rcpp::NumericMatrix rcppZeta(Zeta_);
    Zeta = Rcpp::as<arma::mat>(rcppZeta);
  }
  if(Omega_.isNotNull()) {
    previous = true;
    Rcpp::List rcppOmega(Omega_);
    for(arma::uword i=0; i<nclusters; i++) {
      Omega.slice(i) = Rcpp::as<arma::mat>(rcppOmega[i]);
    }
  }
  if(a_gd_.isNotNull()) {
    previous = true;
    Rcpp::NumericVector rcppa_gd(a_gd_);
    a_gd = Rcpp::as<arma::rowvec>(rcppa_gd);
  }
  if(b_gd_.isNotNull()) {
    previous = true;
    Rcpp::NumericVector rcppb_gd(b_gd_);
    b_gd = Rcpp::as<arma::rowvec>(rcppb_gd);
  }
  if(lw_.isNotNull()) {
    previous = true;
    Rcpp::NumericVector rcpplw(lw_);
    lw = Rcpp::as<arma::rowvec>(rcpplw);
  }
  if(kappa_.isNotNull()) {
    previous = true;
    Rcpp::IntegerVector rcppkappa(kappa_);
    kappa = Rcpp::as<arma::uvec>(rcppkappa);
  }
  
  // initialize (hyper)parameters
  if(!previous) {
    setparam(n, nclusters, m, lambda, nu, Psi, alpha, Omega, Zeta, lw, a_gd, b_gd, kappa);
  }
  
  // return data structures
  Rcpp::List OmegaList(ndpost);  // each is a list of nclusters elements with each element is a dxd mat
  Rcpp::List ZetaList(ndpost);  // each is a nclustersxd mat
  Rcpp::NumericMatrix lwList(ndpost, nclusters);
  Rcpp::IntegerMatrix kappaList(ndpost, n);
  Rcpp::NumericVector alphaList(ndpost*keepevery+nskip);
  Rcpp::NumericMatrix mList(ndpost, d);
  Rcpp::NumericVector lambdaList(ndpost);
  Rcpp::List PsiList(ndpost);  //each is d x d
  Rcpp::NumericMatrix predlDensities(ndpost, ngrid*ngrid);
  
  
  // start mcmc
  for(arma::uword i=0; i<(nskip+ndpost); i++){
    if(i<nskip){
      // update (hyper)parameters
      drawparam(y, n, d, nclusters, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, alpha, m, lambda, nu, Psi, 
                Omega, Zeta, lw, a_gd, b_gd, kappa, yeval, evalDensity, ngrid, false);
      
      if(updateAlpha){
        alphaList[i] = alpha;
      }
    } else {
      // update (hyper)parameters
      for(arma::uword j=0; j<keepevery; j++){
        if(j==(keepevery-1)){
          // evaluate grid points
          drawparam(y, n, d, nclusters, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, alpha, m, lambda, nu, Psi, 
                    Omega, Zeta, lw, a_gd, b_gd, kappa, yeval, evalDensity, ngrid, prediction);
        } else {
          drawparam(y, n, d, nclusters, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, alpha, m, lambda, nu, Psi, 
                    Omega, Zeta, lw, a_gd, b_gd, kappa, yeval, evalDensity, ngrid, false);
        }
        
        if(updateAlpha){
          alphaList[nskip+(i-nskip)*keepevery+j] = alpha;
        }
      }
      
      // keep the posterior sample
      ZetaList[i-nskip] = Rcpp::wrap(Zeta.t());   // nclustersxd mat
      
      Rcpp::List tempOmega(nclusters);
      for(arma::uword j=0; j<nclusters; j++){
        tempOmega[j] = Rcpp::wrap(Omega.slice(j));  // dxd mat
        
        lwList(i-nskip, j) = lw(j);
      }
      OmegaList[i-nskip] = tempOmega;
      
      for(arma::uword j=0; j<n; j++){
        kappaList(i-nskip, j) = kappa(j);
      }
      
      if(useHyperpriors){
        for(arma::uword j=0; j<d; j++){
          mList(i-nskip, j) = m(j);
        }
        
        lambdaList[i-nskip] = lambda;
        
        PsiList[i-nskip] = Rcpp::wrap(Psi);
      }
      
      // predict grid points
      if(prediction){
        for(arma::uword j=0; j<(ngrid*ngrid); j++){
          predlDensities(i-nskip, j) = evalDensity(j);
        }
      }
      
    }
  }
  

  // keep current state
  Rcpp::List state;
  state["updateAlpha"] = updateAlpha;
  state["useHyperpriors"] = useHyperpriors;
  state["nclusters"] = nclusters;
  state["a0"] = a0;
  state["b0"] = b0;
  state["m0"] = Rcpp::wrap(m0);
  state["S0"] = Rcpp::wrap(S0);
  state["gamma1"] = gamma1;
  state["gamma2"] = gamma2;
  state["nu"] = nu;
  state["nu0"] = nu0;
  state["Psi0"] = Rcpp::wrap(Psi0);
  state["alpha"] = alpha;
  state["m"] = Rcpp::wrap(m);
  state["lambda"] = lambda;
  state["Psi"] = Rcpp::wrap(Psi);
  state["Zeta"] = ZetaList[ndpost-1];
  state["Omega"] = OmegaList[ndpost-1];
  state["lw"] = Rcpp::wrap(lw);
  state["a_gd"] = Rcpp::wrap(a_gd);
  state["b_gd"] = Rcpp::wrap(b_gd);
  state["kappa"] = Rcpp::wrap(kappa);
    
    
  // return
  Rcpp::List res;
  res["state"] = state;
  res["cluster.mean"] = ZetaList;
  res["cluster.covariance"] = OmegaList;
  res["cluster.weight"] = lwList;
  res["individual.cluster"] = kappaList;
  if(updateAlpha) {
    res["alpha"] = alphaList;
  }
  if(useHyperpriors) {
    res["m"] = mList;
    res["lambda"] = lambdaList;
    res["Psi"] = PsiList;
  }
  if(prediction){
    res["y.pred"] = yeval;
    res["predict.densities"] = predlDensities;
  }

  return res;
}

















