#include "dpm.h"
  
// [[Rcpp::export]]
Rcpp::List cDPMcdensity (
    const arma::uword n,
    const arma::uword d,
    const arma::mat & data,   // nxd matrix
    const arma::colvec & y,   // vector of length n
    const arma::mat & x,  // nx(d-1) matrix
    const bool status,
    const bool pdf,
    const bool cdf,
    const arma::uword ngrid,
    const arma::uword npred,
    const bool updateAlpha,
    const bool useHyperpriors,
    const double a0,
    const double b0,
    const arma::colvec & m0,  // vector of length d
    const arma::mat & S0,  //dxd matrix
    const double gamma1,
    const double gamma2,
    const int nu0,
    const arma::mat & Psi0, //dxd matrix
    const int nu,
    const arma::uword nclusters,
    const arma::uword nskip,
    const arma::uword ndpost,
    const arma::uword keepevery,
    double alpha,
    double lambda,
    Rcpp::Nullable<Rcpp::NumericVector> m_ = R_NilValue,  // vector of length d
    Rcpp::Nullable<Rcpp::NumericMatrix> Psi_ = R_NilValue,  // dxd matrix
    Rcpp::Nullable<Rcpp::NumericMatrix> Zeta_ = R_NilValue,  // dxnclusters matrix
    Rcpp::Nullable<Rcpp::List> Omega_ = R_NilValue,  // nclusters elements with each is a dxd matrix
    Rcpp::Nullable<Rcpp::NumericVector> a_gd_ = R_NilValue,  // vector of length nclusters-1
    Rcpp::Nullable<Rcpp::NumericVector> b_gd_ = R_NilValue,  // vector of length nclusters-1
    Rcpp::Nullable<Rcpp::NumericVector> lw_ = R_NilValue,  // vector of length nclusters
    Rcpp::Nullable<Rcpp::IntegerVector> kappa_ = R_NilValue,  // vector of length n
    Rcpp::Nullable<Rcpp::NumericVector> grid_ = R_NilValue,  //vector of length ngrid
    Rcpp::Nullable<Rcpp::NumericMatrix> xpred_ = R_NilValue  // npred x (d-1) matrix
) {
  
  //------------------------------------------------------------------
  // process args
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
  arma::colvec grid(ngrid);
  if(grid_.isNotNull()) {
    Rcpp::NumericVector tmp(grid_);
    grid = Rcpp::as<arma::colvec>(tmp);
  }
  arma::mat xpred(npred, (d-1));
  if(xpred_.isNotNull()) {
    Rcpp::NumericMatrix tmp(xpred_);
    xpred = Rcpp::as<arma::mat>(tmp);
  }
  arma::mat evalcPDF(npred, ngrid);  // for each x_i, estimate ngrid conditional density
  arma::mat evalcCDF(npred, ngrid);  // for each x_i, estimate ngrid conditional Cdf
  
  
  //------------------------------------------------------------------
  // initialize hyperparameters
  if(alpha < 0) {
    alpha = arma::randg<double>(arma::distr_param(a0, (1.0/b0)));
  }
  b_gd = b_gd * alpha;
  if(lambda < 0) {
    lambda = arma::randg<double>(arma::distr_param(gamma1, (1.0/gamma2)));
  }
  if(m_.isNull()) {
    m = arma::mvnrnd(m0, S0);
  } else {
    Rcpp::NumericVector rcppm(m_);
    m = Rcpp::as<arma::colvec>(rcppm);
  }
  if(Psi_.isNull()) {
    Psi = arma::wishrnd(Psi0, nu0);
  } else {
    Rcpp::NumericMatrix rcppPsi(Psi_);
    Psi = Rcpp::as<arma::mat>(rcppPsi);
  }
  
  
  Rcpp::Rcout << "Initializing..." << std::endl;
  //------------------------------------------------------------------
  // initialize parameters
  if(!status) {
    Rcpp::NumericMatrix rcppZeta(Zeta_);
    Zeta = Rcpp::as<arma::mat>(rcppZeta);
    Rcpp::List rcppOmega(Omega_);
    for(arma::uword i=0; i<nclusters; i++) {
      arma::mat tmp1 = Rcpp::as<arma::mat>(rcppOmega[i]);
      Omega.slice(i) = tmp1;
      cholOmega.slice(i) = arma::chol(tmp1);
    }
    Rcpp::NumericVector rcppa_gd(a_gd_);
    a_gd = Rcpp::as<arma::rowvec>(rcppa_gd);
    Rcpp::NumericVector rcppb_gd(b_gd_);
    b_gd = Rcpp::as<arma::rowvec>(rcppb_gd);
    Rcpp::NumericVector rcpplw(lw_);
    lw = Rcpp::as<arma::rowvec>(rcpplw);
    Rcpp::IntegerVector rcppkappa(kappa_);
    kappa = Rcpp::as<arma::uvec>(rcppkappa);
  } else {
    setparam(n, nclusters, m, lambda, nu, Psi, alpha, Omega, cholOmega, Zeta, lw, a_gd, b_gd, kappa);
  }
  
  
  //------------------------------------------------------------------
  // return data structures
  Rcpp::List OmegaList(ndpost);  // each is a list of nclusters elements with each element is a dxd mat
  Rcpp::List ZetaList(ndpost);  // each is a nclustersxd mat
  Rcpp::NumericMatrix lwList(ndpost, nclusters);
  Rcpp::IntegerMatrix kappaList(ndpost, n);
  Rcpp::NumericVector alphaList(ndpost);
  Rcpp::NumericMatrix mList(ndpost, d);
  Rcpp::NumericVector lambdaList(ndpost);
  Rcpp::List PsiList(ndpost);  // each is d x d
  Rcpp::List predcPDF(ndpost);  // each is npred x ngrid
  Rcpp::List predcCDF(ndpost);  // each is npred x ngrid
  
  
  Rcpp::Rcout << "MCMC updating..." << std::endl;
  //------------------------------------------------------------------
  // start mcmc
  for(arma::uword i=0; i<(nskip+ndpost); i++){
    if(i<nskip){
      // update (hyper)parameters
      drawparam(n, d, nclusters, data, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, 
                alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, lw, a_gd, b_gd, kappa);
      
    } else {
      // update (hyper)parameters
      for(arma::uword j=0; j<keepevery; j++){
        drawparam(n, d, nclusters, data, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, 
                  alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, lw, a_gd, b_gd, kappa);
      }
      
      // prediction
      if(pdf || cdf) {
        predict_conditional(ngrid, npred, d, nclusters, grid, xpred, Zeta, Omega, lw, pdf, cdf, evalcPDF, evalcCDF);
        if(pdf)
          predcPDF[i-nskip] = Rcpp::wrap(evalcPDF);
        if(cdf)
          predcCDF[i-nskip] = Rcpp::wrap(evalcCDF);
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
      if(updateAlpha){
        alphaList[i-nskip] = alpha;
      }
      if(useHyperpriors){
        for(arma::uword j=0; j<d; j++){
          mList(i-nskip, j) = m(j);
        }
        lambdaList[i-nskip] = lambda;
        PsiList[i-nskip] = Rcpp::wrap(Psi);
      }
      
    }
  }
  
  
  Rcpp::Rcout << "Collecting returns..." << std::endl;
  //------------------------------------------------------------------
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
  
  
  //------------------------------------------------------------------
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
  if(pdf)
    res["predict.pdf"] = predcPDF;
  if(cdf)
    res["predict.cdf"] = predcCDF;
  
  return res;
}

















