#include "dpm.h"

// [[Rcpp::export]]
Rcpp::List cDPMdensity (
    const arma::uword n,
    const arma::uword d,
    const arma::mat & y,   //nxd matrix
    const bool status,
    const bool diag,
    const bool prediction,
    const arma::uword ngrid,
    const bool updateAlpha,
    const bool useHyperpriors,
    const double a0,
    const double b0,
    const arma::colvec & m0,  //vector of length d
    const arma::mat & S0, //dxd matrix
    const double gamma1,
    const double gamma2,
    const int nu0,
    const arma::mat & Psi0, //dxd matrix
    const int nu,
    const arma::uword nclusters,
    const arma::uword nskip,
    const arma::uword ndpost,
    const arma::uword keepevery,
    const arma::uword printevery,
    double & alpha,
    double & lambda,
    arma::colvec & m,
    arma::mat & Psi,
    arma::rowvec & a_gd,
    arma::rowvec & b_gd,
    Rcpp::Nullable<Rcpp::NumericMatrix> Zeta_ = R_NilValue,
    Rcpp::Nullable<Rcpp::List> Omega_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> lw_ = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> kappa_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> grid1_ = R_NilValue,  //vector of length ngrid
    Rcpp::Nullable<Rcpp::NumericVector> grid2_ = R_NilValue  //vector of length ngrid
) {
  //------------------------------------------------------------------
  // initialize parameters and clusters
  arma::mat invS0 = arma::inv_sympd(S0);
  arma::colvec invS0m0 = invS0 * m0;
  arma::mat invPsi0 = arma::inv_sympd(Psi0);
  
  arma::mat Zeta(d, nclusters); // each column is of length d
  arma::cube Omega(d, d, nclusters); // each slice is dxd
  arma::cube cholOmega(d, d, nclusters); // cholOmega.slice(i) = arma::chol(Omega.slice(i))
  arma::cube icholOmega(d, d, nclusters); // icholOmega.slice(i) = arma::inv(cholOmega.slice(i))
  arma::colvec othersOmega(nclusters); // terms excluding data in the log pdf of Normal(Zeta.col(i), Omega.slice(i))
  arma::rowvec lw(nclusters);  // log(weight)
  arma::uvec kappa(n);  // support: 0 ~ nclusters-1
  
  double lmpp;  // log marginal partition posterior
  
  if(!status) {
    Rcpp::NumericMatrix rcppZeta(Zeta_);
    Zeta = Rcpp::as<arma::mat>(rcppZeta);
    Rcpp::List rcppOmega(Omega_);
    for(arma::uword i=0; i<nclusters; i++) {
      arma::mat tmp1 = Rcpp::as<arma::mat>(rcppOmega[i]);
      Omega.slice(i) = tmp1;
      cholOmega.slice(i) = arma::chol(tmp1);
    }
    Rcpp::NumericVector rcpplw(lw_);
    lw = Rcpp::as<arma::rowvec>(rcpplw);
    Rcpp::IntegerVector rcppkappa(kappa_);
    kappa = Rcpp::as<arma::uvec>(rcppkappa);
  } else {
    setparam(n, nclusters, m, lambda, nu, Psi, alpha, Omega, cholOmega, Zeta, lw, a_gd, b_gd, kappa);
  }
  
  
  //------------------------------------------------------------------
  // process prediction args
  arma::colvec grid1(ngrid);
  arma::colvec grid2(ngrid);
  arma::mat ypred(ngrid*ngrid, d);  // arma::mat in column major order
  if(grid1_.isNotNull()) {
    Rcpp::NumericVector tmp1(grid1_);
    Rcpp::NumericVector tmp2(grid2_);
    grid1 = Rcpp::as<arma::colvec>(tmp1);
    grid2 = Rcpp::as<arma::colvec>(tmp2);
    arma::uword tmp_idx = 0;
    for(arma::uword i=0; i<ngrid; i++){
      for(arma::uword j=0; j<ngrid; j++){
        ypred(tmp_idx, 0) = grid1(j);
        ypred(tmp_idx, 1) = grid2(i);
        tmp_idx = tmp_idx + 1;
      }
    }
  }
  arma::mat evalPDF(ngrid, ngrid);
  arma::mat evalPDFm(ngrid, ngrid, arma::fill::zeros);
  

  //------------------------------------------------------------------
  // return data structures
  Rcpp::List OmegaList(ndpost);  // each is a list of nclusters elements with each element is a dxd mat
  Rcpp::List ZetaList(ndpost);  // each is a nclustersxd mat
  Rcpp::NumericMatrix lwList(ndpost, nclusters);
  Rcpp::IntegerMatrix kappaList(ndpost, n);
  Rcpp::NumericVector alphaList(ndpost);
  Rcpp::NumericMatrix mList(ndpost, d);
  Rcpp::NumericVector lambdaList(ndpost);
  Rcpp::List PsiList(ndpost);  // each is dxd
  Rcpp::NumericVector lmpps(ndpost);  // log marginal partition posteriors
  
  Rcpp::List predPDFs(ndpost);  // each is a ngridxngrid mat
  Rcpp::NumericMatrix predPDFm(ngrid, ngrid);
  
  
  //------------------------------------------------------------------
  // start mcmc
  arma::uword nmcmc = nskip + ndpost*keepevery;
  
  for(arma::uword i=0; i<(nskip+ndpost); i++){
    if(i<nskip){
      // update (hyper)parameters
      drawparam(n, d, nclusters, y, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, 
                alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, lw, a_gd, b_gd, kappa, diag, lmpp);
      if(((i+1)%printevery) == 0)
        Rcpp::Rcout << "-------MCMC scan " << i+1 << " of " << nmcmc << std::endl;
      
    } else {
      // update (hyper)parameters
      for(arma::uword j=0; j<keepevery; j++){
        drawparam(n, d, nclusters, y, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, 
                  alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, lw, a_gd, b_gd, kappa, diag, lmpp);
        if(((nskip+(i-nskip)*keepevery+j+1)%printevery) == 0)
          Rcpp::Rcout << "-------MCMC scan " << nskip+(i-nskip)*keepevery+j+1 << " of " << nmcmc << std::endl;
      }
      
      // prediction
      if(prediction) {
        predict_joint(ngrid, d, nclusters, ypred, Zeta, icholOmega, othersOmega, lw, evalPDF);
        
        evalPDFm = evalPDFm + evalPDF;
        predPDFs[i-nskip] = Rcpp::wrap(evalPDF);
      }
      
      // log marginal partition posterior
      if(diag) 
        lmpps[i-nskip] = lmpp;
      
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
  

  Rcpp::Rcout << "*****Collecting returns..." << std::endl;
  //------------------------------------------------------------------
  // keep current state
  Rcpp::List state;
  state["method"] = "truncated";
  
  state["updateAlpha"] = updateAlpha;
  state["useHyperpriors"] = useHyperpriors;
  state["nclusters"] = nclusters;
  state["a0"] = a0;
  state["b0"] = b0;
  state["m0"] = Rcpp::as<std::vector<double>>(Rcpp::wrap(m0));
  state["S0"] = Rcpp::wrap(S0);
  state["gamma1"] = gamma1;
  state["gamma2"] = gamma2;
  state["nu"] = nu;
  state["nu0"] = nu0;
  state["Psi0"] = Rcpp::wrap(Psi0);
  state["alpha"] = alpha;
  state["m"] = Rcpp::as<std::vector<double>>(Rcpp::wrap(m));
  state["lambda"] = lambda;
  state["Psi"] = Rcpp::wrap(Psi);
  state["Zeta"] = ZetaList[ndpost-1];
  state["Omega"] = OmegaList[ndpost-1];
  state["lw"] = Rcpp::as<std::vector<double>>(Rcpp::wrap(lw));
  state["a_gd"] = Rcpp::as<std::vector<double>>(Rcpp::wrap(a_gd));
  state["b_gd"] = Rcpp::as<std::vector<double>>(Rcpp::wrap(b_gd));
  state["kappa"] = Rcpp::as<std::vector<double>>(Rcpp::wrap(kappa));
    
  
  //------------------------------------------------------------------
  // return
  Rcpp::List res;
  
  res["method"] = "truncated";
  res["updateAlpha"] = updateAlpha;
  res["useHyperpriors"] = useHyperpriors;
  res["status"] = status;
  res["state"] = state;
  
  Rcpp::List posterior;
  posterior["Zeta"] = ZetaList;
  posterior["Omega"] = OmegaList;
  posterior["lw"] = lwList;
  posterior["kappa"] = kappaList;
  if(updateAlpha) {
    posterior["alpha"] = alphaList;
  } else {
    posterior["alpha"] = alpha;
  }
  if(useHyperpriors) {
    posterior["m"] = mList;
    posterior["lambda"] = lambdaList;
    posterior["Psi"] = PsiList;
  } else {
    posterior["m"] = Rcpp::as<std::vector<double>>(Rcpp::wrap(m));
    posterior["lambda"] = lambda;
    posterior["Psi"] = Rcpp::wrap(Psi);
  }
  if(diag)
    posterior["logMPPs"] = lmpps;
  res["posterior"] = posterior;
  
  res["prediction"] = prediction;
  if(prediction){
    res["predict.pdfs"] = predPDFs;
    
    evalPDFm = evalPDFm / ndpost;
    predPDFm = Rcpp::wrap(evalPDFm);
    res["predict.pdf.avg"] = predPDFm;
    
    res["grid1"] = Rcpp::as<std::vector<double>>(Rcpp::wrap(grid1));
    res["grid2"] = Rcpp::as<std::vector<double>>(Rcpp::wrap(grid2));
  }

  return res;
}

















