#include "dpmNeal.h"

// [[Rcpp::export]]
Rcpp::List cDPMdensityNeal (
    const arma::uword n,
    const arma::uword d,
    const arma::mat & y,   //nxd matrix
    const bool status,
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
    arma::uword & nclusters, // number of USED clusters
    const arma::uword nskip,
    const arma::uword ndpost,
    const arma::uword keepevery,
    const arma::uword printevery,
    double & alpha,
    double & lambda,
    arma::colvec & m,  // vector of length d
    arma::mat & Psi,   // dxd matrix
    Rcpp::Nullable<Rcpp::NumericMatrix> Zeta_ = R_NilValue,
    Rcpp::Nullable<Rcpp::List> Omega_ = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> kappa_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> grid1_ = R_NilValue,  //vector of length ngrid
    Rcpp::Nullable<Rcpp::NumericVector> grid2_ = R_NilValue  //vector of length ngrid
) {
  //------------------------------------------------------------------
  // initialize parameters and clusters
  arma::mat invS0 = arma::inv_sympd(S0);
  arma::colvec invS0m0 = invS0 * m0;
  arma::mat invPsi0 = arma::inv_sympd(Psi0);
  
  arma::mat Zeta(d, n+2); // cluster-wise: column
  arma::cube Omega(d, d, n+2); // cluster-wise: slice
  arma::cube cholOmega(d, d, n+2); // cholOmega.slice(i) = arma::chol(Omega.slice(i))
  arma::cube icholOmega(d, d, n+2); // icholOmega.slice(i) = arma::inv(cholOmega.slice(i))
  arma::colvec othersOmega(n+2); // -(d/2)*log2pi - log(det(Omega.slice(i)))/2, terms in the log pdf of Normal(Zeta.col(i), Omega.slice(i))
  arma::uvec kappa(n);  // cluster indices
  arma::urowvec clusterSize(n+2, arma::fill::zeros);
  
  if(!status) {
    // use previous analysis
    Rcpp::NumericMatrix rcppZeta(Zeta_);
    Zeta = Rcpp::as<arma::mat>(rcppZeta);
    Rcpp::List rcppOmega(Omega_);
    for(arma::uword i=0; i<nclusters; i++) { // nclusters is given in R
      arma::mat tmp1 = Rcpp::as<arma::mat>(rcppOmega[i]);
      arma::mat tmp2 = arma::chol(tmp1);
      arma::mat tmp3 = arma::inv(arma::trimatu(tmp2));
      Omega.slice(i) = tmp1;
      cholOmega.slice(i) = tmp2;
      icholOmega.slice(i) = tmp3;
      othersOmega(i) = arma::sum(log(tmp3.diag())) - (double)d/2.0 * log2pi;
    }
    Rcpp::IntegerVector rcppkappa(kappa_);
    kappa = Rcpp::as<arma::uvec>(rcppkappa);
    for(arma::uword i=0; i<n; i++){
      clusterSize(kappa(i)) = clusterSize(kappa(i)) + 1;
    }
  } else {
    setparamNeal(n, d, nclusters, m, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, kappa, clusterSize);
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
  Rcpp::NumericVector nclusterList(ndpost);
  Rcpp::IntegerMatrix kappaList(ndpost, n);
  Rcpp::NumericVector alphaList(ndpost);
  Rcpp::NumericMatrix mList(ndpost, d);
  Rcpp::NumericVector lambdaList(ndpost);
  Rcpp::List PsiList(ndpost);  // each is dxd
  
  Rcpp::List predPDFs(ndpost);  // each is a ngridxngrid mat
  Rcpp::NumericMatrix predPDFm(ngrid, ngrid);
  
  
  //------------------------------------------------------------------
  // start mcmc
  arma::uword nmcmc = nskip + ndpost*keepevery;
  
  for(arma::uword i=0; i<(nskip+ndpost); i++){
    if(i<nskip){
      // update (hyper)parameters
      drawparamNeal(n, d, nclusters, y, updateAlpha, useHyperpriors, a0, b0, m0, S0, invS0, invS0m0, gamma1, gamma2, nu0, Psi0, invPsi0,
                    alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, kappa, clusterSize);
      if(((i+1)%printevery) == 0)
        Rcpp::Rcout << "-------MCMC scan " << i+1 << " of " << nmcmc << std::endl;
      
    } else {
      // update (hyper)parameters
      for(arma::uword j=0; j<keepevery; j++){
        drawparamNeal(n, d, nclusters, y, updateAlpha, useHyperpriors, a0, b0, m0, S0, invS0, invS0m0, gamma1, gamma2, nu0, Psi0, invPsi0,
                      alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, kappa, clusterSize);
        if(((nskip+(i-nskip)*keepevery+j+1)%printevery) == 0)
          Rcpp::Rcout << "-------MCMC scan " << nskip+(i-nskip)*keepevery+j+1 << " of " << nmcmc << std::endl;
      }
      
      // prediction
      if(prediction) {
        predict_joint_Neal(ngrid, n, d, nclusters, ypred, Zeta, icholOmega, othersOmega, alpha, m, lambda, nu, Psi,
                           clusterSize, evalPDF);
        
        evalPDFm = evalPDFm + evalPDF;
        predPDFs[i-nskip] = Rcpp::wrap(evalPDF);
      }
      
      
      // keep the posterior sample
      nclusterList[i-nskip] = nclusters;
      ZetaList[i-nskip] = Rcpp::wrap(Zeta.t());   // nclustersxd mat
      Rcpp::List tmpOmega(nclusters);
      for(arma::uword j=0; j<nclusters; j++){
        tmpOmega[j] = Rcpp::wrap(Omega.slice(j));  // dxd mat
      }
      OmegaList[i-nskip] = tmpOmega;
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
  state["method"] = "neal";
  
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
  state["kappa"] = Rcpp::as<std::vector<double>>(Rcpp::wrap(kappa));
  
  
  //------------------------------------------------------------------
  // return
  Rcpp::List res;
  
  res["method"] = "neal";
  res["updateAlpha"] = updateAlpha;
  res["useHyperpriors"] = useHyperpriors;
  res["status"] = status;
  res["state"] = state;
  
  Rcpp::List posterior;
  posterior["Zeta"] = ZetaList;
  posterior["Omega"] = OmegaList;
  posterior["nclusters"] = nclusterList;
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

















