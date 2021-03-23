#include "dpmNeal.h"

// [[Rcpp::export]]
Rcpp::List cDPMcdensityNeal(
    const arma::uword n,
    const arma::uword d,
    const arma::mat & data,   // nxd matrix
    const arma::colvec & y,   // vector of length n
    const arma::mat & x,  // nx(d-1) matrix
    const bool status,
    const bool pdf,
    const bool cdf,
    const bool meanReg,
    const arma::uword ngrid,
    const arma::uword npred,
    const bool hpd,
    const bool bci,
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
    arma::uword & nclusters,
    const arma::uword nskip,
    const arma::uword ndpost,
    const arma::uword keepevery,
    const arma::uword printevery,
    double & alpha,
    double & lambda,
    arma::colvec & m,  // vector of length d
    arma::mat & Psi,   // dxd matrix
    Rcpp::Nullable<Rcpp::NumericMatrix> Zeta_ = R_NilValue,  // dxnclusters matrix
    Rcpp::Nullable<Rcpp::List> Omega_ = R_NilValue,  // nclusters elements with each is a dxd matrix
    Rcpp::Nullable<Rcpp::IntegerVector> kappa_ = R_NilValue,  // vector of length n
    Rcpp::Nullable<Rcpp::NumericVector> grid_ = R_NilValue,  //vector of length ngrid
    Rcpp::Nullable<Rcpp::NumericMatrix> xpred_ = R_NilValue  // npred x (d-1) matrix
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
  arma::cube evalcPDFs(npred, ngrid, ndpost);
  arma::cube evalcCDFs(npred, ngrid, ndpost);
  arma::mat evalcMeans(npred, ndpost);
  
  
  //------------------------------------------------------------------
  // return data structures
  Rcpp::List OmegaList(ndpost);  // each is a list of nclusters elements with each element is a dxd mat
  Rcpp::List ZetaList(ndpost);  // each is a nclustersxd mat
  Rcpp::NumericVector nclusterList(ndpost);
  Rcpp::IntegerMatrix kappaList(ndpost, n);
  Rcpp::NumericVector alphaList(ndpost);
  Rcpp::NumericMatrix mList(ndpost, d);
  Rcpp::NumericVector lambdaList(ndpost);
  Rcpp::List PsiList(ndpost);  // each is d x d
  
  Rcpp::List predcPDFs(ndpost);
  Rcpp::NumericMatrix predcPDFl(npred, ngrid);
  Rcpp::NumericMatrix predcPDFm(npred, ngrid);
  Rcpp::NumericMatrix predcPDFh(npred, ngrid);
  
  Rcpp::List predcCDFs(ndpost);
  Rcpp::NumericMatrix predcCDFl(npred, ngrid);
  Rcpp::NumericMatrix predcCDFm(npred, ngrid);
  Rcpp::NumericMatrix predcCDFh(npred, ngrid);
  
  Rcpp::NumericMatrix predcMeans(ndpost, npred);
  Rcpp::NumericVector predcMeanl(npred);
  Rcpp::NumericVector predcMeanm(npred);
  Rcpp::NumericVector predcMeanh(npred);
  
  
  //------------------------------------------------------------------
  // start mcmc
  arma::uword nmcmc = nskip + ndpost*keepevery;
  
  for(arma::uword i=0; i<(nskip+ndpost); i++){
    if(i<nskip){
      Rcpp::checkUserInterrupt();
      
      // update (hyper)parameters
      drawparamNeal(n, d, nclusters, data, updateAlpha, useHyperpriors, a0, b0, m0, S0, invS0, invS0m0, gamma1, gamma2, nu0, Psi0, invPsi0,
                    alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, kappa, clusterSize);
      
      if(((i+1)%printevery) == 0)
        Rcpp::Rcout << "-------MCMC scan " << i+1 << " of " << nmcmc << std::endl << std::flush;
      
    } else {
      // update (hyper)parameters
      for(arma::uword j=0; j<keepevery; j++){
        Rcpp::checkUserInterrupt();
        
        drawparamNeal(n, d, nclusters, data, updateAlpha, useHyperpriors, a0, b0, m0, S0, invS0, invS0m0, gamma1, gamma2, nu0, Psi0, invPsi0,
                      alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, kappa, clusterSize);
        
        if(((nskip+(i-nskip)*keepevery+j+1)%printevery) == 0)
          Rcpp::Rcout << "-------MCMC scan " << nskip+(i-nskip)*keepevery+j+1 << " of " << nmcmc << std::endl << std::flush;
      }
      
      // prediction
      if(pdf || cdf || meanReg) {
        arma::mat tmp_pdf(npred, ngrid, arma::fill::zeros);
        arma::mat tmp_cdf(npred, ngrid, arma::fill::zeros);
        arma::colvec tmp_mean(npred, arma::fill::zeros);
        
        predict_conditional_Neal(ngrid, npred, n, d, nclusters, grid, xpred, Zeta, Omega, alpha,  m, lambda, nu, Psi, clusterSize,
                                 pdf, cdf, meanReg, tmp_pdf, tmp_cdf, tmp_mean);
        
        if(meanReg)
          evalcMeans.col(i-nskip) = tmp_mean;
        if(pdf) {
          evalcPDFs.slice(i-nskip) = tmp_pdf;
          predcPDFs[i-nskip] = Rcpp::wrap(tmp_pdf);
        }
        if(cdf) {
          evalcCDFs.slice(i-nskip) = tmp_cdf;
          predcCDFs[i-nskip] = Rcpp::wrap(tmp_cdf);
        }
      }
      
      
      // keep the posterior sample
      nclusterList[i-nskip] = nclusters;
      
      arma::mat tmpZeta = Zeta.cols(0, (nclusters-1)).t();
      ZetaList[i-nskip] = Rcpp::wrap(tmpZeta);   // nclustersxd mat
      
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
  
  if(meanReg) {
    predcMeans = Rcpp::wrap(evalcMeans.t());  //ndpostxnpred
    res["predict.meanRegs"] = predcMeans;
    
    arma::vec tmp_avg = arma::mean(evalcMeans, 1);
    predcMeanm = Rcpp::as<std::vector<double>>(Rcpp::wrap(tmp_avg));
    res["predict.meanReg.avg"] = predcMeanm;
    
    if(hpd || bci) {
      for(arma::uword i=0; i<npred; i++) {
        arma::vec lower(2);
        arma::vec upper(2);
        arma::rowvec tmp_x =  evalcMeans.row(i);
        
        credible_interval(ndpost, 0.05, tmp_x, lower, upper);
        
        if(hpd) {
          predcMeanl[i] = lower(1);
          predcMeanh[i] = upper(1);
        }
        if(bci) {
          predcMeanl[i] = lower(0);
          predcMeanh[i] = upper(0);
        }
      }
      
      res["predict.meanReg.lower"] = predcMeanl;
      res["predict.meanReg.upper"] = predcMeanh;
    }
  }
  
  if(pdf) {
    res["predict.pdfs"] = predcPDFs;
    
    arma::mat tmp_avg = arma::mean(evalcPDFs, 2);
    predcPDFm = Rcpp::wrap(tmp_avg);
    res["predict.pdf.avg"] = predcPDFm;
    
    if(hpd || bci) {
      for(arma::uword i=0; i<npred; i++) {
        for(arma::uword j=0; j<ngrid; j++) {
          arma::vec lower(2);
          arma::vec upper(2);
          arma::rowvec tmp_x = evalcPDFs.tube(i, j);
          
          credible_interval(ndpost, 0.05, tmp_x, lower, upper);
          
          if(hpd) {
            predcPDFl(i, j) = lower(1);
            predcPDFh(i, j) = upper(1);
          }
          if(bci) {
            predcPDFl(i, j) = lower(0);
            predcPDFh(i, j) = upper(0);
          }
        }
      }
      
      res["predict.pdf.lower"] = predcPDFl;
      res["predict.pdf.upper"] = predcPDFh;
    }
  }
  
  if(cdf) {
    res["predict.cdfs"] = predcCDFs;
    
    arma::mat tmp_avg = arma::mean(evalcCDFs, 2);
    predcCDFm = Rcpp::wrap(tmp_avg);
    res["predict.cdf.avg"] = predcCDFm;
    
    if(hpd || bci) {
      for(arma::uword i=0; i<npred; i++) {
        for(arma::uword j=0; j<ngrid; j++) {
          arma::vec lower(2);
          arma::vec upper(2);
          arma::rowvec tmp_x = evalcCDFs.tube(i, j); 
          
          credible_interval(ndpost, 0.05, tmp_x, lower, upper);
          
          if(hpd) {
            predcCDFl(i, j) = lower(1);
            predcCDFh(i, j) = upper(1);
          }
          if(bci) {
            predcCDFl(i, j) = lower(0);
            predcCDFh(i, j) = upper(0);
          }
        }
      }
      
      res["predict.cdf.lower"] = predcCDFl;
      res["predict.cdf.upper"] = predcCDFh;
    }
  }
  
  return res;
}

















