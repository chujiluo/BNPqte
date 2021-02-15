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
    const arma::uword nclusters,
    const arma::uword nskip,
    const arma::uword ndpost,
    const arma::uword keepevery,
    const arma::uword printevery,
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
  arma::cube evalcPDFs(npred, ngrid, ndpost);
  arma::cube evalcCDFs(npred, ngrid, ndpost);
  arma::mat evalcMeans(npred, ndpost);
  
  
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
  
  
  Rcpp::Rcout << "MCMC updating..." << std::endl;
  arma::uword nmcmc = nskip + ndpost*keepevery;
  //------------------------------------------------------------------
  // start mcmc
  for(arma::uword i=0; i<(nskip+ndpost); i++){
    if(i<nskip){
      // update (hyper)parameters
      auto start = std::chrono::high_resolution_clock::now();
      
      drawparam(n, d, nclusters, data, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, 
                alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, lw, a_gd, b_gd, kappa);
      
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
      Rcpp::Rcout << "drawparam (before ndpost): " << duration.count() << std::endl;
        
      if(((i+1)%printevery) == 0)
        Rcpp::Rcout << " - MCMC scan " << i+1 << " of " << nmcmc << std::endl;
    } else {
      // update (hyper)parameters
      for(arma::uword j=0; j<keepevery; j++){
        auto start = std::chrono::high_resolution_clock::now();
        
        drawparam(n, d, nclusters, data, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, 
                  alpha, m, lambda, nu, Psi, Omega, cholOmega, icholOmega, othersOmega, Zeta, lw, a_gd, b_gd, kappa);
        
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
        Rcpp::Rcout << "drawparam (after ndpost): " << duration.count() << std::endl;
        
        if(((nskip+(i-nskip)*keepevery+j+1)%printevery) == 0)
          Rcpp::Rcout << " - MCMC scan " << nskip+(i-nskip)*keepevery+j+1 << " of " << nmcmc << std::endl;
      }
      
      // prediction
      if(pdf || cdf || meanReg) {
        arma::mat tmp_pdf(npred, ngrid);
        arma::mat tmp_cdf(npred, ngrid);
        arma::colvec tmp_mean(npred);
        
        auto start = std::chrono::high_resolution_clock::now();
        
        predict_conditional(ngrid, npred, d, nclusters, grid, xpred, Zeta, Omega, lw, pdf, cdf, meanReg, tmp_pdf, tmp_cdf, tmp_mean);
        
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
        Rcpp::Rcout << "drawppredict_conditional: " << duration.count() << std::endl;
        
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
  res["posterior"] = posterior;
  
  if(meanReg) {
    auto start = std::chrono::high_resolution_clock::now();
    
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
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
    Rcpp::Rcout << "meanReg: " << duration.count() << std::endl;
  }
  
  if(pdf) {
    auto start = std::chrono::high_resolution_clock::now();
    
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
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
    Rcpp::Rcout << "pdf: " << duration.count() << std::endl;
  }
  
  if(cdf) {
    auto start = std::chrono::high_resolution_clock::now();
    
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
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
    Rcpp::Rcout << "cdf: " << duration.count() << std::endl;
  }
  
  return res;
}

















