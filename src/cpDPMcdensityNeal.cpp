#include "dpmNeal.h"

// [[Rcpp::export]]
Rcpp::List cpDPMcdensityNeal(
    const arma::uword ngrid,
    const arma::uword npred,
    arma::colvec & grid,
    arma::mat & xpred,
    const arma::uword n,
    const arma::uword d,
    const bool updateAlpha,
    const bool useHyperpriors,
    const Rcpp::List & ZetaList,
    const Rcpp::List & OmegaList,
    const arma::mat & kappaList,
    const arma::colvec & nclusterList,
    const arma::colvec & alphaList,
    const arma::colvec & lambdaList,
    const arma::mat & mList,
    const Rcpp::List & PsiList,
    const int nu,
    const arma::uword ndpost,
    const arma::uword printevery,
    const bool pdf,
    const bool cdf,
    const bool meanReg,
    const bool hpd,
    const bool bci
) {
  //------------------------------------------------------------------
  // process args
  double alpha;
  double lambda;
  arma::colvec m;
  arma::mat Psi;
  if(!updateAlpha)
    alpha = alphaList(0);
  if(!useHyperpriors) {
    lambda = lambdaList(0);
    m = mList.row(0).t();
    Psi = Rcpp::as<arma::mat>(PsiList[0]);
  }
  
  arma::cube evalcPDFs(npred, ngrid, ndpost);
  arma::cube evalcCDFs(npred, ngrid, ndpost);
  arma::mat evalcMeans(npred, ndpost);
  
  
  //------------------------------------------------------------------
  // return data structures
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
  // start evaluation
  for(arma::uword i=0; i<ndpost; i++) {
    Rcpp::checkUserInterrupt();
    
    // nclusters
    arma::uword nclusters = nclusterList(i);
    
    // alpha, lambda, m, Psi
    if(updateAlpha)
      alpha = alphaList(i);
    if(useHyperpriors) {
      lambda = lambdaList(i);
      m = mList.row(i).t();
      Psi = Rcpp::as<arma::mat>(PsiList[i]);
    }
    
    // Zeta and Omegas
    arma::mat tZeta = Rcpp::as<arma::mat>(ZetaList[i]);
    arma::mat Zeta = tZeta.t();
    Rcpp::List rcppOmega = OmegaList[i];
    arma::cube Omega(d, d, nclusters);
    for(arma::uword j=0; j<nclusters; j++) {
      Omega.slice(j) = Rcpp::as<arma::mat>(rcppOmega[j]);
    }
    
    // clusterSize
    arma::urowvec clusterSize(n+2, arma::fill::zeros);
    for(arma::uword j=0; j<n; j++){
      clusterSize(kappaList(i, j)) = clusterSize(kappaList(i, j)) + 1;
    }
    
    // evaluation
    arma::mat tmp_pdf(npred, ngrid, arma::fill::zeros);
    arma::mat tmp_cdf(npred, ngrid, arma::fill::zeros);
    arma::colvec tmp_mean(npred, arma::fill::zeros);
    
    predict_conditional_Neal(ngrid, npred, n, d, nclusters, grid, xpred, Zeta, Omega, alpha, m, lambda, nu, Psi, clusterSize, 
                             pdf, cdf, meanReg, tmp_pdf, tmp_cdf, tmp_mean);
    
    // keep results
    if(meanReg)
      evalcMeans.col(i) = tmp_mean;
    if(pdf) {
      evalcPDFs.slice(i) = tmp_pdf;
      predcPDFs[i] = Rcpp::wrap(tmp_pdf);
    }
    if(cdf) {
      evalcCDFs.slice(i) = tmp_cdf;
      predcCDFs[i] = Rcpp::wrap(tmp_cdf);
    }
    
    if(((i+1)%printevery) == 0)
      Rcpp::Rcout << ".";
  }
  Rcpp::Rcout << std::endl;
  
  
  //------------------------------------------------------------------
  // returns
  Rcpp::List res;
  
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