#include "dpm.h"

// [[Rcpp::export]]
Rcpp::List cpDPMcdensity(
    const arma::uword ngrid,
    const arma::uword npred,
    arma::colvec & grid,
    arma::mat & xpred,
    const arma::uword d,
    const arma::uword nclusters,
    const arma::uword ndpost,
    const Rcpp::List & ZetaList,
    const Rcpp::List & OmegaList,
    const arma::mat & lwList,
    const bool pdf,
    const bool cdf,
    const bool meanReg,
    const bool hpd,
    const bool bci
) {
  
  //------------------------------------------------------------------
  // process args
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
    arma::mat tZeta = Rcpp::as<arma::mat>(ZetaList[i]);
    arma::mat Zeta = tZeta.t();
    
    Rcpp::List rcppOmega = OmegaList[i];
    arma::cube Omega(d, d, nclusters);
    arma::rowvec lw(nclusters);
    for(arma::uword j=0; j<nclusters; j++) {
      Omega.slice(j) = Rcpp::as<arma::mat>(rcppOmega[j]);
      lw(j) = lwList(i, j);
    }
    
    arma::mat tmp_pdf(npred, ngrid);
    arma::mat tmp_cdf(npred, ngrid);
    arma::colvec tmp_mean(npred);
    
    predict_conditional(ngrid, npred, d, nclusters, grid, xpred, Zeta, Omega, lw, pdf, cdf, meanReg, tmp_pdf, tmp_cdf, tmp_mean);
    
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
  }
  
  
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
