#include "dpm.h"

// [[Rcpp::export]]
Rcpp::List cpDPMdensity(
    const arma::uword ngrid,
    const arma::colvec & grid1,
    const arma::colvec & grid2,
    const arma::uword d,
    const arma::uword nclusters,
    const arma::uword ndpost,
    const Rcpp::List & ZetaList,
    const Rcpp::List & OmegaList,
    const arma::mat & lwList,
    const arma::uword printevery
) {
  
  //------------------------------------------------------------------
  // process args
  arma::mat ypred(ngrid*ngrid, d);  // arma::mat in column major order
  arma::uword tmp_idx = 0;
  for(arma::uword i=0; i<ngrid; i++){
    for(arma::uword j=0; j<ngrid; j++){
      ypred(tmp_idx, 0) = grid1(j);
      ypred(tmp_idx, 1) = grid2(i);
      tmp_idx = tmp_idx + 1;
    }
  }
  arma::mat evalPDF(ngrid, ngrid);
  arma::mat evalPDFm(ngrid, ngrid, arma::fill::zeros);
  
  
  //------------------------------------------------------------------
  // return data structures
  Rcpp::List predPDFs(ndpost);  // each is a ngridxngrid mat
  Rcpp::NumericMatrix predPDFm(ngrid, ngrid);
  
  
  //------------------------------------------------------------------
  // start evaluation
  for(arma::uword i=0; i<ndpost; i++) {
    arma::mat tZeta = Rcpp::as<arma::mat>(ZetaList[i]);
    arma::mat Zeta = tZeta.t();
    
    Rcpp::List tmp_omegalist = OmegaList[i];
    arma::cube icholOmega(d, d, nclusters);
    arma::colvec othersOmega(nclusters);
    arma::rowvec lw(nclusters);
    
    for(arma::uword j=0; j<nclusters; j++) {
      arma::mat tmp_omega = Rcpp::as<arma::mat>(tmp_omegalist[j]);
      arma::mat tmp_icholomega = arma::inv(arma::trimatu(arma::chol(tmp_omega)));
      double rootisum = arma::sum(log(tmp_icholomega.diag())), constants = -(double)d/2.0 * log2pi;
      
      icholOmega.slice(j) = tmp_icholomega;
      othersOmega(j) = rootisum + constants;
      
      lw(j) = lwList(i, j);
    }
    
    predict_joint(ngrid, d, nclusters, ypred, Zeta, icholOmega, othersOmega, lw, evalPDF);
    
    evalPDFm = evalPDFm + evalPDF;
    predPDFs[i] = Rcpp::wrap(evalPDF);
    
    if(((i+1)%printevery) == 0)
      Rcpp::Rcout << ".";
  }
  Rcpp::Rcout << std::endl;
  
  //------------------------------------------------------------------
  // returns
  Rcpp::List res;
  
  res["prediction"] = true;
  
  res["predict.pdfs"] = predPDFs;
  
  evalPDFm = evalPDFm / ndpost;
  predPDFm = Rcpp::wrap(evalPDFm);
  res["predict.pdf.avg"] = predPDFm;
  
  res["grid1"] = Rcpp::as<std::vector<double>>(Rcpp::wrap(grid1));
  res["grid2"] = Rcpp::as<std::vector<double>>(Rcpp::wrap(grid2));
  
  return res;
}
