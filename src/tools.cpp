#include "common.h"

// generate n samples from dirichlet distribution
// [[Rcpp::export]]
arma::mat rdirichlet(arma::uword n, arma::colvec & probs) {
  arma::uword k = probs.n_elem;
  arma::mat y(n, k);
  arma::colvec rowsum(n, arma::fill::zeros);
  
  for(arma::uword j=0; j<k; j++) {
    y.col(j) = arma::randg(n, arma::distr_param(probs(j), 1.0));
    rowsum = rowsum + y.col(j);
  }
  for(arma::uword i=0; i<n; i++) {
    y.row(i) = y.row(i) / rowsum(i);
  }
  
  return y;
}


// calculate credible intervals of x
// [[Rcpp::export]]
arma::mat credible_interval(arma::uword n, arma::colvec & x, arma::colvec & alphas, std::string type) {
  
  arma::uword nalphas = alphas.n_elem;
  
  // return
  arma::mat band(2, nalphas);
  
  // sort x from low to high
  arma::colvec sorted_x = arma::sort(x);
  
  for(arma::uword i=0; i<nalphas; i++) {
    int nq1 = (int)round((n*alphas(i)/2.0));
    int nq2 = (int)round((n*(1.0-alphas(i)/2.0)));
    int nq = nq2 - nq1;
    
    if(type == "BCI") {
      // Bayesian credible interval
      band(0, i) = sorted_x(std::max(nq1-1, 0));
      band(1, i) = sorted_x(std::max(nq2-1, 0));
    } else {
      // Highest posterior interval
      double hpd_width = 0.0;
      double hpd_lower = 0.0;
      double hpd_upper = 0.0;
      
      for(arma::uword j=0; j<(n-nq); j++) {
        double tmp_lower = sorted_x(j);
        double tmp_upper = sorted_x(j+nq);
        double tmp_width = tmp_upper - tmp_lower;
        if(j==0) {
          hpd_width = tmp_width;
          hpd_lower = tmp_lower;
          hpd_upper = tmp_upper;
        } else {
          if(hpd_width > tmp_width) {
            hpd_width = tmp_width;
            hpd_lower = tmp_lower;
            hpd_upper = tmp_upper;
          }
        }
      }
      
      band(0, i) = hpd_lower;
      band(1, i) = hpd_upper;
    }
  }
  
  return band;
}
