#ifndef GUARD_random_h
#define GUARD_random_h

#include "common.h"

// return log(sum(exp(v)))
double log_sum_exp(arma::rowvec & v, const bool logd = true){
  double mx = v.max();
  double sm=0.;
  for(arma::uword i=0; i<v.n_elem; i++){
    sm += exp(v(i) - mx);
  }
  if(logd){
    return mx+log(sm);
  }
  return exp(mx)*sm;
}

// Multivariate distributions

// Generate a sample from Normal-Inverse-Wishart distribution
void rNormalInverseWishartArma(arma::colvec & m, double lambda, int nu, arma::mat & Psi, arma::mat & omega, arma::colvec & zeta){
  arma::mat tmp1 = arma::iwishrnd(Psi, nu);
  arma::colvec tmp2 = arma::mvnrnd(m, tmp1 / lambda);
  omega = tmp1;
  zeta = tmp2;
}
// Generate k i.i.d samples from Normal-Inverse-Wishart distribution
void rNormalInverseWishartArma(arma::uword k, arma::colvec & m, double lambda, int nu, arma::mat & Psi, arma::mat & R, arma::cube & Omega, arma::mat & Zeta){
  for(arma::uword i=0; i<k; i++){
    arma::mat tmp1 = arma::iwishrnd(Psi, nu, R);
    arma::colvec tmp2 = arma::mvnrnd(m, tmp1 / lambda);
    Omega.slice(i) = tmp1;
    Zeta.col(i) = tmp2;
  }
}

// Generate a sample from generalized Dirichlet distribution in log scale
void rGeneralizedDirichletArma(arma::uword k, arma::rowvec & a, arma::rowvec & b, arma::rowvec & lw) {
  arma::rowvec tmp(k-1, arma::fill::zeros);
  
  for(arma::uword i=0; i<(k-1); i++){
    double v1 = arma::randg<double>(arma::distr_param(a(i), 1.0));
    double v2 = arma::randg<double>(arma::distr_param(b(i), 1.0));
    if(i==0){
      lw(i) = log(v1) - log(v1+v2);
    }else{
      lw(i) = arma::sum(tmp) + log(v1) - log(v1+v2);
    }
    tmp(i) = log(v2) - log(v1+v2);
  }
  lw(k-1) = arma::sum(tmp);
}

// Generate a sample from discrete distribution with support: 0 ~ k-1; lw is unnormalized log weights
// Gumbel-max trick
void rCat(arma::uword k, arma::rowvec & lw, arma::uword kappai){
  arma::rowvec u = arma::randu<arma::rowvec>(k);  // k samples from U[0,1]
  arma::rowvec g = -log(-log(u)) + lw;
  kappai = g.index_max();
}
// Generate n samples from discrete distribution with support: 0 ~ k-1; lw is unnormalized log weights
void rCat(arma::uword k, arma::uword n, arma::rowvec & lw, arma::uvec & kappa){
  for(arma::uword i=0; i<n; i++){
    arma::rowvec u = arma::randu<arma::rowvec>(k);  // k samples from U[0,1]
    arma::rowvec g = -log(-log(u)) + lw;
    kappa(i) = g.index_max();
  }
}

// Calculate probability densities of multivariate normal
void inplace_tri_mat_mult(arma::rowvec & y, arma::uword d, arma::mat & trimat){
  for(unsigned j=d; j-- > 0;){
    double tmp(0.);
    for(unsigned i=0; i<=j; ++i)
      tmp += trimat.at(i, j) * y[i];
    y[j] = tmp;
  }
}

arma::colvec dMvnormArma(const arma::mat & y, arma::uword n, arma::uword d, arma::colvec & zeta, arma::mat & omega, arma::mat & rooti, double & other_terms, const bool logd = true) { 
  arma::colvec out(n);
  
  rooti = arma::inv(arma::trimatu(arma::chol(omega)));
  double rootisum = arma::sum(log(rooti.diag())), constants = -(double)d/2.0 * log2pi;
  other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (arma::uword i=0; i<n; i++) {
    z = y.row(i) - zeta.t();
    inplace_tri_mat_mult(z, d, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}

arma::colvec dMvnormArma(const arma::mat & y, arma::uword n, arma::uword d, arma::colvec & zeta, arma::mat & rooti, double & other_terms, const bool logd = true) { 
  arma::colvec out(n);
  
  arma::rowvec z;
  for (arma::uword i=0; i<n; i++) {
    z = y.row(i) - zeta.t();
    inplace_tri_mat_mult(z, d, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}




#endif