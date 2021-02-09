#ifndef GUARD_random_h
#define GUARD_random_h

#include "common.h"

// return log(sum(exp(v)))
static double log_sum_exp(arma::rowvec & v, const bool logd = true){
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

// Generate a sample from a Multivariate Normal distribution: N(mu, sigma), cholsigma = chol(sigma)
static arma::colvec rMvnormArma(arma::colvec & mu, arma::mat & cholsigma) {
  int d = mu.n_elem;
  arma::colvec z = arma::randn<arma::colvec>(d);
  return mu + cholsigma.t() * z;
}

// Generate a sample from a Normal-Inverse-Wishart distribution: N(Zeta.col(i) | m, Omega.slice(i)/lambda)xIW(Omega.slice(i) | nu, Psi)
static void rNormalInverseWishartArma(arma::colvec & m, double lambda, int nu, arma::mat & Psi, arma::mat & omega, arma::mat & cholomega, arma::colvec & zeta){
  omega = arma::iwishrnd(Psi, nu);
  cholomega = arma::chol(omega);
  arma::mat tmp = cholomega / sqrt(lambda);
  zeta = rMvnormArma(m, tmp);
}

// Generate k samples from a Normal-Inverse-Wishart distribution: N(Zeta.col(i) | m, Omega.slice(i)/lambda)xIW(Omega.slice(i) | nu, Psi), R = chol(inv(Psi))
static void rNormalInverseWishartArma(arma::uword k, arma::colvec & m, double lambda, int nu, arma::mat & Psi, arma::mat & R, arma::cube & Omega, arma::cube & cholOmega, arma::mat & Zeta){
  for(arma::uword i=0; i<k; i++){
    arma::mat tmp1 = arma::iwishrnd(Psi, nu, R);
    Omega.slice(i) = tmp1;
    
    arma::mat tmp2 = arma::chol(tmp1);
    cholOmega.slice(i) = tmp2;
    
    arma::mat tmp3 = tmp2 / sqrt(lambda);
    Zeta.col(i) = rMvnormArma(m, tmp3);
  }
}

// Generate a sample from a Generalized Dirichlet distribution in log scale: (w_1, ..., w_k) ~ GD(a, b)
static void rGeneralizedDirichletArma(arma::uword k, arma::rowvec & a, arma::rowvec & b, arma::rowvec & lw) {
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

// Gumbel-max trick
// Generate a sample from a discrete distribution with support: 0 ~ k-1, lw is unnormalized log weights
static void rCat(arma::uword k, arma::rowvec & lw, arma::uword i, arma::uvec & kappa){
  arma::rowvec u = arma::randu<arma::rowvec>(k);  // k samples from U[0,1]
  arma::rowvec g = -log(-log(u)) + lw;
  kappa(i) = g.index_max();
}

// Gumbel-max trick
// Generate n samples from discrete distribution with support: 0 ~ k-1， lw is unnormalized log weights
static void rCat(arma::uword k, arma::uword n, arma::rowvec & lw, arma::uvec & kappa){
  for(arma::uword i=0; i<n; i++){
    arma::rowvec u = arma::randu<arma::rowvec>(k);  // k samples from U[0,1]
    arma::rowvec g = -log(-log(u)) + lw;
    kappa(i) = g.index_max();
  }
}

// Calculate probability densities of multivariate normal
static void inplace_tri_mat_mult(arma::rowvec & y, arma::uword d, arma::mat & trimat){
  for(unsigned j=d; j-- > 0;){
    double tmp(0.);
    for(unsigned i=0; i<=j; ++i)
      tmp += trimat.at(i, j) * y[i];
    y[j] = tmp;
  }
}

static arma::colvec dMvnormArma(const arma::mat & y, arma::uword n, arma::uword d, arma::colvec & zeta, arma::mat & cholomega, arma::mat & rooti, double & other_terms, const bool logd = true) { 
  arma::colvec out(n);
  
  rooti = arma::inv(arma::trimatu(cholomega));
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

static arma::colvec dMvnormArma(const arma::mat & y, arma::uword n, arma::uword d, arma::colvec & zeta, arma::mat & rooti, double & other_terms, const bool logd = true) { 
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