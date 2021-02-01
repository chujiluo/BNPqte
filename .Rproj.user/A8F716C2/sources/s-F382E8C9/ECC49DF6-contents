#ifndef GUARD_dpm_h
#define GUARD_dpm_h

#include "common.h"
#include "random.h"

// initialize parameters
void setparam(arma::uword n, arma::uword nclusters, arma::colvec & m, double lambda, int nu, arma::mat & Psi, double alpha, arma::cube & Omega, arma::mat & Zeta, arma::rowvec & lw, arma::rowvec & a_gd, arma::rowvec & b_gd, arma::uvec & kappa){
  arma::mat R = arma::chol(arma::inv(Psi));
  rNormalInverseWishartArma(nclusters, m, lambda, nu, Psi, R, Omega, Zeta);
  rGeneralizedDirichletArma(nclusters, a_gd, b_gd, lw);
  rCat(nclusters, n, lw, kappa);
}

// update (hyper)parameters and evaluate grid points
void drawparam(const arma::mat & y, arma::uword n, arma::uword d, arma::uword nclusters, bool updateAlpha, bool useHyperpriors, double a0, double b0, const arma::colvec & m0, const arma::mat & S0, double gamma1, double gamma2, int nu0, const arma::mat & Psi0, double & alpha, arma::colvec &m, double & lambda, int nu, arma::mat & Psi, arma::cube & Omega, arma::mat & Zeta, arma::rowvec & lw, arma::rowvec & a_gd, arma::rowvec & b_gd, arma::uvec & kappa, arma::mat & yeval, arma::colvec & evalDensity, arma::uword ngrid, bool prediction){
  
  // calculate sufficient statistics
  std::vector<std::vector<unsigned>> clusterMembers(nclusters);
  std::vector<int> clusterSize(nclusters, 0);
  arma::mat clusterSumy(nclusters, d, arma::fill::zeros);
  arma::cube clusterSumyy(d, d, nclusters, arma::fill::zeros);
  
  for(arma::uword i=0; i<n; i++){
    clusterMembers[kappa(i)].push_back(i);
    clusterSize[kappa(i)] = clusterSize[kappa(i)] + 1;
    clusterSumy.row(kappa(i)) = clusterSumy.row(kappa(i)) + y.row(i);
    clusterSumyy.slice(kappa(i)) = clusterSumyy.slice(kappa(i)) + y.row(i).t() * y.row(i);
  }
  
  
  // prepare to update b_gd
  b_gd(0) = alpha + n;
  
  for(arma::uword i=0; i<nclusters; i++){
    // update Omega and Zeta
    if(clusterSize[i]==0){
      arma::mat tmp_omega(d, d);
      arma::colvec tmp_zeta(d);
      
      rNormalInverseWishartArma(m, lambda, nu, Psi, tmp_omega, tmp_zeta);
      
      Omega.slice(i) = tmp_omega;
      Zeta.col(i) = tmp_zeta;
      
    } else{
      double lambda_new = lambda + clusterSize[i];
      arma::colvec m_new = (lambda*m + clusterSumy.row(i).t()) / lambda_new;
      int nu_new = nu + clusterSize[i];
      arma::mat Psi_new = Psi + clusterSumyy.slice(i) + 
        ((lambda*clusterSize[i]/lambda_new - clusterSize[i])/pow(clusterSize[i],2))*clusterSumy.row(i).t()*clusterSumy.row(i) +
        (lambda*clusterSize[i]/lambda_new)*(m*m.t() - clusterSumy.row(i).t()*m.t()/clusterSize[i] - m*clusterSumy.row(i)/clusterSize[i]);
      
      arma::mat tmp_omega(d, d);
      arma::colvec tmp_zeta(d);
      
      rNormalInverseWishartArma(m_new, lambda_new, nu_new, Psi_new, tmp_omega, tmp_zeta); 
      
      Omega.slice(i) = tmp_omega;
      Zeta.col(i) = tmp_zeta;
      
    }
    
    // update a_gd and b_gd
    if(i<(nclusters-1)) {
      a_gd(i) = 1.0 + clusterSize[i];
      
      if(i==0){
        b_gd(i) = b_gd(i) - clusterSize[i];
      } else{
        b_gd(i) = b_gd(i-1) - clusterSize[i];
      }
    }
    
  }
  
  // update w (weight)
  rGeneralizedDirichletArma(nclusters, a_gd, b_gd, lw);
  
  // update kappa
  arma::mat loglik(n, nclusters);
  // evaluate grid points
  arma::mat evalloglik(ngrid*ngrid, nclusters);
  
  for(arma::uword j=0; j<nclusters; j++){
    arma::mat tmp_omega = Omega.slice(j);
    arma::colvec tmp_zeta = Zeta.col(j);
    arma::mat rooti;
    double other_terms;
    
    loglik.col(j) = dMvnormArma(y, n, d, tmp_zeta, tmp_omega, rooti, other_terms, true);
    
    if(prediction){
      evalloglik.col(j) = dMvnormArma(yeval, ngrid*ngrid, d, tmp_zeta, rooti, other_terms, true);
    }
    
  }
  
  for(arma::uword i=0; i<n; i++){
    arma::rowvec lw_tmp = lw + loglik.row(i);
    rCat(nclusters, lw_tmp, kappa(i));
  }
  
  if(prediction) {
    for(arma::uword i=0; i<(ngrid*ngrid); i++){
      arma::rowvec v = lw + evalloglik.row(i);
      evalDensity(i) = log_sum_exp(v, false);
    }
  }

  
  // update alpha
  if(updateAlpha){
    double a0_new = nclusters + a0 - 1.0;
    double b0_new = 1.0 / (b0 - lw(nclusters-1));
    
    alpha = arma::randg<double>(arma::distr_param(a0_new, b0_new));
  }
  
  if(useHyperpriors){
    arma::mat sumInvOmega(d, d, arma::fill::zeros);
    arma::colvec sumInvOmegaZeta(d, arma::fill::zeros);
    double sumZetaInvOmegaZeta = 0.0;
    
    for(arma::uword i=0; i<nclusters; i++){
      arma::mat tmp_mat = arma::inv_sympd(Omega.slice(i));
      arma::colvec tmp_vec = tmp_mat * Zeta.col(i);
      double tmp = arma::as_scalar(tmp_vec.t() * Zeta.col(i));

      sumInvOmega = sumInvOmega + tmp_mat;
      sumInvOmegaZeta = sumInvOmegaZeta + tmp_vec;
      sumZetaInvOmegaZeta = sumZetaInvOmegaZeta + tmp;
    }
    
    // update m
    arma::mat invS0 = arma::inv_sympd(S0);
    arma::mat m_cov = arma::inv_sympd(lambda * sumInvOmega + invS0);
    arma::colvec m_mean = m_cov * (lambda * sumInvOmegaZeta + invS0 * m0);
    m = arma::mvnrnd(m_mean, m_cov);
    
    // update lambda
    double gamma1_new = nclusters + gamma1;
    double gamma2_new = sumZetaInvOmegaZeta - 2.0 * arma::as_scalar(m.t() * sumInvOmegaZeta) + arma::as_scalar(m.t() * sumInvOmega * m);
    gamma2_new = 1.0 / gamma2 + gamma2_new / 2.0;
    lambda = arma::randg<double>(arma::distr_param(gamma1_new, gamma2_new));
    
    // update Psi
    double nu0_new = nclusters * nu + nu0;
    arma::mat Psi0_new = arma::inv_sympd(sumInvOmega + arma::inv_sympd(Psi0));
    Psi = arma::wishrnd(Psi0_new, nu0_new);
  }
  
}


#endif