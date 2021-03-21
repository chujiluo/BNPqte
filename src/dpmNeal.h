#ifndef GUARD_dpmNeal_h
#define GUARD_dpmNeal_h

#include "common.h"
#include "random.h"

//------------------------------------------------------------------
// initialize parameters (Algorithm 8 with m = 1 in Neal, 2000)
static void setparamNeal(arma::uword n, arma::uword d, arma::uword & nclusters, arma::colvec & m, arma::mat & Psi, arma::cube & Omega, arma::cube & cholOmega, arma::cube & icholOmega, arma::colvec & othersOmega, arma::mat & Zeta, arma::uvec & kappa, arma::urowvec & clusterSize){
  
  // initialize clusters
  kappa.zeros();
  nclusters = 1;
  clusterSize(0) = n;
  
  // initialize Zeta, Omega, cholOmega, icholOmega, othersOmega
  Zeta.col(0) = m;
  Omega.slice(0) = Psi;
  arma::mat tmp1 = arma::chol(Psi);
  cholOmega.slice(0) = tmp1;
  arma::mat tmp2 = arma::inv(arma::trimatu(tmp1));
  icholOmega.slice(0) = tmp2;
  othersOmega(0) = arma::sum(log(tmp2.diag())) - (double)d/2.0 * log2pi;
}


//------------------------------------------------------------------
// update (hyper)parameters (Algorithm 8 with m = 1 in Neal, 2000)
static void drawparamNeal(arma::uword n, arma::uword d, arma::uword & nclusters, const arma::mat & y, bool updateAlpha, bool useHyperpriors, double a0, double b0, const arma::colvec & m0, const arma::mat & S0, const arma::mat & invS0, const arma::colvec & invS0m0, double gamma1, double gamma2, int nu0, const arma::mat & Psi0, const arma::mat & invPsi0, double & alpha, arma::colvec & m, double & lambda, int nu, arma::mat & Psi, arma::cube & Omega, arma::cube & cholOmega, arma::cube & icholOmega, arma::colvec & othersOmega, arma::mat & Zeta, arma::uvec & kappa, arma::urowvec & clusterSize){

  // update kappa, nclusters, clusterSize
  for(arma::uword i=0; i<n; i++) {
    
    arma::uword kMinus;
    
    if (clusterSize(kappa(i)) > 1) {
      // kappa_i = some kappa_j, for j != i
      kMinus = nclusters;
      
      arma::mat tmp_omega(d, d);
      arma::mat tmp_cholomega(d, d);
      arma::colvec tmp_zeta(d);
      
      // draw auxiliary variables
      rNormalInverseWishartArma(m, lambda, nu, Psi, tmp_omega, tmp_cholomega, tmp_zeta);
      
      arma::mat tmp_icholomega = arma::inv(arma::trimatu(tmp_cholomega));
      double tmp_othersomega = arma::sum(log(tmp_icholomega.diag())) - (double)d/2.0 * log2pi;
      
      Omega.slice(kMinus) = tmp_omega;
      cholOmega.slice(kMinus) = tmp_cholomega;
      icholOmega.slice(kMinus) = tmp_icholomega;
      othersOmega(kMinus) = tmp_othersomega;
      Zeta.col(kMinus) = tmp_zeta;
      
    } else {
      // kappa_i != any kappa_j, for j != i
      kMinus = nclusters - 1;

      // swap cluster kappa(i) and the last cluster nclusters-1
      std::swap(clusterSize(kappa(i)), clusterSize(nclusters-1));
      Omega.slice(kappa(i)).swap(Omega.slice(nclusters-1));
      cholOmega.slice(kappa(i)).swap(cholOmega.slice(nclusters-1));
      icholOmega.slice(kappa(i)).swap(icholOmega.slice(nclusters-1));
      std::swap(othersOmega(kappa(i)), othersOmega(nclusters-1));
      Zeta.swap_cols(kappa(i), nclusters-1);
      for (arma::uword j=0; j<n; j++){
        if (kappa(j) == (nclusters-1)){
          kappa(j) = kappa(i);
        }
      }
      kappa(i) = nclusters - 1;

      nclusters--;
      
      /* 
      // relabel as DPpackage (not efficient)
      arma::uword since = kappa(i);
      
      arma::colvec muwork = Zeta.col(since);
      arma::mat sigmawork = Omega.slice(since);
      arma::mat cholOmegaWork = cholOmega.slice(since);
      arma::mat icholOmegaWork = icholOmega.slice(since);
      double othersOmegaWork = othersOmega(since);
      
      if (since < (nclusters - 1)) {
        for (arma::uword iii=since+1; iii<nclusters; iii++){
          for (arma::uword j=0; j < n; j++){
            if (kappa(j) == iii){
              kappa(j) = iii - 1;
            }  
          }
          Zeta.col(iii-1) = Zeta.col(iii);
          Omega.slice(iii-1) = Omega.slice(iii);
          cholOmega.slice(iii-1) = cholOmega.slice(iii);
          icholOmega.slice(iii-1) = icholOmega.slice(iii);
          othersOmega(iii-1) = othersOmega(iii);
          
          clusterSize(iii-1) = clusterSize(iii);
        }
        kappa(i) = nclusters-1;
        
        Zeta.col(nclusters-1) = muwork;
        Omega.slice(nclusters-1) = sigmawork;
        cholOmega.slice(nclusters-1) = cholOmegaWork;
        icholOmega.slice(nclusters-1) = icholOmegaWork;
        othersOmega(nclusters-1) = othersOmegaWork;
      
        clusterSize(nclusters-1) = 1;
        
      }
      nclusters--;
       */
    }
    
    clusterSize(kappa(i))--;
    
    // calculate the probs for sampling a new kappa_i
    arma::rowvec log_weights(kMinus+1);
    arma::rowvec y_i = y.row(i);
    
    for (arma::uword k=0; k<=kMinus; k++){
      arma::colvec tmp_zeta = Zeta.col(k);
      arma::mat tmp_icholomega = icholOmega.slice(k);
      double tmp_othersomega = othersOmega(k);
      
      if(k < kMinus) {
        log_weights(k) = log(clusterSize(k)) + dMvnormArma(y_i, d, tmp_zeta, tmp_icholomega, tmp_othersomega, true);
      } else {
        log_weights(k) = log(alpha) + dMvnormArma(y_i, d, tmp_zeta, tmp_icholomega, tmp_othersomega, true);
      }
    }
    
    // sampling a new kappa_i from 0, ..., kMinus=nclusters
    arma::rowvec work_weights = exp(log_weights);
    //kappa(i) = simdisc(work_weights, (kMinus + 1), (kMinus + 1));
    rCat((kMinus + 1), log_weights, i, kappa);
    
    // update nclusters and clusterSize
    if (kappa(i) == nclusters){
      nclusters++;
      clusterSize(kappa(i)) = 0;
    }
    clusterSize(kappa(i))++;
  }
  
  // calculate sufficient statistics
  arma::mat clusterSumy(nclusters, d, arma::fill::zeros);
  arma::cube clusterSumyy(d, d, nclusters, arma::fill::zeros);
  
  for(arma::uword i=0; i<n; i++){
    clusterSumy.row(kappa(i)) = clusterSumy.row(kappa(i)) + y.row(i);
    clusterSumyy.slice(kappa(i)) = clusterSumyy.slice(kappa(i)) + y.row(i).t() * y.row(i);
  }
  
  // update Zeta, Omega, cholOmega, icholOmega and othersOmega
  for(arma::uword i=0; i<nclusters; i++){
    double lambda_new = lambda + 1.0*clusterSize(i);
    arma::colvec m_new = (lambda*m + clusterSumy.row(i).t()) / lambda_new;
    int nu_new = nu + 1.0*clusterSize(i);
    arma::mat tmp_mat = 1.0*clusterSize(i)*m*m.t() - clusterSumy.row(i).t()*clusterSumy.row(i)/lambda - m*clusterSumy.row(i) - clusterSumy.row(i).t()*m.t();
    arma::mat Psi_new = Psi + clusterSumyy.slice(i) + (lambda/lambda_new)*tmp_mat;
    
    arma::mat tmp_omega(d, d);
    arma::mat tmp_cholomega(d, d);
    arma::colvec tmp_zeta(d);
    
    rNormalInverseWishartArma(m_new, lambda_new, nu_new, Psi_new, tmp_omega, tmp_cholomega, tmp_zeta);
    
    arma::mat tmp_icholomega = arma::inv(arma::trimatu(tmp_cholomega));
    double tmp_othersomega = arma::sum(log(tmp_icholomega.diag())) - (double)d/2.0 * log2pi;
    
    Omega.slice(i) = tmp_omega;
    cholOmega.slice(i) = tmp_cholomega;
    icholOmega.slice(i) = tmp_icholomega;
    othersOmega(i) = tmp_othersomega;
    Zeta.col(i) = tmp_zeta;
  }
  
  
  // update alpha: Escobar & West, 1995
  if(updateAlpha){
    // sample eta
    double tmp1 = arma::randg<double>(arma::distr_param(alpha+1.0, 1.0));
    double tmp2 = arma::randg<double>(arma::distr_param(n*1.0, 1.0));
    double eta = tmp1 / (tmp1 + tmp2);
    
    // compute weights
    double tmp_w1 = a0 + nclusters - 1.0;
    double tmp_w2 = b0 - log(eta);
    
    // sample new alpha
    double tmp_u = arma::randu<double>();
    if(tmp_u < (tmp_w1 / (tmp_w1 + n*tmp_w2))) {
      tmp_w1 = tmp_w1 + 1.0;
    }
    alpha = arma::randg<double>(arma::distr_param(tmp_w1, 1.0/tmp_w2));
  }
  
  
  if(useHyperpriors){
    arma::mat sumInvOmega(d, d, arma::fill::zeros);
    arma::colvec sumInvOmegaZeta(d, arma::fill::zeros);
    double sumZetaInvOmegaZeta = 0.0;
    
    for(arma::uword i=0; i<nclusters; i++){
      arma::mat tmp_mat = icholOmega.slice(i) * icholOmega.slice(i).t();
      arma::colvec tmp_vec = tmp_mat * Zeta.col(i);
      double tmp = arma::as_scalar(tmp_vec.t() * Zeta.col(i));
      
      sumInvOmega = sumInvOmega + tmp_mat;
      sumInvOmegaZeta = sumInvOmegaZeta + tmp_vec;
      sumZetaInvOmegaZeta = sumZetaInvOmegaZeta + tmp;
    }
    
    // update m
    arma::mat m_cov = arma::inv_sympd(lambda * sumInvOmega + invS0);
    arma::colvec m_mean = m_cov * (lambda * sumInvOmegaZeta + invS0m0);
    m = arma::mvnrnd(m_mean, m_cov);
    
    // update lambda
    double gamma1_new = 0.5*nclusters*d + gamma1;
    double gamma2_new = sumZetaInvOmegaZeta - 2.0 * arma::as_scalar(m.t() * sumInvOmegaZeta) + arma::as_scalar(m.t() * sumInvOmega * m);
    gamma2_new = gamma2 + 0.5*gamma2_new;
    lambda = arma::randg<double>(arma::distr_param(gamma1_new, 1.0/gamma2_new));
    
    // update Psi
    double nu0_new = 1.0*nclusters*nu + nu0;
    arma::mat Psi0_new = arma::inv_sympd(sumInvOmega + invPsi0);
    Psi = arma::wishrnd(Psi0_new, nu0_new);
  }
}


//------------------------------------------------------------------
// calculate joint density
static void predict_joint_Neal(arma::uword ngrid, arma::uword n, arma::uword d, arma::uword nclusters, arma::mat & ypred, arma::mat & Zeta, arma::cube & icholOmega, arma::colvec & othersOmega, double alpha, arma::colvec & m, double lambda, int nu, arma::mat & Psi, arma::urowvec & clusterSize, arma::mat & evalPDF) {
  
  // generate a sample from NIW prior
  arma::mat prior_omega(d, d);
  arma::mat prior_cholomega(d, d);
  arma::colvec prior_zeta(d);
  
  rNormalInverseWishartArma(m, lambda, nu, Psi, prior_omega, prior_cholomega, prior_zeta);
  
  // calculate likelihoods
  arma::mat yloglik(ngrid*ngrid, nclusters+1);
  for(arma::uword i=0; i<nclusters; i++) {
    arma::colvec tmp_zeta = Zeta.col(i);
    arma::mat rooti = icholOmega.slice(i);
    double other_terms = othersOmega(i);
    
    yloglik.col(i) = dMvnormArma(ypred, ngrid*ngrid, d, tmp_zeta, rooti, other_terms, true);
  }
  arma::mat prior_rooti(d, d);
  double prior_others;
  yloglik.col(nclusters) = dMvnormArma(ypred, ngrid*ngrid, d, prior_zeta, prior_cholomega, prior_rooti, prior_others, true);
  
  // calculate predictive densities
  for(arma::uword i=0; i<(ngrid*ngrid); i++){
    arma::rowvec v = yloglik.row(i) + log(clusterSize(arma::span(0, nclusters)));
    v(nclusters) = yloglik(i, nclusters) + log(alpha);
    
    evalPDF(i) = log_sum_exp(v, false) / (alpha + 1.0*n);
  }
  
}

#endif



