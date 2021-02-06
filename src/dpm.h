#ifndef GUARD_dpm_h
#define GUARD_dpm_h

#include "common.h"
#include "random.h"


//------------------------------------------------------------------
// initialize parameters
static void setparam(arma::uword n, arma::uword nclusters, arma::colvec & m, double lambda, int nu, arma::mat & Psi, double alpha, arma::cube & Omega, arma::cube & cholOmega, arma::mat & Zeta, arma::rowvec & lw, arma::rowvec & a_gd, arma::rowvec & b_gd, arma::uvec & kappa){
  // initialize Zeta, Omega and cholOmega
  arma::mat R = arma::chol(arma::inv(Psi));
  rNormalInverseWishartArma(nclusters, m, lambda, nu, Psi, R, Omega, cholOmega, Zeta);
  
  // initialize lw
  rGeneralizedDirichletArma(nclusters, a_gd, b_gd, lw);
  
  // initialize kappa
  rCat(nclusters, n, lw, kappa);
}


//------------------------------------------------------------------
// update (hyper)parameters
static void drawparam(arma::uword n, arma::uword d, arma::uword nclusters, const arma::mat & y, bool updateAlpha, bool useHyperpriors, double a0, double b0, const arma::colvec & m0, const arma::mat & S0, double gamma1, double gamma2, int nu0, const arma::mat & Psi0, double & alpha, arma::colvec & m, double & lambda, int nu, arma::mat & Psi, arma::cube & Omega, arma::cube & cholOmega, arma::cube & icholOmega, arma::colvec & othersOmega, arma::mat & Zeta, arma::rowvec & lw, arma::rowvec & a_gd, arma::rowvec & b_gd, arma::uvec & kappa){
  
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
    // update Zeta, Omega and cholOmega
    if(clusterSize[i]==0){
      arma::mat tmp_omega(d, d);
      arma::mat tmp_cholomega(d, d);
      arma::colvec tmp_zeta(d);
      
      rNormalInverseWishartArma(m, lambda, nu, Psi, tmp_omega, tmp_cholomega, tmp_zeta);
      
      Omega.slice(i) = tmp_omega;
      cholOmega.slice(i) = tmp_cholomega;
      Zeta.col(i) = tmp_zeta;
    } else{
      double lambda_new = lambda + clusterSize[i];
      arma::colvec m_new = (lambda*m + clusterSumy.row(i).t()) / lambda_new;
      int nu_new = nu + clusterSize[i];
      arma::mat Psi_new = Psi + clusterSumyy.slice(i) + 
        ((lambda*clusterSize[i]/lambda_new - clusterSize[i])/pow(clusterSize[i],2))*clusterSumy.row(i).t()*clusterSumy.row(i) +
        (lambda*clusterSize[i]/lambda_new)*(m*m.t() - clusterSumy.row(i).t()*m.t()/clusterSize[i] - m*clusterSumy.row(i)/clusterSize[i]);
      
      arma::mat tmp_omega(d, d);
      arma::mat tmp_cholomega(d, d);
      arma::colvec tmp_zeta(d);
      
      rNormalInverseWishartArma(m_new, lambda_new, nu_new, Psi_new, tmp_omega, tmp_cholomega, tmp_zeta); 
      
      Omega.slice(i) = tmp_omega;
      cholOmega.slice(i) = tmp_cholomega;
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
  
  for(arma::uword i=0; i<nclusters; i++){
    arma::mat tmp_cholomega = cholOmega.slice(i);
    arma::colvec tmp_zeta = Zeta.col(i);
    arma::mat rooti;
    double other_terms;
    
    loglik.col(i) = dMvnormArma(y, n, d, tmp_zeta, tmp_cholomega, rooti, other_terms, true);
    
    icholOmega.slice(i) = rooti;
    othersOmega(i) = other_terms;
  }
  
  for(arma::uword i=0; i<n; i++){
    arma::rowvec lw_tmp = lw + loglik.row(i);
    rCat(nclusters, lw_tmp, kappa(i));
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
      arma::mat tmp_mat = icholOmega.slice(i) * icholOmega.slice(i).t();
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


//------------------------------------------------------------------
// calculate joint density
static void predict_joint(arma::uword ngrid, arma::uword d, arma::uword nclusters, arma::mat & ypred, arma::mat & Zeta, arma::cube icholOmega, arma::colvec othersOmega, arma::rowvec & lw, arma::mat & evalPDF) {
  
  arma::mat yloglik(ngrid*ngrid, nclusters);
  
  for(arma::uword i=0; i<nclusters; i++) {
    arma::colvec tmp_zeta = Zeta.col(i);
    arma::mat rooti = icholOmega.slice(i);
    double other_terms = othersOmega(i);
    
    yloglik.col(i) = dMvnormArma(ypred, ngrid*ngrid, d, tmp_zeta, rooti, other_terms, true);
  }
  
  for(arma::uword i=0; i<(ngrid*ngrid); i++){
    arma::rowvec v = lw + yloglik.row(i);
    evalPDF(i) = log_sum_exp(v, false);  // treating arma::mat as flat and in column major order
  }
  
}
  
  
//------------------------------------------------------------------
// calculate conditional density or cdf
static void predict_conditional(arma::uword ngrid, arma::uword npred, arma::uword d, arma::uword nclusters, arma::colvec & grid, const arma::mat & xpred, arma::mat & Zeta, arma::cube & Omega, arma::rowvec & lw, bool pdf, bool cdf, arma::mat & evalcPDF, arma::mat & evalcCDF) {
  
  if(d == 2) {
    // x is univariate
    
    arma::mat logwx(npred, nclusters);
    arma::rowvec beta0(nclusters);
    arma::rowvec beta(nclusters);
    arma::rowvec sigma2(nclusters);
    
    for(arma::uword i=0; i<nclusters; i++) {
      arma::mat omega = Omega.slice(i);
      arma::colvec zeta = Zeta.col(i);
      
      // beta
      beta(i) = omega.at(0, 1) / omega.at(1, 1);
      
      // beta0
      beta0(i) = zeta(0) - beta(i) * zeta(1);
      
      // sigma2
      sigma2(i) = omega.at(0, 0) - beta(i) * omega.at(1, 0);
      
      // log(pdf(x|zeta2, omega22))
      logwx.col(i) = arma::log_normpdf(xpred, zeta(1), omega.at(1, 1));
    }
    
    // log(w) + log(pdf(x|zeta2, omega22))
    logwx = logwx + arma::repmat(lw, npred, 1);
    
    // sum(exp(logwx.row(i)))
    arma::colvec wx_norm(npred);
    for(arma::uword i=0; i<npred; i++) {
      arma::rowvec tmp = logwx.row(i);
      wx_norm(i) = log_sum_exp(tmp, false);
    }
    
    // mean and var for the conditional kernels
    arma::mat cmean = xpred * beta + arma::repmat(beta0, npred, 1);  // npred x nclusters
    arma::mat cvar = arma::repmat(sigma2, npred, 1);  // npred x nclusters
    
    for(arma::uword j=0; j<ngrid; j++) {
      if(pdf) {
        arma::mat gloglik = arma::log_normpdf(grid(j), cmean, cvar);  // output is npred x nclusters, for grid(j) only
        
        for(arma::uword i=0; i<npred; i++) {
          arma::rowvec tmp_vec = logwx.row(i) + gloglik.row(i);
          evalcPDF(i, j) = log_sum_exp(tmp_vec, false) / wx_norm(i);
        }
      }
      
      if(cdf) {
        arma::mat glogcdf = arma::normcdf(grid(j), cmean, cvar);
        glogcdf = log(glogcdf);
        
        for(arma::uword i=0; i<npred; i++) {
          arma::rowvec tmp_vec = logwx.row(i) + glogcdf.row(i);
          evalcCDF(i, j) = log_sum_exp(tmp_vec, false) / wx_norm(i);
        }
      }
    }
    
  } else {
    // x is multivariate
    
    arma::mat logwx(npred, nclusters);
    arma::rowvec beta0(nclusters);
    arma::mat beta((d-1), nclusters);
    arma::rowvec sigma2(nclusters);
    
    for(arma::uword i=0; i<nclusters; i++) {
      double zeta1 = Zeta.at(0, i);
      arma::colvec zeta2 = Zeta.submat(1, i, (d-1), i);
      arma::mat omega = Omega.slice(i);
      arma::rowvec omega12 = omega.submat(0, 1, 0, (d-1));
      arma::mat omega22 = omega.submat(1, 1, (d-1), (d-1));
      
      arma::mat rooti22((d-1), (d-1));
      double others22;
      
      // log(pdf(x|zeta2, omega22))
      logwx.col(i) = dMvnormArma(xpred, npred, (d-1), zeta2, omega22, rooti22, others22, true);
      
      // beta
      arma::colvec tmp_vec = rooti22 * rooti22.t() * omega12.t();
      beta.col(i) = tmp_vec;
      
      // beta0
      beta0(i) = arma::as_scalar(zeta1 - tmp_vec.t() * zeta2);
      
      // sigma2
      sigma2(i) = omega.at(0, 0) - arma::as_scalar(omega12 * tmp_vec);
    }
    
    // log(w) + log(pdf(x|zeta2, omega22))
    logwx = logwx + arma::repmat(lw, npred, 1);
    
    // sum(exp(logwx.row(i)))
    arma::colvec wx_norm(npred);
    for(arma::uword i=0; i<npred; i++) {
      arma::rowvec tmp = logwx.row(i);
      wx_norm(i) = log_sum_exp(tmp, false);
    }
    
    // mean and var for the conditional kernels
    arma::mat cmean = xpred * beta + arma::repmat(beta0, npred, 1);  // npred x nclusters
    arma::mat cvar = arma::repmat(sigma2, npred, 1);  // npred x nclusters
    
    for(arma::uword j=0; j<ngrid; j++) {
      if(pdf) {
        arma::mat gloglik = arma::log_normpdf(grid(j), cmean, cvar);  // output is npred x nclusters, for grid(j) only
        for(arma::uword i=0; i<npred; i++) {
          arma::rowvec tmp_vec = logwx.row(i) + gloglik.row(i);
          evalcPDF(i, j) = log_sum_exp(tmp_vec, false) / wx_norm(i);
        }
      }
      
      if(cdf) {
        arma::mat glogcdf = arma::normcdf(grid(j), cmean, cvar);
        glogcdf = log(glogcdf);
        for(arma::uword i=0; i<npred; i++) {
          arma::rowvec tmp_vec = logwx.row(i) + glogcdf.row(i);
          evalcCDF(i, j) = log_sum_exp(tmp_vec, false) / wx_norm(i);
        }
      }
    }
  }
  
}


#endif