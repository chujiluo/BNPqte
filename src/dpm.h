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
  b_gd(0) = alpha + 1.0*n;
  
  
  //Rcpp::Rcout << "Updating Zeta and Omega.." << std::endl;
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
      double lambda_new = lambda + 1.0*clusterSize[i];
      arma::colvec m_new = (lambda*m + clusterSumy.row(i).t()) / lambda_new;
      int nu_new = nu + 1.0*clusterSize[i];
      arma::mat tmp_mat = 1.0*clusterSize[i]*m*m.t() - clusterSumy.row(i).t()*clusterSumy.row(i)/lambda - m*clusterSumy.row(i) - clusterSumy.row(i).t()*m.t();
      arma::mat Psi_new = Psi + clusterSumyy.slice(i) + (lambda/lambda_new)*tmp_mat;
      
      //Rcpp::Rcout << "lambda: " << lambda << ", lambda_new: " << lambda_new << std::endl;
      //Rcpp::Rcout << "clusterSumyy: " << clusterSumyy(0, 0, i) << " " << clusterSumyy(0, 1, i) << " " << clusterSumyy(1, 1, i) << std::endl;
      //Rcpp::Rcout << "tmp_mat: " << tmp_mat(0, 0) << " " << tmp_mat(0, 1) << " " << tmp_mat(1, 1) << std::endl;
      //Rcpp::Rcout << "Psi: " << Psi(0, 0) << " " << Psi(0, 1) << " " << Psi(1, 1) << std::endl;
      
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
      a_gd(i) = 1.0 + 1.0*clusterSize[i];
      
      if(i==0){
        b_gd(i) = b_gd(i) - 1.0*clusterSize[i];
      } else{
        b_gd(i) = b_gd(i-1) - 1.0*clusterSize[i];
      }
    }
  }
  
  //Rcpp::Rcout << "Updating lw.." << std::endl;
  // update w (weight)
  rGeneralizedDirichletArma(nclusters, a_gd, b_gd, lw);
  
  //Rcpp::Rcout << "Updating kappa.." << std::endl;
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
    rCat(nclusters, lw_tmp, i, kappa);
  }
  
  //Rcpp::Rcout << "Updating alpha.." << std::endl;
  // update alpha
  if(updateAlpha){
    double a0_new = a0 + 1.0*nclusters - 1.0;
    double b0_new = b0 - lw(nclusters-1);
    alpha = arma::randg<double>(arma::distr_param(a0_new, 1.0/b0_new));
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
    
    //Rcpp::Rcout << "Updating m.." << std::endl;
    // update m
    arma::mat invS0 = arma::inv_sympd(S0);
    arma::mat m_cov = arma::inv_sympd(lambda * sumInvOmega + invS0);
    arma::colvec m_mean = m_cov * (lambda * sumInvOmegaZeta + invS0 * m0);
    m = arma::mvnrnd(m_mean, m_cov);
    
    //Rcpp::Rcout << "Updating lambda.." << std::endl;
    // update lambda
    double gamma1_new = 0.5*nclusters*d + gamma1;
    //double gamma1_new = 0.5*nclusters + gamma1;
    double gamma2_new = sumZetaInvOmegaZeta - 2.0 * arma::as_scalar(m.t() * sumInvOmegaZeta) + arma::as_scalar(m.t() * sumInvOmega * m);
    gamma2_new = gamma2 + 0.5*gamma2_new;
    
    //Rcpp::Rcout << "gamma1: " << gamma1_new << ", gamma2: " << gamma2_new << std::endl;
    lambda = arma::randg<double>(arma::distr_param(gamma1_new, 1.0/gamma2_new));
    
    //Rcpp::Rcout << "Updating Psi.." << std::endl;
    // update Psi
    double nu0_new = 1.0*nclusters*nu + nu0;
    arma::mat Psi0_new = arma::inv_sympd(sumInvOmega + arma::inv_sympd(Psi0));
    Psi = arma::wishrnd(Psi0_new, nu0_new);
  }
}


//------------------------------------------------------------------
// calculate joint density
static void predict_joint(arma::uword ngrid, arma::uword d, arma::uword nclusters, arma::mat & ypred, arma::mat & Zeta, arma::cube & icholOmega, arma::colvec & othersOmega, arma::rowvec & lw, arma::mat & evalPDF) {
  
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
static void predict_conditional(arma::uword ngrid, arma::uword npred, arma::uword d, arma::uword nclusters, arma::colvec & grid, const arma::mat & xpred, arma::mat & Zeta, arma::cube & Omega, arma::rowvec & lw, bool pdf, bool cdf, arma::mat & evalcPDF, arma::mat & evalcCDF, arma::colvec & evalcMean) {
  
  if(d == 2) {
    // x is univariate
    
    arma::mat xloglik(npred, nclusters);
    arma::rowvec beta0(nclusters);
    arma::rowvec beta(nclusters);
    arma::rowvec sigma(nclusters);
    
    for(arma::uword i=0; i<nclusters; i++) {
      arma::mat omega = Omega.slice(i);
      arma::colvec zeta = Zeta.col(i);
      
      beta(i) = omega(0, 1) / omega(1, 1);
      beta0(i) = zeta(0) - beta(i) * zeta(1);
      
      double tmp1 = omega(0, 0) - beta(i) * omega(1, 0);
      sigma(i) = sqrt(tmp1);
      
      double tmp2 = omega(1, 1);
      tmp2 = sqrt(tmp2);
      xloglik.col(i) = arma::log_normpdf(xpred, zeta(1), tmp2);
    }
    
    // log_sum_exp(lw+xloglik)
    
    arma::colvec lwx_norm(npred);
    for(arma::uword i=0; i<npred; i++) {
      arma::rowvec tmp_vec = lw + xloglik.row(i);
      lwx_norm(i) = log_sum_exp(tmp_vec, true);
    }
    
    // mean and variance for the conditional kernels
    arma::mat cmean = xpred * beta + arma::repmat(beta0, npred, 1);  // npred x nclusters
    arma::mat csd = arma::repmat(sigma, npred, 1);  // npred x nclusters
    
    //Rcpp::Rcout << "max and min beta: " << beta.max() << " " << beta.min() << std::endl;
    //Rcpp::Rcout << "max and min beta0: " << beta0.max() << " " << beta0.min() << std::endl;
    //Rcpp::Rcout << "max and min sigma2: " << sigma2.max() << " " << sigma2.min() << std::endl;
    
    for(arma::uword i=0; i<npred; i++) {
      // evalcMean
      arma::rowvec tmp1 = lw + xloglik.row(i);
      tmp1 = exp(tmp1);
      double tmp2 = arma::as_scalar(tmp1 * cmean.row(i).t());
      evalcMean(i) = tmp2/exp(lwx_norm(i));
      
      /*
      arma::rowvec tmp1 = lw + xloglik.row(i) + log(cmean.row(i));
      double tmp2 = log_sum_exp(tmp1, true);
      
      arma::rowvec tmp3 = log(cmean.row(i));
      Rcpp::Rcout << "cmean(1): " << cmean(i,1) << " log(cmean(1)): " << tmp3(1) << std::endl;
      evalcMean(i) = exp(tmp2 - lwx_norm(i));
      
      Rcpp::Rcout << "tmp2: " << tmp2 << " evalcMean: " << evalcMean(i) << std::endl;
      */
      
      /* 
      //Rcpp::Rcout << "here: " << exp(xloglik(i,0)-lwx_norm(i)) << " " << exp(xloglik(i,1)-lwx_norm(i)) << std::endl;
      arma::rowvec tmp1 = lw + log(cmean.row(i));
      double tmp2 = log_sum_exp(tmp1, true);
      evalcMean(i) = exp(tmp2);
      */
      for(arma::uword j=0; j<ngrid; j++) {
        // evalcPDF
        if(pdf) {
          arma::mat gloglik = arma::log_normpdf(grid(j), cmean, csd);  // output is npred x nclusters, for grid(j) only
          
          arma::rowvec tmp_vec = lw + xloglik.row(i) + gloglik.row(i);
          double tmp = log_sum_exp(tmp_vec, true);
          evalcPDF(i, j) = exp(tmp - lwx_norm(i));
          
          /*
          arma::rowvec tmp_vec = lw + gloglik.row(i);
          double tmp = log_sum_exp(tmp_vec, true);
          evalcPDF(i, j) = exp(tmp);
           */
        }
        
        // evalcCDF
        if(cdf) {
          arma::mat glogcdf = arma::normcdf(grid(j), cmean, csd);
          glogcdf = log(glogcdf);
          
          arma::rowvec tmp_vec = lw + xloglik.row(i) + glogcdf.row(i);
          double tmp = log_sum_exp(tmp_vec, true);
          evalcCDF(i, j) = exp(tmp - lwx_norm(i));
           
          /*
          arma::rowvec tmp_vec = lw + glogcdf.row(i);
          double tmp = log_sum_exp(tmp_vec, true);
          evalcCDF(i, j) = exp(tmp);
           */
        }
      }
    }
    
  } else {
    // x is multivariate
    
    arma::mat xloglik(npred, nclusters);
    arma::rowvec beta0(nclusters);
    arma::mat beta((d-1), nclusters);
    arma::rowvec sigma(nclusters);
    
    for(arma::uword i=0; i<nclusters; i++) {
      double zeta1 = Zeta(0, i);
      arma::colvec zeta2 = Zeta.submat(1, i, (d-1), i);
      
      arma::mat omega = Omega.slice(i);
      arma::rowvec omega12 = omega.submat(0, 1, 0, (d-1));
      arma::mat omega22 = omega.submat(1, 1, (d-1), (d-1));
      arma::mat cholomega22 = arma::chol(omega22);
      arma::mat rooti22((d-1), (d-1));
      double others22;
      
      // xloglik
      xloglik.col(i) = dMvnormArma(xpred, npred, (d-1), zeta2, cholomega22, rooti22, others22, true);
      
      // beta
      arma::colvec tmp_vec = rooti22 * rooti22.t() * omega12.t();
      beta.col(i) = tmp_vec;
      
      // beta0
      beta0(i) = arma::as_scalar(zeta1 - tmp_vec.t() * zeta2);
      
      // sigma2
      double tmp = omega(0, 0) - arma::as_scalar(omega12 * tmp_vec);
      sigma(i) = sqrt(tmp);
    }
    
    // log_sum_exp(lw+xloglik)
    arma::colvec lwx_norm(npred);
    for(arma::uword i=0; i<npred; i++) {
      arma::rowvec tmp_vec = lw + xloglik.row(i);
      lwx_norm(i) = log_sum_exp(tmp_vec, true);
    }
    
    // mean and variance for the conditional kernels
    arma::mat cmean = xpred * beta + arma::repmat(beta0, npred, 1);  // npred x nclusters
    arma::mat csd = arma::repmat(sigma, npred, 1);  // npred x nclusters
    
    for(arma::uword i=0; i<npred; i++) {
      // evalcMean
      arma::rowvec tmp1 = lw + xloglik.row(i) + log(cmean.row(i));
      double tmp2 = log_sum_exp(tmp1, true);
      evalcMean(i) = exp(tmp2 - lwx_norm(i));
      
      for(arma::uword j=0; j<ngrid; j++) {
        // evalcPDF
        if(pdf) {
          arma::mat gloglik = arma::log_normpdf(grid(j), cmean, csd);  // output is npred x nclusters, for grid(j) only
          
          arma::rowvec tmp_vec = lw + xloglik.row(i) + gloglik.row(i);
          double tmp = log_sum_exp(tmp_vec, true);
          evalcPDF(i, j) = exp(tmp - lwx_norm(i));
        }
        
        // evalcCDF
        if(cdf) {
          arma::mat glogcdf = arma::normcdf(grid(j), cmean, csd);
          glogcdf = log(glogcdf);
          
          arma::rowvec tmp_vec = lw + xloglik.row(i) + glogcdf.row(i);
          double tmp = log_sum_exp(tmp_vec, true);
          evalcCDF(i, j) = exp(tmp - lwx_norm(i));
        }
      }
      
    }
    
  }
  
}


#endif