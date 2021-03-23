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
      
      // // swap cluster kappa(i) and the last cluster nclusters-1
      // std::swap(clusterSize(kappa(i)), clusterSize(nclusters-1));
      // Omega.slice(kappa(i)).swap(Omega.slice(nclusters-1));
      // cholOmega.slice(kappa(i)).swap(cholOmega.slice(nclusters-1));
      // icholOmega.slice(kappa(i)).swap(icholOmega.slice(nclusters-1));
      // std::swap(othersOmega(kappa(i)), othersOmega(nclusters-1));
      // Zeta.swap_cols(kappa(i), nclusters-1);
      // for (arma::uword j=0; j<n; j++){
      //   if (kappa(j) == (nclusters-1)){
      //     kappa(j) = kappa(i);
      //   }
      // }
      // kappa(i) = nclusters - 1;
      // 
      // nclusters--;
      
      
      // relabel the clusters
      arma::uword old_cluster = kappa(i);
      
      arma::colvec ZetaWork = Zeta.col(old_cluster);
      arma::mat OmegaWork = Omega.slice(old_cluster);
      arma::mat cholOmegaWork = cholOmega.slice(old_cluster);
      arma::mat icholOmegaWork = icholOmega.slice(old_cluster);
      double othersOmegaWork = othersOmega(old_cluster);
      
      if (old_cluster < (nclusters - 1)) {
        
        for (arma::uword k=old_cluster+1; k<nclusters; k++){
          for (arma::uword j=0; j < n; j++){
            if (kappa(j) == k){
              kappa(j) = k - 1;
            }  
          }
          Zeta.col(k-1) = Zeta.col(k);
          Omega.slice(k-1) = Omega.slice(k);
          cholOmega.slice(k-1) = cholOmega.slice(k);
          icholOmega.slice(k-1) = icholOmega.slice(k);
          othersOmega(k-1) = othersOmega(k);
          
          clusterSize(k-1) = clusterSize(k);
        }
        
        kappa(i) = nclusters - 1;
        
        Zeta.col(nclusters-1) = ZetaWork;
        Omega.slice(nclusters-1) = OmegaWork;
        cholOmega.slice(nclusters-1) = cholOmegaWork;
        icholOmega.slice(nclusters-1) = icholOmegaWork;
        othersOmega(nclusters-1) = othersOmegaWork;
        
        clusterSize(nclusters-1) = 1;
        
      }
      nclusters--;
      
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
    //arma::rowvec work_weights = exp(log_weights);
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


//------------------------------------------------------------------
// calculate conditional density or cdf
static void predict_conditional_Neal(arma::uword ngrid, arma::uword npred, arma::uword n, arma::uword d, arma::uword nclusters, arma::colvec & grid, arma::mat & xpred, arma::mat & Zeta, arma::cube & Omega, double alpha, arma::colvec & m, double lambda, int nu, arma::mat & Psi, arma::urowvec & clusterSize, bool pdf, bool cdf, bool meanReg, arma::mat & evalcPDF, arma::mat & evalcCDF, arma::colvec & evalcMean) {
  
  // generate stick-breaking weights with epsilon = 0.01 tolerance
  std::vector<double> std_weights;
  double current_weight = 0.0;
  double total_weight = 0.0;
  double remaining_weight = 1.0;
  arma::uword cnt = 0;
  
  while((1.0 - total_weight) > 0.01) {
    double tmp = rBetaArma(1.0, alpha + (double)n);
    current_weight = remaining_weight * tmp;
    total_weight = total_weight + current_weight;
    remaining_weight = remaining_weight * (1.0 - tmp);
    std_weights.push_back(current_weight);
    cnt++;
  }
  arma::rowvec weights(std_weights);
  
  // sample cnt clusters from 0 ~ nclusters (nclusters+1 different clusters)
  arma::rowvec probs = arma::conv_to<arma::rowvec>::from(clusterSize(arma::span(0, nclusters))) / (alpha + 1.0*n);
  probs(nclusters) = alpha / (alpha + 1.0*n);
  
  arma::urowvec clusters(cnt);
  arma::uword lastCnt;
  simdisc(cnt, (nclusters+1), (nclusters+1), probs, lastCnt, clusters);
  
  // working variables
  arma::mat xlikWork(npred, lastCnt);
  arma::mat cmeanWork(npred, lastCnt);
  arma::mat csdWork(npred, lastCnt);
  arma::cube glikWork(npred, lastCnt, ngrid);
  arma::cube gcdfWork(npred, lastCnt, ngrid);
  
  arma::mat xlik(npred, nclusters);
  arma::mat cmean(npred, nclusters);
  arma::mat csd(npred, nclusters);
  arma::cube glik(npred, nclusters, ngrid);
  arma::cube gcdf(npred, nclusters, ngrid);
  
  if(d == 2) {
    // x is univariate
    
    // sample lastCnt Zetas and Omegas from NIW prior, and calculate related densities given these Zetas and Omegas
    if(lastCnt > 0) {
      // generate lastCnt samples from NIW prior
      arma::mat R = arma::chol(arma::inv(Psi));
      arma::cube OmegaWork(d, d, lastCnt);
      arma::cube cholOmegaWork(d, d, lastCnt);
      arma::mat ZetaWork(d, lastCnt);
      
      rNormalInverseWishartArma(lastCnt, m, lambda, nu, Psi, R, OmegaWork, cholOmegaWork, ZetaWork);
      
      // calculate densitites
      arma::rowvec beta0Work(lastCnt);
      arma::rowvec betaWork(lastCnt);
      arma::rowvec sigmaWork(lastCnt);
      
      for(arma::uword k=0; k<lastCnt; k++) {
        arma::mat omegaWork = OmegaWork.slice(k);
        arma::colvec zetaWork = ZetaWork.col(k);
        
        betaWork(k) = omegaWork(0, 1) / omegaWork(1, 1);
        beta0Work(k) = zetaWork(0) - betaWork(k) * zetaWork(1);
        sigmaWork(k) = sqrt((omegaWork(0, 0) - betaWork(k) * omegaWork(1, 0)));
        
        double tmpWork = sqrt(omegaWork(1, 1));
        xlikWork.col(k) = arma::log_normpdf(xpred, zetaWork(1), tmpWork);
      }
      xlikWork = exp(xlikWork);
      
      cmeanWork = xpred * betaWork + arma::repmat(beta0Work, npred, 1);
      csdWork = arma::repmat(sigmaWork, npred, 1);
      
      if(pdf) {
        for(arma::uword j=0; j<ngrid; j++)
          glikWork.slice(j) = arma::log_normpdf(grid(j), cmeanWork, csdWork);  // (i, k): log f(y_j | x_i*b_k, s_k)
        
        glikWork = exp(glikWork);
      }
      
      if(cdf) {
        for(arma::uword j=0; j<ngrid; j++) {
          gcdfWork.slice(j) = arma::normcdf(grid(j), cmeanWork, csdWork);
        }
      }
    }
    
    // calculate log marginal likelihood for xpred and conditional log likelihood / cdf for grid, using non-empty clusters' Zeta and Omega
    arma::rowvec beta0(nclusters);
    arma::rowvec beta(nclusters);
    arma::rowvec sigma(nclusters);
    
    for(arma::uword k=0; k<nclusters; k++) {
      arma::mat omega = Omega.slice(k);
      arma::colvec zeta = Zeta.col(k);
      
      beta(k) = omega(0, 1) / omega(1, 1);
      beta0(k) = zeta(0) - beta(k) * zeta(1);
      sigma(k) = sqrt((omega(0, 0) - beta(k) * omega(1, 0)));
      
      double tmp = sqrt(omega(1, 1));
      xlik.col(k) = arma::log_normpdf(xpred, zeta(1), tmp);
    }
    xlik = exp(xlik);
    
    cmean = xpred * beta + arma::repmat(beta0, npred, 1);  // npred x nclusters
    csd = arma::repmat(sigma, npred, 1);  // npred x nclusters
    
    if(pdf) {
      for(arma::uword j=0; j<ngrid; j++)
        glik.slice(j) = arma::log_normpdf(grid(j), cmean, csd);  // (i, k): log f(y_j | x_i*b_k, s_k)
    }
    glik = exp(glik);
    
    if(cdf) {
      for(arma::uword j=0; j<ngrid; j++) {
        gcdf.slice(j) = arma::normcdf(grid(j), cmean, csd);
      }
    }
    
  } else {
    // x is multivariate
    
    // sample lastCnt Zetas and Omegas from NIW prior, and calculate related densities given these Zetas and Omegas
    if(lastCnt > 0) {
      // generate lastCnt samples from NIW prior
      arma::mat R = arma::chol(arma::inv(Psi));
      arma::cube OmegaWork(d, d, lastCnt);
      arma::cube cholOmegaWork(d, d, lastCnt);
      arma::mat ZetaWork(d, lastCnt);
      
      rNormalInverseWishartArma(lastCnt, m, lambda, nu, Psi, R, OmegaWork, cholOmegaWork, ZetaWork);
      
      // calculate densitites
      arma::rowvec beta0Work(lastCnt);
      arma::mat betaWork((d-1), lastCnt);
      arma::rowvec sigmaWork(lastCnt);
      
      for(arma::uword k=0; k<lastCnt; k++) {
        double zeta1Work = ZetaWork(0, k);
        arma::colvec zeta2Work = ZetaWork.submat(1, k, (d-1), k);
        
        arma::mat omegaWork = OmegaWork.slice(k);
        arma::rowvec omega12Work = omegaWork.submat(0, 1, 0, (d-1));
        arma::mat omega22Work = omegaWork.submat(1, 1, (d-1), (d-1));
        arma::mat cholomega22Work = arma::chol(omega22Work);
        arma::mat rooti22Work((d-1), (d-1));
        double others22Work;
        
        xlikWork.col(k) = dMvnormArma(xpred, npred, (d-1), zeta2Work, cholomega22Work, rooti22Work, others22Work, true);
        
        arma::colvec tmp_vec = rooti22Work * rooti22Work.t() * omega12Work.t();
        betaWork.col(k) = tmp_vec;
        beta0Work(k) = arma::as_scalar(zeta1Work - tmp_vec.t() * zeta2Work);
        sigmaWork(k) = sqrt(omegaWork(0, 0) - arma::as_scalar(omega12Work * tmp_vec));
      }
      xlikWork = exp(xlikWork);
      
      cmeanWork = xpred * betaWork + arma::repmat(beta0Work, npred, 1);
      csdWork = arma::repmat(sigmaWork, npred, 1);
      
      if(pdf) {
        for(arma::uword j=0; j<ngrid; j++)
          glikWork.slice(j) = arma::log_normpdf(grid(j), cmeanWork, csdWork);  // (i, k): log f(y_j | x_i*b_k, s_k)
        
        glikWork = exp(glikWork);
      }
      
      if(cdf) {
        for(arma::uword j=0; j<ngrid; j++) {
          gcdfWork.slice(j) = arma::normcdf(grid(j), cmeanWork, csdWork);
        }
      }
    }
    
    // calculate log marginal likelihood for xpred and conditional log likelihood / cdf for grid, using non-empty clusters' Zeta and Omega
    arma::rowvec beta0(nclusters);
    arma::mat beta((d-1), nclusters);
    arma::rowvec sigma(nclusters);
    
    for(arma::uword k=0; k<nclusters; k++) {
      double zeta1 = Zeta(0, k);
      arma::colvec zeta2 = Zeta.submat(1, k, (d-1), k);
      
      arma::mat omega = Omega.slice(k);
      arma::rowvec omega12 = omega.submat(0, 1, 0, (d-1));
      arma::mat omega22 = omega.submat(1, 1, (d-1), (d-1));
      arma::mat cholomega22 = arma::chol(omega22);
      arma::mat rooti22((d-1), (d-1));
      double others22;
      
      xlik.col(k) = dMvnormArma(xpred, npred, (d-1), zeta2, cholomega22, rooti22, others22, true);
      
      arma::colvec tmp_vec = rooti22 * rooti22.t() * omega12.t();
      beta.col(k) = tmp_vec;
      beta0(k) = arma::as_scalar(zeta1 - tmp_vec.t() * zeta2);
      sigma(k) = sqrt(omega(0, 0) - arma::as_scalar(omega12 * tmp_vec));
    }
    xlik = exp(xlik);
    
    cmean = xpred * beta + arma::repmat(beta0, npred, 1);  // npred x nclusters
    csd = arma::repmat(sigma, npred, 1);  // npred x nclusters
    
    if(pdf) {
      for(arma::uword j=0; j<ngrid; j++)
        glik.slice(j) = arma::log_normpdf(grid(j), cmean, csd);  // (i, k): log f(y_j | x_i*b_k, s_k)
    }
    glik = exp(glik);
    
    if(cdf) {
      for(arma::uword j=0; j<ngrid; j++) {
        gcdf.slice(j) = arma::normcdf(grid(j), cmean, csd);
      }
    }  
  }
  
  // denominator of the conditional expectation / pdf / cdf
  arma::colvec denom(npred, arma::fill::zeros);
  
  // evaluation
  arma::uword WorkCnt = 0;
  for(arma::uword k=0; k<cnt; k++) {
    if(clusters(k) < nclusters) {
      denom = denom + weights(k) * xlik.col(clusters(k));
      
      // evalcMean
      if(meanReg) {
        evalcMean = evalcMean + weights(k) * xlik.col(clusters(k)) % cmean.col(clusters(k));
      }
      
      // evalcPDF
      if(pdf) {
        for(arma::uword j=0; j<ngrid; j++) {
          evalcPDF.col(j) = evalcPDF.col(j) + weights(k) * xlik.col(clusters(k)) % glik.slice(j).col(clusters(k));
        }
      }
      
      // evalcCDF
      if(cdf) {
        for(arma::uword j=0; j<ngrid; j++) {
          evalcCDF.col(j) = evalcCDF.col(j) + weights(k) * xlik.col(clusters(k)) % gcdf.slice(j).col(clusters(k));
        }
      }
    } else {
      denom = denom + weights(k) * xlikWork.col(WorkCnt);
      
      // evalcMean
      if(meanReg) {
        evalcMean = evalcMean + weights(k) * xlikWork.col(WorkCnt) % cmeanWork.col(WorkCnt);
      }
      
      // evalcPDF
      if(pdf) {
        for(arma::uword j=0; j<ngrid; j++) {
          evalcPDF.col(j) = evalcPDF.col(j) + weights(k) * xlikWork.col(WorkCnt) % glikWork.slice(j).col(WorkCnt);
        }
      }
      
      // evalcCDF
      if(cdf) {
        for(arma::uword j=0; j<ngrid; j++) {
          evalcCDF.col(j) = evalcCDF.col(j) + weights(k) * xlikWork.col(WorkCnt) % gcdfWork.slice(j).col(WorkCnt);
        }
      }
      
      WorkCnt++;
    }
  }
  
  if(meanReg)
    evalcMean = evalcMean / denom;
  
  if(pdf) {
    for(arma::uword j=0; j<ngrid; j++) {
      evalcPDF.col(j) = evalcPDF.col(j) / denom;
    }
  }
  
  if(cdf) {
    for(arma::uword j=0; j<ngrid; j++) {
      evalcCDF.col(j) = evalcCDF.col(j) / denom;
    }
  }
}


//------------------------------------------------------------------
// calculate marginal density or cdf of y
static void predict_marginal_Neal(arma::uword ngrid, arma::uword npred, arma::uword n, arma::uword d, arma::uword nclusters, arma::uword nprobs, arma::colvec & grid, arma::mat & xpred, const arma::colvec & probs, arma::mat & Zeta, arma::cube & Omega, double alpha, arma::colvec & m, double lambda, int nu, arma::mat & Psi, arma::urowvec & clusterSize, bool pdf, bool cdf, const arma::rowvec & diri, arma::rowvec & evalpdf, arma::rowvec & evalcdf, arma::rowvec & evalqtls) {
  
  // calculate conditional density or cdf of y
  arma::mat evalcPDF(npred, ngrid, arma::fill::zeros);
  arma::mat evalcCDF(npred, ngrid, arma::fill::zeros);
  arma::colvec evalcMean(npred);
  
  predict_conditional_Neal(ngrid, npred, n, d, nclusters, grid, xpred, Zeta, Omega, alpha, m, lambda, nu, Psi, clusterSize, 
                           pdf, cdf, false, evalcPDF, evalcCDF, evalcMean);
  

  // calculate marginal density from conditional density using Bayesian boostrap
  if(pdf) {
    evalpdf = diri * evalcPDF;
  }
  
  if(cdf) {
    evalcdf = diri * evalcCDF;
    quantile_fun(ngrid, nprobs, probs, evalcdf, grid, evalqtls);
  }
  
}


#endif



