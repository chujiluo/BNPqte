// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cDPMcdensity
Rcpp::List cDPMcdensity(const arma::uword n, const arma::uword d, const arma::mat& data, const arma::colvec& y, const arma::mat& x, const bool status, const bool diag, const bool pdf, const bool cdf, const bool meanReg, const arma::uword ngrid, const arma::uword npred, const bool hpd, const bool bci, const bool updateAlpha, const bool useHyperpriors, const double a0, const double b0, const arma::colvec& m0, const arma::mat& S0, const double gamma1, const double gamma2, const int nu0, const arma::mat& Psi0, const int nu, const arma::uword nclusters, const arma::uword nskip, const arma::uword ndpost, const arma::uword keepevery, const arma::uword printevery, double alpha, double lambda, Rcpp::Nullable<Rcpp::NumericVector> m_, Rcpp::Nullable<Rcpp::NumericMatrix> Psi_, Rcpp::Nullable<Rcpp::NumericMatrix> Zeta_, Rcpp::Nullable<Rcpp::List> Omega_, Rcpp::Nullable<Rcpp::NumericVector> a_gd_, Rcpp::Nullable<Rcpp::NumericVector> b_gd_, Rcpp::Nullable<Rcpp::NumericVector> lw_, Rcpp::Nullable<Rcpp::IntegerVector> kappa_, Rcpp::Nullable<Rcpp::NumericVector> grid_, Rcpp::Nullable<Rcpp::NumericMatrix> xpred_);
RcppExport SEXP _BNPqte_cDPMcdensity(SEXP nSEXP, SEXP dSEXP, SEXP dataSEXP, SEXP ySEXP, SEXP xSEXP, SEXP statusSEXP, SEXP diagSEXP, SEXP pdfSEXP, SEXP cdfSEXP, SEXP meanRegSEXP, SEXP ngridSEXP, SEXP npredSEXP, SEXP hpdSEXP, SEXP bciSEXP, SEXP updateAlphaSEXP, SEXP useHyperpriorsSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP m0SEXP, SEXP S0SEXP, SEXP gamma1SEXP, SEXP gamma2SEXP, SEXP nu0SEXP, SEXP Psi0SEXP, SEXP nuSEXP, SEXP nclustersSEXP, SEXP nskipSEXP, SEXP ndpostSEXP, SEXP keepeverySEXP, SEXP printeverySEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP m_SEXP, SEXP Psi_SEXP, SEXP Zeta_SEXP, SEXP Omega_SEXP, SEXP a_gd_SEXP, SEXP b_gd_SEXP, SEXP lw_SEXP, SEXP kappa_SEXP, SEXP grid_SEXP, SEXP xpred_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const bool >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const bool >::type diag(diagSEXP);
    Rcpp::traits::input_parameter< const bool >::type pdf(pdfSEXP);
    Rcpp::traits::input_parameter< const bool >::type cdf(cdfSEXP);
    Rcpp::traits::input_parameter< const bool >::type meanReg(meanRegSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type ngrid(ngridSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type npred(npredSEXP);
    Rcpp::traits::input_parameter< const bool >::type hpd(hpdSEXP);
    Rcpp::traits::input_parameter< const bool >::type bci(bciSEXP);
    Rcpp::traits::input_parameter< const bool >::type updateAlpha(updateAlphaSEXP);
    Rcpp::traits::input_parameter< const bool >::type useHyperpriors(useHyperpriorsSEXP);
    Rcpp::traits::input_parameter< const double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< const double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< const double >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< const double >::type gamma2(gamma2SEXP);
    Rcpp::traits::input_parameter< const int >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Psi0(Psi0SEXP);
    Rcpp::traits::input_parameter< const int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nclusters(nclustersSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nskip(nskipSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type ndpost(ndpostSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type keepevery(keepeverySEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type printevery(printeverySEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type m_(m_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Psi_(Psi_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Zeta_(Zeta_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type Omega_(Omega_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type a_gd_(a_gd_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type b_gd_(b_gd_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type lw_(lw_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type kappa_(kappa_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type grid_(grid_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type xpred_(xpred_SEXP);
    rcpp_result_gen = Rcpp::wrap(cDPMcdensity(n, d, data, y, x, status, diag, pdf, cdf, meanReg, ngrid, npred, hpd, bci, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, nu, nclusters, nskip, ndpost, keepevery, printevery, alpha, lambda, m_, Psi_, Zeta_, Omega_, a_gd_, b_gd_, lw_, kappa_, grid_, xpred_));
    return rcpp_result_gen;
END_RCPP
}
// cDPMdensity
Rcpp::List cDPMdensity(const arma::uword n, const arma::uword d, const arma::mat& y, const bool status, const bool diag, const bool prediction, const arma::uword ngrid, const bool updateAlpha, const bool useHyperpriors, const double a0, const double b0, const arma::colvec& m0, const arma::mat& S0, const double gamma1, const double gamma2, const int nu0, const arma::mat& Psi0, const int nu, const arma::uword nclusters, const arma::uword nskip, const arma::uword ndpost, const arma::uword keepevery, const arma::uword printevery, double& alpha, double& lambda, arma::colvec& m, arma::mat& Psi, arma::rowvec& a_gd, arma::rowvec& b_gd, Rcpp::Nullable<Rcpp::NumericMatrix> Zeta_, Rcpp::Nullable<Rcpp::List> Omega_, Rcpp::Nullable<Rcpp::NumericVector> lw_, Rcpp::Nullable<Rcpp::IntegerVector> kappa_, Rcpp::Nullable<Rcpp::NumericVector> grid1_, Rcpp::Nullable<Rcpp::NumericVector> grid2_);
RcppExport SEXP _BNPqte_cDPMdensity(SEXP nSEXP, SEXP dSEXP, SEXP ySEXP, SEXP statusSEXP, SEXP diagSEXP, SEXP predictionSEXP, SEXP ngridSEXP, SEXP updateAlphaSEXP, SEXP useHyperpriorsSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP m0SEXP, SEXP S0SEXP, SEXP gamma1SEXP, SEXP gamma2SEXP, SEXP nu0SEXP, SEXP Psi0SEXP, SEXP nuSEXP, SEXP nclustersSEXP, SEXP nskipSEXP, SEXP ndpostSEXP, SEXP keepeverySEXP, SEXP printeverySEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP mSEXP, SEXP PsiSEXP, SEXP a_gdSEXP, SEXP b_gdSEXP, SEXP Zeta_SEXP, SEXP Omega_SEXP, SEXP lw_SEXP, SEXP kappa_SEXP, SEXP grid1_SEXP, SEXP grid2_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const bool >::type diag(diagSEXP);
    Rcpp::traits::input_parameter< const bool >::type prediction(predictionSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type ngrid(ngridSEXP);
    Rcpp::traits::input_parameter< const bool >::type updateAlpha(updateAlphaSEXP);
    Rcpp::traits::input_parameter< const bool >::type useHyperpriors(useHyperpriorsSEXP);
    Rcpp::traits::input_parameter< const double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< const double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< const double >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< const double >::type gamma2(gamma2SEXP);
    Rcpp::traits::input_parameter< const int >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Psi0(Psi0SEXP);
    Rcpp::traits::input_parameter< const int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nclusters(nclustersSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nskip(nskipSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type ndpost(ndpostSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type keepevery(keepeverySEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type printevery(printeverySEXP);
    Rcpp::traits::input_parameter< double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< arma::rowvec& >::type a_gd(a_gdSEXP);
    Rcpp::traits::input_parameter< arma::rowvec& >::type b_gd(b_gdSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Zeta_(Zeta_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type Omega_(Omega_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type lw_(lw_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type kappa_(kappa_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type grid1_(grid1_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type grid2_(grid2_SEXP);
    rcpp_result_gen = Rcpp::wrap(cDPMdensity(n, d, y, status, diag, prediction, ngrid, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, nu, nclusters, nskip, ndpost, keepevery, printevery, alpha, lambda, m, Psi, a_gd, b_gd, Zeta_, Omega_, lw_, kappa_, grid1_, grid2_));
    return rcpp_result_gen;
END_RCPP
}
// cDPMdensityNeal
Rcpp::List cDPMdensityNeal(const arma::uword n, const arma::uword d, const arma::mat& y, const bool status, const bool prediction, const arma::uword ngrid, const bool updateAlpha, const bool useHyperpriors, const double a0, const double b0, const arma::colvec& m0, const arma::mat& S0, const double gamma1, const double gamma2, const int nu0, const arma::mat& Psi0, const int nu, arma::uword& nclusters, const arma::uword nskip, const arma::uword ndpost, const arma::uword keepevery, const arma::uword printevery, double& alpha, double& lambda, arma::colvec& m, arma::mat& Psi, Rcpp::Nullable<Rcpp::NumericMatrix> Zeta_, Rcpp::Nullable<Rcpp::List> Omega_, Rcpp::Nullable<Rcpp::IntegerVector> kappa_, Rcpp::Nullable<Rcpp::NumericVector> grid1_, Rcpp::Nullable<Rcpp::NumericVector> grid2_);
RcppExport SEXP _BNPqte_cDPMdensityNeal(SEXP nSEXP, SEXP dSEXP, SEXP ySEXP, SEXP statusSEXP, SEXP predictionSEXP, SEXP ngridSEXP, SEXP updateAlphaSEXP, SEXP useHyperpriorsSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP m0SEXP, SEXP S0SEXP, SEXP gamma1SEXP, SEXP gamma2SEXP, SEXP nu0SEXP, SEXP Psi0SEXP, SEXP nuSEXP, SEXP nclustersSEXP, SEXP nskipSEXP, SEXP ndpostSEXP, SEXP keepeverySEXP, SEXP printeverySEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP mSEXP, SEXP PsiSEXP, SEXP Zeta_SEXP, SEXP Omega_SEXP, SEXP kappa_SEXP, SEXP grid1_SEXP, SEXP grid2_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const bool >::type prediction(predictionSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type ngrid(ngridSEXP);
    Rcpp::traits::input_parameter< const bool >::type updateAlpha(updateAlphaSEXP);
    Rcpp::traits::input_parameter< const bool >::type useHyperpriors(useHyperpriorsSEXP);
    Rcpp::traits::input_parameter< const double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< const double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< const double >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< const double >::type gamma2(gamma2SEXP);
    Rcpp::traits::input_parameter< const int >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Psi0(Psi0SEXP);
    Rcpp::traits::input_parameter< const int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::uword& >::type nclusters(nclustersSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nskip(nskipSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type ndpost(ndpostSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type keepevery(keepeverySEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type printevery(printeverySEXP);
    Rcpp::traits::input_parameter< double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Zeta_(Zeta_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type Omega_(Omega_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type kappa_(kappa_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type grid1_(grid1_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type grid2_(grid2_SEXP);
    rcpp_result_gen = Rcpp::wrap(cDPMdensityNeal(n, d, y, status, prediction, ngrid, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, nu, nclusters, nskip, ndpost, keepevery, printevery, alpha, lambda, m, Psi, Zeta_, Omega_, kappa_, grid1_, grid2_));
    return rcpp_result_gen;
END_RCPP
}
// cDPMmdensity
Rcpp::List cDPMmdensity(const arma::uword n, const arma::uword d, const arma::mat& data, const arma::colvec& y, const arma::mat& x, const bool diag, const bool pdf, const bool cdf, const arma::uword ngrid, const arma::colvec& grid, const arma::uword npred, const arma::mat& xpred, const bool updateAlpha, const bool useHyperpriors, double alpha, const double a0, const double b0, double lambda, const double gamma1, const double gamma2, const int nu0, const int nu, const arma::uword nclusters, const arma::uword nskip, const arma::uword ndpost, const arma::uword keepevery, const arma::rowvec& diri, const arma::colvec& probs, const arma::uword nprobs);
RcppExport SEXP _BNPqte_cDPMmdensity(SEXP nSEXP, SEXP dSEXP, SEXP dataSEXP, SEXP ySEXP, SEXP xSEXP, SEXP diagSEXP, SEXP pdfSEXP, SEXP cdfSEXP, SEXP ngridSEXP, SEXP gridSEXP, SEXP npredSEXP, SEXP xpredSEXP, SEXP updateAlphaSEXP, SEXP useHyperpriorsSEXP, SEXP alphaSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP lambdaSEXP, SEXP gamma1SEXP, SEXP gamma2SEXP, SEXP nu0SEXP, SEXP nuSEXP, SEXP nclustersSEXP, SEXP nskipSEXP, SEXP ndpostSEXP, SEXP keepeverySEXP, SEXP diriSEXP, SEXP probsSEXP, SEXP nprobsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const bool >::type diag(diagSEXP);
    Rcpp::traits::input_parameter< const bool >::type pdf(pdfSEXP);
    Rcpp::traits::input_parameter< const bool >::type cdf(cdfSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type ngrid(ngridSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type npred(npredSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xpred(xpredSEXP);
    Rcpp::traits::input_parameter< const bool >::type updateAlpha(updateAlphaSEXP);
    Rcpp::traits::input_parameter< const bool >::type useHyperpriors(useHyperpriorsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< const double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< const double >::type gamma2(gamma2SEXP);
    Rcpp::traits::input_parameter< const int >::type nu0(nu0SEXP);
    Rcpp::traits::input_parameter< const int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nclusters(nclustersSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nskip(nskipSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type ndpost(ndpostSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type keepevery(keepeverySEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type diri(diriSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nprobs(nprobsSEXP);
    rcpp_result_gen = Rcpp::wrap(cDPMmdensity(n, d, data, y, x, diag, pdf, cdf, ngrid, grid, npred, xpred, updateAlpha, useHyperpriors, alpha, a0, b0, lambda, gamma1, gamma2, nu0, nu, nclusters, nskip, ndpost, keepevery, diri, probs, nprobs));
    return rcpp_result_gen;
END_RCPP
}
// cpDPMcdensity
Rcpp::List cpDPMcdensity(const arma::uword ngrid, const arma::uword npred, arma::colvec& grid, arma::mat& xpred, const arma::uword d, const arma::uword nclusters, const arma::uword ndpost, const Rcpp::List& ZetaList, const Rcpp::List& OmegaList, const arma::mat& lwList, const bool pdf, const bool cdf, const bool meanReg, const bool hpd, const bool bci, const arma::uword printevery);
RcppExport SEXP _BNPqte_cpDPMcdensity(SEXP ngridSEXP, SEXP npredSEXP, SEXP gridSEXP, SEXP xpredSEXP, SEXP dSEXP, SEXP nclustersSEXP, SEXP ndpostSEXP, SEXP ZetaListSEXP, SEXP OmegaListSEXP, SEXP lwListSEXP, SEXP pdfSEXP, SEXP cdfSEXP, SEXP meanRegSEXP, SEXP hpdSEXP, SEXP bciSEXP, SEXP printeverySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword >::type ngrid(ngridSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type npred(npredSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xpred(xpredSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nclusters(nclustersSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type ndpost(ndpostSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type ZetaList(ZetaListSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type OmegaList(OmegaListSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lwList(lwListSEXP);
    Rcpp::traits::input_parameter< const bool >::type pdf(pdfSEXP);
    Rcpp::traits::input_parameter< const bool >::type cdf(cdfSEXP);
    Rcpp::traits::input_parameter< const bool >::type meanReg(meanRegSEXP);
    Rcpp::traits::input_parameter< const bool >::type hpd(hpdSEXP);
    Rcpp::traits::input_parameter< const bool >::type bci(bciSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type printevery(printeverySEXP);
    rcpp_result_gen = Rcpp::wrap(cpDPMcdensity(ngrid, npred, grid, xpred, d, nclusters, ndpost, ZetaList, OmegaList, lwList, pdf, cdf, meanReg, hpd, bci, printevery));
    return rcpp_result_gen;
END_RCPP
}
// cpDPMdensity
Rcpp::List cpDPMdensity(const arma::uword ngrid, const arma::colvec& grid1, const arma::colvec& grid2, const arma::uword d, const arma::uword nclusters, const arma::uword ndpost, const Rcpp::List& ZetaList, const Rcpp::List& OmegaList, const arma::mat& lwList, const arma::uword printevery);
RcppExport SEXP _BNPqte_cpDPMdensity(SEXP ngridSEXP, SEXP grid1SEXP, SEXP grid2SEXP, SEXP dSEXP, SEXP nclustersSEXP, SEXP ndpostSEXP, SEXP ZetaListSEXP, SEXP OmegaListSEXP, SEXP lwListSEXP, SEXP printeverySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword >::type ngrid(ngridSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type grid1(grid1SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type grid2(grid2SEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nclusters(nclustersSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type ndpost(ndpostSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type ZetaList(ZetaListSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type OmegaList(OmegaListSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lwList(lwListSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type printevery(printeverySEXP);
    rcpp_result_gen = Rcpp::wrap(cpDPMdensity(ngrid, grid1, grid2, d, nclusters, ndpost, ZetaList, OmegaList, lwList, printevery));
    return rcpp_result_gen;
END_RCPP
}
// cpDPMdensityNeal
Rcpp::List cpDPMdensityNeal(const arma::uword ngrid, const arma::colvec& grid1, const arma::colvec& grid2, const arma::uword d, const arma::uword n, const bool updateAlpha, const bool useHyperpriors, const Rcpp::List& ZetaList, const Rcpp::List& OmegaList, const arma::mat& kappaList, const arma::colvec& nclusterList, const arma::colvec& alphaList, const arma::colvec& lambdaList, const arma::mat& mList, const Rcpp::List& PsiList, const int nu, const arma::uword ndpost, const arma::uword printevery);
RcppExport SEXP _BNPqte_cpDPMdensityNeal(SEXP ngridSEXP, SEXP grid1SEXP, SEXP grid2SEXP, SEXP dSEXP, SEXP nSEXP, SEXP updateAlphaSEXP, SEXP useHyperpriorsSEXP, SEXP ZetaListSEXP, SEXP OmegaListSEXP, SEXP kappaListSEXP, SEXP nclusterListSEXP, SEXP alphaListSEXP, SEXP lambdaListSEXP, SEXP mListSEXP, SEXP PsiListSEXP, SEXP nuSEXP, SEXP ndpostSEXP, SEXP printeverySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword >::type ngrid(ngridSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type grid1(grid1SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type grid2(grid2SEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< const bool >::type updateAlpha(updateAlphaSEXP);
    Rcpp::traits::input_parameter< const bool >::type useHyperpriors(useHyperpriorsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type ZetaList(ZetaListSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type OmegaList(OmegaListSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type kappaList(kappaListSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type nclusterList(nclusterListSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type alphaList(alphaListSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type lambdaList(lambdaListSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mList(mListSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type PsiList(PsiListSEXP);
    Rcpp::traits::input_parameter< const int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type ndpost(ndpostSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type printevery(printeverySEXP);
    rcpp_result_gen = Rcpp::wrap(cpDPMdensityNeal(ngrid, grid1, grid2, d, n, updateAlpha, useHyperpriors, ZetaList, OmegaList, kappaList, nclusterList, alphaList, lambdaList, mList, PsiList, nu, ndpost, printevery));
    return rcpp_result_gen;
END_RCPP
}
// rnorm_cpp
Rcpp::List rnorm_cpp(double s, int N);
RcppExport SEXP _BNPqte_rnorm_cpp(SEXP sSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(rnorm_cpp(s, N));
    return rcpp_result_gen;
END_RCPP
}
// simdisc
void simdisc(arma::rowvec& prob, arma::uword n, arma::uword m, arma::uword& val);
RcppExport SEXP _BNPqte_simdisc(SEXP probSEXP, SEXP nSEXP, SEXP mSEXP, SEXP valSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::uword& >::type val(valSEXP);
    simdisc(prob, n, m, val);
    return R_NilValue;
END_RCPP
}
// test
void test(arma::rowvec& prob, arma::uword n, arma::uword m);
RcppExport SEXP _BNPqte_test(SEXP probSEXP, SEXP nSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type m(mSEXP);
    test(prob, n, m);
    return R_NilValue;
END_RCPP
}
// test_randi
void test_randi(arma::uword m);
RcppExport SEXP _BNPqte_test_randi(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type m(mSEXP);
    test_randi(m);
    return R_NilValue;
END_RCPP
}
// rdirichlet
arma::mat rdirichlet(arma::uword n, arma::colvec& probs);
RcppExport SEXP _BNPqte_rdirichlet(SEXP nSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichlet(n, probs));
    return rcpp_result_gen;
END_RCPP
}
// credible_interval
arma::mat credible_interval(arma::uword n, arma::colvec& x, arma::colvec& alphas, std::string type);
RcppExport SEXP _BNPqte_credible_interval(SEXP nSEXP, SEXP xSEXP, SEXP alphasSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(credible_interval(n, x, alphas, type));
    return rcpp_result_gen;
END_RCPP
}
// quantile_fun
arma::mat quantile_fun(arma::uword ngrid, arma::uword nprobs, arma::uword ndpost, const arma::colvec& grid, const arma::colvec& probs, arma::mat& cdfs);
RcppExport SEXP _BNPqte_quantile_fun(SEXP ngridSEXP, SEXP nprobsSEXP, SEXP ndpostSEXP, SEXP gridSEXP, SEXP probsSEXP, SEXP cdfsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type ngrid(ngridSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type nprobs(nprobsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type ndpost(ndpostSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type cdfs(cdfsSEXP);
    rcpp_result_gen = Rcpp::wrap(quantile_fun(ngrid, nprobs, ndpost, grid, probs, cdfs));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP clbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP cpbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP cwbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP mc_cores_openmp();

static const R_CallMethodDef CallEntries[] = {
    {"_BNPqte_cDPMcdensity", (DL_FUNC) &_BNPqte_cDPMcdensity, 42},
    {"_BNPqte_cDPMdensity", (DL_FUNC) &_BNPqte_cDPMdensity, 35},
    {"_BNPqte_cDPMdensityNeal", (DL_FUNC) &_BNPqte_cDPMdensityNeal, 31},
    {"_BNPqte_cDPMmdensity", (DL_FUNC) &_BNPqte_cDPMmdensity, 29},
    {"_BNPqte_cpDPMcdensity", (DL_FUNC) &_BNPqte_cpDPMcdensity, 16},
    {"_BNPqte_cpDPMdensity", (DL_FUNC) &_BNPqte_cpDPMdensity, 10},
    {"_BNPqte_cpDPMdensityNeal", (DL_FUNC) &_BNPqte_cpDPMdensityNeal, 18},
    {"_BNPqte_rnorm_cpp", (DL_FUNC) &_BNPqte_rnorm_cpp, 2},
    {"_BNPqte_simdisc", (DL_FUNC) &_BNPqte_simdisc, 4},
    {"_BNPqte_test", (DL_FUNC) &_BNPqte_test, 3},
    {"_BNPqte_test_randi", (DL_FUNC) &_BNPqte_test_randi, 1},
    {"_BNPqte_rdirichlet", (DL_FUNC) &_BNPqte_rdirichlet, 2},
    {"_BNPqte_credible_interval", (DL_FUNC) &_BNPqte_credible_interval, 4},
    {"_BNPqte_quantile_fun", (DL_FUNC) &_BNPqte_quantile_fun, 6},
    {"clbart",                    (DL_FUNC) &clbart,                    24},
    {"cpbart",                    (DL_FUNC) &cpbart,                    27},
    {"cwbart",                    (DL_FUNC) &cwbart,                    31},
    {"mc_cores_openmp",           (DL_FUNC) &mc_cores_openmp,            0},
    {NULL, NULL, 0}
};

RcppExport void R_init_BNPqte(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
