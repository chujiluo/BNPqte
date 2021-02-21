// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cDPMcdensity
Rcpp::List cDPMcdensity(const arma::uword n, const arma::uword d, const arma::mat& data, const arma::colvec& y, const arma::mat& x, const bool status, const bool pdf, const bool cdf, const bool meanReg, const arma::uword ngrid, const arma::uword npred, const bool hpd, const bool bci, const bool updateAlpha, const bool useHyperpriors, const double a0, const double b0, const arma::colvec& m0, const arma::mat& S0, const double gamma1, const double gamma2, const int nu0, const arma::mat& Psi0, const int nu, const arma::uword nclusters, const arma::uword nskip, const arma::uword ndpost, const arma::uword keepevery, const arma::uword printevery, double alpha, double lambda, Rcpp::Nullable<Rcpp::NumericVector> m_, Rcpp::Nullable<Rcpp::NumericMatrix> Psi_, Rcpp::Nullable<Rcpp::NumericMatrix> Zeta_, Rcpp::Nullable<Rcpp::List> Omega_, Rcpp::Nullable<Rcpp::NumericVector> a_gd_, Rcpp::Nullable<Rcpp::NumericVector> b_gd_, Rcpp::Nullable<Rcpp::NumericVector> lw_, Rcpp::Nullable<Rcpp::IntegerVector> kappa_, Rcpp::Nullable<Rcpp::NumericVector> grid_, Rcpp::Nullable<Rcpp::NumericMatrix> xpred_);
RcppExport SEXP _BNPqte_cDPMcdensity(SEXP nSEXP, SEXP dSEXP, SEXP dataSEXP, SEXP ySEXP, SEXP xSEXP, SEXP statusSEXP, SEXP pdfSEXP, SEXP cdfSEXP, SEXP meanRegSEXP, SEXP ngridSEXP, SEXP npredSEXP, SEXP hpdSEXP, SEXP bciSEXP, SEXP updateAlphaSEXP, SEXP useHyperpriorsSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP m0SEXP, SEXP S0SEXP, SEXP gamma1SEXP, SEXP gamma2SEXP, SEXP nu0SEXP, SEXP Psi0SEXP, SEXP nuSEXP, SEXP nclustersSEXP, SEXP nskipSEXP, SEXP ndpostSEXP, SEXP keepeverySEXP, SEXP printeverySEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP m_SEXP, SEXP Psi_SEXP, SEXP Zeta_SEXP, SEXP Omega_SEXP, SEXP a_gd_SEXP, SEXP b_gd_SEXP, SEXP lw_SEXP, SEXP kappa_SEXP, SEXP grid_SEXP, SEXP xpred_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const bool >::type status(statusSEXP);
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
    rcpp_result_gen = Rcpp::wrap(cDPMcdensity(n, d, data, y, x, status, pdf, cdf, meanReg, ngrid, npred, hpd, bci, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, nu, nclusters, nskip, ndpost, keepevery, printevery, alpha, lambda, m_, Psi_, Zeta_, Omega_, a_gd_, b_gd_, lw_, kappa_, grid_, xpred_));
    return rcpp_result_gen;
END_RCPP
}
// cDPMdensity
Rcpp::List cDPMdensity(const arma::uword n, const arma::uword d, const arma::mat& y, const bool status, const bool prediction, const arma::uword ngrid, const bool updateAlpha, const bool useHyperpriors, const double a0, const double b0, const arma::colvec& m0, const arma::mat& S0, const double gamma1, const double gamma2, const int nu0, const arma::mat& Psi0, const int nu, const arma::uword nclusters, const arma::uword nskip, const arma::uword ndpost, const arma::uword keepevery, const arma::uword printevery, double alpha, double lambda, Rcpp::Nullable<Rcpp::NumericVector> m_, Rcpp::Nullable<Rcpp::NumericMatrix> Psi_, Rcpp::Nullable<Rcpp::NumericMatrix> Zeta_, Rcpp::Nullable<Rcpp::List> Omega_, Rcpp::Nullable<Rcpp::NumericVector> a_gd_, Rcpp::Nullable<Rcpp::NumericVector> b_gd_, Rcpp::Nullable<Rcpp::NumericVector> lw_, Rcpp::Nullable<Rcpp::IntegerVector> kappa_, Rcpp::Nullable<Rcpp::NumericVector> grid1_, Rcpp::Nullable<Rcpp::NumericVector> grid2_);
RcppExport SEXP _BNPqte_cDPMdensity(SEXP nSEXP, SEXP dSEXP, SEXP ySEXP, SEXP statusSEXP, SEXP predictionSEXP, SEXP ngridSEXP, SEXP updateAlphaSEXP, SEXP useHyperpriorsSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP m0SEXP, SEXP S0SEXP, SEXP gamma1SEXP, SEXP gamma2SEXP, SEXP nu0SEXP, SEXP Psi0SEXP, SEXP nuSEXP, SEXP nclustersSEXP, SEXP nskipSEXP, SEXP ndpostSEXP, SEXP keepeverySEXP, SEXP printeverySEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP m_SEXP, SEXP Psi_SEXP, SEXP Zeta_SEXP, SEXP Omega_SEXP, SEXP a_gd_SEXP, SEXP b_gd_SEXP, SEXP lw_SEXP, SEXP kappa_SEXP, SEXP grid1_SEXP, SEXP grid2_SEXP) {
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
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type grid1_(grid1_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type grid2_(grid2_SEXP);
    rcpp_result_gen = Rcpp::wrap(cDPMdensity(n, d, y, status, prediction, ngrid, updateAlpha, useHyperpriors, a0, b0, m0, S0, gamma1, gamma2, nu0, Psi0, nu, nclusters, nskip, ndpost, keepevery, printevery, alpha, lambda, m_, Psi_, Zeta_, Omega_, a_gd_, b_gd_, lw_, kappa_, grid1_, grid2_));
    return rcpp_result_gen;
END_RCPP
}
// cpDPMcdensity
Rcpp::List cpDPMcdensity(const arma::uword ngrid, const arma::uword npred, arma::colvec& grid, arma::mat& xpred, const arma::uword d, const arma::uword nclusters, const arma::uword ndpost, const Rcpp::List& ZetaList, const Rcpp::List& OmegaList, const arma::mat& lwList, const bool pdf, const bool cdf, const bool meanReg, const bool hpd, const bool bci);
RcppExport SEXP _BNPqte_cpDPMcdensity(SEXP ngridSEXP, SEXP npredSEXP, SEXP gridSEXP, SEXP xpredSEXP, SEXP dSEXP, SEXP nclustersSEXP, SEXP ndpostSEXP, SEXP ZetaListSEXP, SEXP OmegaListSEXP, SEXP lwListSEXP, SEXP pdfSEXP, SEXP cdfSEXP, SEXP meanRegSEXP, SEXP hpdSEXP, SEXP bciSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(cpDPMcdensity(ngrid, npred, grid, xpred, d, nclusters, ndpost, ZetaList, OmegaList, lwList, pdf, cdf, meanReg, hpd, bci));
    return rcpp_result_gen;
END_RCPP
}
// cpDPMdensity
Rcpp::List cpDPMdensity(const arma::uword ngrid, const arma::colvec& grid1, const arma::colvec& grid2, const arma::uword d, const arma::uword nclusters, const arma::uword ndpost, const Rcpp::List& ZetaList, const Rcpp::List& OmegaList, const arma::mat& lwList);
RcppExport SEXP _BNPqte_cpDPMdensity(SEXP ngridSEXP, SEXP grid1SEXP, SEXP grid2SEXP, SEXP dSEXP, SEXP nclustersSEXP, SEXP ndpostSEXP, SEXP ZetaListSEXP, SEXP OmegaListSEXP, SEXP lwListSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(cpDPMdensity(ngrid, grid1, grid2, d, nclusters, ndpost, ZetaList, OmegaList, lwList));
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
double timesTwo(arma::mat& x, arma::uword i);
RcppExport SEXP _BNPqte_timesTwo(SEXP xSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x, i));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BNPqte_cDPMcdensity", (DL_FUNC) &_BNPqte_cDPMcdensity, 41},
    {"_BNPqte_cDPMdensity", (DL_FUNC) &_BNPqte_cDPMdensity, 34},
    {"_BNPqte_cpDPMcdensity", (DL_FUNC) &_BNPqte_cpDPMcdensity, 15},
    {"_BNPqte_cpDPMdensity", (DL_FUNC) &_BNPqte_cpDPMdensity, 9},
    {"_BNPqte_timesTwo", (DL_FUNC) &_BNPqte_timesTwo, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_BNPqte(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
