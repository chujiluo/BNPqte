#ifndef GUARD_lambda
#define GUARD_lambda

#include "common.h"

using R::pnorm;

RcppExport SEXP cdraw_lambda_i(SEXP lambda, SEXP mean, SEXP kmax, SEXP thin);

//draw lambda from its (infinite mixture) prior
double draw_lambda_prior(double *psii, int kmax, rn& gen);

//Metropolis-Hastings algorithm for drawing lambda from its full conditional -- uses proposals from the prior
double draw_lambda_i(double lambda_old, double xbeta, int kmax, int thin, rn& gen);

#endif
