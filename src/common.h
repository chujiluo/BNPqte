#ifndef GUARD_common_h
#define GUARD_common_h

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <iostream>
#include <vector>
#include <cmath>
#include <typeinfo>
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// log(2*pi)
#define log2pi 1.837877066409345483560659472811

// log(pi)
#define logpi 1.1447298858494001741434273513530587

#endif