#ifndef GUARD_common_h
#define GUARD_common_h

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstddef>
#include <typeinfo>
#include <map>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#define printf Rprintf

// log(2*pi)
#define log2pi 1.837877066409345483560659472811

// log(pi)
#define logpi 1.1447298858494001741434273513530587

// sqrt(2*pi)
#define RTPI 2.506628274631000502415765284811

#include "rn.h"

#endif