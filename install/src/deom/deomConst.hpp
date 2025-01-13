#ifndef DEOMCONST_H_
#define DEOMCONST_H_

#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_LAPACK
#include "armadillo"

using namespace arma;

static const double    deom_pi = datum::pi; // armadillo's const
static const cx_double deom_ci = cx_double(0.0,1.0); 
static const cx_double deom_c1 = cx_double(1.0,0.0); 

static const double    deom_cm2unit = 1.0;
static const double    deom_kt2unit = 1.0;
static const double    deom_fs2unit = 1.0;

#endif
