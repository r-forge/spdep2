#ifndef SPDEP2_H
#define SPDEP2_H

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Applic.h>

SEXP R_ml_sse_env(SEXP env, SEXP lambda);
SEXP R_ml_Jac_env(SEXP env, SEXP lambda);

#endif /* SPDEP2_H */



