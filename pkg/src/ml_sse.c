#include "spdep2.h"

static int c__1 = 1;


SEXP R_ml_sse_env(SEXP env, SEXP lambda) {

  SEXP res;
  double sse, coef;
  int pc=0;

  coef = NUMERIC_POINTER(lambda)[0];

  sse = ml_sse_env(env, coef);

  PROTECT(res=NEW_NUMERIC(1)); pc++;
  NUMERIC_POINTER(res)[0] = sse;
  UNPROTECT(pc);

  return(res);
}

double ml_sse_env(SEXP env, double coef) {
  SEXP y, x, wy, WX;
  int i, j, k, n, p;
  double *yl, *xlq, *xlqyl, *qy;
  double tol=1e-7, cyl, cxlqyl, sse;
  int *jpvt;
  double *work, *qraux;
  char *trans = "T";
  double one = 1.0, zero = 0.0;

  n = INTEGER_POINTER(findVarInFrame(env, install("n")))[0];
  y = findVarInFrame(env, install("y"));
  x = findVarInFrame(env, install("x"));
  p = length(x) / n;
  wy = findVarInFrame(env, install("wy"));
  WX = findVarInFrame(env, install("WX"));

  yl = (double *) R_alloc(n, sizeof(double));
  xlq = (double *) R_alloc(n*p, sizeof(double));
  qy = (double *) R_alloc(n*p, sizeof(double));
  xlqyl = (double *) R_alloc(p, sizeof(double));
  jpvt = (int *) R_alloc(p, sizeof(int));
  work = (double *) R_alloc(p*2, sizeof(double));
  qraux = (double *) R_alloc(p, sizeof(double));

  for (i=0; i<n; i++) yl[i] = NUMERIC_POINTER(y)[i] - 
    coef * NUMERIC_POINTER(wy)[i];
  for (i=0; i<n*p; i++) xlq[i] = NUMERIC_POINTER(x)[i] - 
    coef * NUMERIC_POINTER(WX)[i];

  F77_CALL(dqrdc2)(xlq, &n, &n, &p, &tol, &k, qraux, jpvt, work);
  if (p != k) warning("Q looses full rank");

  for (i=0; i<n*k; i++) qy[i] = 0.0;
  for (i=0; i<k; i++) qy[(i +(n*i))] = 1.0;

  F77_CALL(dqrqy)(xlq, &n, &k, qraux, qy, &k, qy);

  F77_CALL(dgemv)(trans, &n, &k, &one, qy, &n, yl, &c__1, &zero, xlqyl, &c__1);

  cyl = F77_CALL(ddot)(&n, yl, &c__1, yl, &c__1);

  cxlqyl = F77_CALL(ddot)(&k, xlqyl, &c__1, xlqyl, &c__1);

  sse = cyl - cxlqyl;

  return(sse);

}


#include <R_ext/Rdynload.h>
#include "Matrix.h"

SEXP R_ml_Jac_env(SEXP env, SEXP lambda) {

  SEXP res;
  double jac, coef;
  int pc=0;

  coef = NUMERIC_POINTER(lambda)[0];

  jac = ml_Jac_env(env, coef);

  PROTECT(res=NEW_NUMERIC(1)); pc++;
  NUMERIC_POINTER(res)[0] = jac;
  UNPROTECT(pc);

  return(res);
}

double ml_Jac_env(SEXP env, double coef) {
  SEXP dWd, ndWd, C1p, C1n;
  int n;
  double a, b, mult, nmult, jac;

  n = INTEGER_POINTER(findVarInFrame(env, install("n")))[0];
  a = NUMERIC_POINTER(findVarInFrame(env, install("a")))[0];
  b = NUMERIC_POINTER(findVarInFrame(env, install("b")))[0];
  
  dWd = findVarInFrame(env, install("dWd"));
  ndWd = findVarInFrame(env, install("ndWd"));
  C1p = findVarInFrame(env, install("C1p"));
  C1n = findVarInFrame(env, install("C1n"));
  mult = 1/coef;
  nmult = 1/(-coef);
  if (coef > b) jac = n * log(coef) +
    M_chm_factor_ldetL2(M_chm_factor_update(AS_CHM_FR(C1p),
      AS_CHM_SP(ndWd), mult));
  else if (coef < a) jac = n * log(-(coef)) +
    M_chm_factor_ldetL2(M_chm_factor_update(AS_CHM_FR(C1n),
      AS_CHM_SP(dWd), nmult));
  else jac = 0;

  return(jac);

}

double f_esar_ll(double alpha, SEXP env) {
  double loglik, s2, SSE, Jacobian;
  int n;
  int verbose;

  n = INTEGER_POINTER(findVarInFrame(env, install("n")))[0];
  verbose = LOGICAL_POINTER(findVarInFrame(env, install("verbose")))[0];

  SSE = ml_sse_env(env, alpha);
  s2 = SSE/n;

  Jacobian = ml_Jac_env(env, alpha);

  loglik = (Jacobian - ((n/2) * log(2 * PI)) -
    (n/2) * log(s2) - (1/(2 * (s2))) * SSE);

  if (verbose) Rprintf("coef: %9.6f, SSE: %9.3f, Jacobian: %9.3f, LL: %10.3f\n",
    alpha, SSE, Jacobian, loglik);
  
  return -loglik;
}

SEXP do_LL(SEXP env, SEXP interval, SEXP tol) {
  double low, up, toler;
  double max, obj;
  SEXP res, resnames;
  int pc=0;

  low = NUMERIC_POINTER(interval)[0];
  up = NUMERIC_POINTER(interval)[1];
  toler = NUMERIC_POINTER(tol)[0];
 
  max = Brent_fmin(low, up, (double (*)(double, void*)) f_esar_ll,
    env, toler);
  obj = f_esar_ll(max, env);

  PROTECT(res = NEW_LIST(2)); pc++;
  PROTECT(resnames = NEW_CHARACTER(2)); pc++;
  SET_VECTOR_ELT(res, 0, NEW_NUMERIC(1));
  SET_VECTOR_ELT(res, 1, NEW_NUMERIC(1));
  SET_STRING_ELT(resnames, 0, COPY_TO_USER_STRING("maximum"));
  SET_STRING_ELT(resnames, 1, COPY_TO_USER_STRING("objective"));
  setAttrib(res, R_NamesSymbol, resnames);
  NUMERIC_POINTER(VECTOR_ELT(res, 0))[0] = max;
  NUMERIC_POINTER(VECTOR_ELT(res, 1))[0] = obj;
  UNPROTECT(pc);

  return(res);
}


/*

library(spdep2)
data(columbus)
listw <- nb2listw(col.gal.nb)
res <- ml_sse_setup(CRIME ~ HOVAL + INC, columbus, listw, verbose=TRUE)

optimize(do_LL, interval=c(-1, 1), env=res, interp=TRUE, maximum=TRUE, tol=.Machine$double.eps^0.5)

optimize(do_LL, interval=c(-1, 1), env=res, interp=FALSE, maximum=TRUE, tol=.Machine$double.eps^0.5)

res <- ml_sse_setup(CRIME ~ HOVAL + INC, columbus, listw, verbose=FALSE)

system.time(for (i in 1:1000) optimize(do_LL, interval=c(-1, 1), env=res, interp=TRUE, maximum=TRUE, tol=.Machine$double.eps^0.5))

system.time(for (i in 1:1000) optimize(do_LL, interval=c(-1, 1), env=res, interp=FALSE, maximum=TRUE, tol=.Machine$double.eps^0.5))


data(house)
listw <- nb2listw(LO_nb)
form <- formula(log(price) ~ age + I(age^2) + I(age^3) + log(lotsize) +
   rooms + log(TLA) + beds + syear)
res <- ml_sse_setup(form, house, listw, verbose=FALSE)

optimize(do_LL, interval=c(-1, 1), env=res, interp=TRUE, maximum=TRUE, tol=.Machine$double.eps^0.5)

optimize(do_LL, interval=c(-1, 1), env=res, interp=FALSE, maximum=TRUE, tol=.Machine$double.eps^0.5)

system.time(for (i in 1:10) optimize(do_LL, interval=c(-1, 1), env=res, interp=TRUE, maximum=TRUE, tol=.Machine$double.eps^0.5))

system.time(for (i in 1:10) optimize(do_LL, interval=c(-1, 1), env=res, interp=FALSE, maximum=TRUE, tol=.Machine$double.eps^0.5))


library(spdep)
data(columbus)
listw <- nb2listw(col.gal.nb)
#W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
source("../R/ml_sse.R")
res <- ml_sse_setup(CRIME ~ HOVAL + INC, columbus, listw)

dyn.load("ml_sse.so")

optimize(do_LL, interval=c(-1, 1), env=res, interp=FALSE, maximum=TRUE)

.Call("do_LL", res, c(-1, 1), 1e-9)

ml_sse(0.5, res)
dyn.load("ml_sse.so")
.Call("R_ml_sse_env", res, 0.5)

library(spdep)
source("../R/ml_sse.R")
dyn.load("ml_sse.so")
load("../data/house.RData")

listw <- nb2listw(LO_nb)
W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
form <- formula(log(price) ~ age + I(age^2) + I(age^3) + log(lotsize) +
   rooms + log(TLA) + beds + syear)
res <- ml_sse_setup(form, house, listw)
ml_sse(0.5, res)
.Call("R_ml_sse_env", res, 0.5)

system.time(for (i in 1:100) ml_sse(0.5, res))
system.time(for (i in 1:100) .Call("R_ml_sse_env", res, 0.5))

rho <- seq(-0.9, 0.99, 0.01)
system.time(interp <- sapply(rho, function(x) ml_sse(x, res)))
system.time(compil <- sapply(rho, function(x) .Call("R_ml_sse_env", res, x)))
all.equal(interp, compil)

*/
