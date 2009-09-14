#include "spdep2.h"

/** */
static int c__1 = 1;

/** */
typedef struct {
  int n, p, verbose;
  SEXP y, x, wy, WX;
  double *yl, *xlq, *xlqyl, *qy;
  int *jpvt;
  double *work, *qraux;
  double a, b;
  CHM_FR C1p, C1n;
  CHM_SP dWd, ndWd;
} Jac_SSE_info;


/**
 * Calculate the sum of squared errors term for spatial regression
 * using an environment to hold data
 *
 * @param env pointer to an SEXP environment
 * @param coef current value of coefficient being optimzed
 * 
 * @return double, value of SSE for current coef
 *
 */
static double ml_sse_env(SEXP env, double coef) {
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


/**
 * Calculate the sum of squared errors term for spatial regression
 * using a structure to hold data
 *
 * @param INFO pointer to a Jac_SSE_info structure 
 * @param coef current value of coefficient being optimzed
 * 
 * @return double, value of SSE for current coef
 * 
 */
static double ml_sse_info(Jac_SSE_info *INFO, double coef) {
  int i, j, k;
  double tol=1e-7, cyl, cxlqyl, sse;
  char *trans = "T";
  double one = 1.0, zero = 0.0;

  for (i=0; i<INFO->n; i++)
    INFO->yl[i] = NUMERIC_POINTER(INFO->y)[i] - 
      coef * NUMERIC_POINTER(INFO->wy)[i];
  for (i=0; i<INFO->n*INFO->p; i++) 
    INFO->xlq[i] = NUMERIC_POINTER(INFO->x)[i] - 
      coef * NUMERIC_POINTER(INFO->WX)[i];

  F77_CALL(dqrdc2)(INFO->xlq, &INFO->n, &INFO->n, &INFO->p, &tol, &k,
    INFO->qraux, INFO->jpvt, INFO->work);
  if (INFO->p != k) warning("Q looses full rank");

  for (i=0; i<INFO->n*k; i++) INFO->qy[i] = 0.0;
  for (i=0; i<k; i++) INFO->qy[(i +(INFO->n*i))] = 1.0;

  F77_CALL(dqrqy)(INFO->xlq, &INFO->n, &k, INFO->qraux, INFO->qy, &k,
    INFO->qy);

  F77_CALL(dgemv)(trans, &INFO->n, &k, &one, INFO->qy, &INFO->n, INFO->yl,
    &c__1, &zero, INFO->xlqyl, &c__1);

  cyl = F77_CALL(ddot)(&INFO->n, INFO->yl, &c__1, INFO->yl, &c__1);

  cxlqyl = F77_CALL(ddot)(&k, INFO->xlqyl, &c__1, INFO->xlqyl, &c__1);

  sse = cyl - cxlqyl;

  return(sse);
}


/**
 * SEXP wrapper for the SSE term for spatial regression
 * using an environment to hold data
 *
 * @param env pointer to a SEXP environment
 * @param lambda pointer to a SEXP numeric coefficient value
 * 
 * @return pointer to a SEXP numeric SSE value
 *
 */
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


/**
 * Calculate the Jacobian term for spatial regression
 * using an environment to hold data
 *
 * @param env pointer to a SEXP environment
 * @param coef current value of coefficient being optimzed
 * 
 * @return double, value of the Jacobian for current coef
 * 
 */
static double ml_Jac_env(SEXP env, double coef) {
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

/**
 * Calculate the Jacobian term for spatial regression
 * using a structure to hold data
 *
 * @param INFO pointer to a Jac_SSE_info structure 
 * @param coef current value of coefficient being optimzed
 * 
 * @return double, value of the Jacobian for current coef
 * 
 */
static double ml_Jac_info(Jac_SSE_info *INFO, double coef) {
  double mult, nmult, jac;
  
  mult = 1/coef;
  nmult = 1/(-coef);
  if (coef > INFO->b) jac = INFO->n * log(coef) +
    M_chm_factor_ldetL2(M_chm_factor_update(INFO->C1p,
      INFO->ndWd, mult));
  else if (coef < INFO->a) jac = INFO->n * log(-(coef)) +
    M_chm_factor_ldetL2(M_chm_factor_update(INFO->C1n,
      INFO->dWd, nmult));
  else jac = 0;

  return(jac);

}


/**
 * SEXP wrapper for the Jacobian term for spatial regression
 * using an environment to hold data
 *
 * @param env pointer to a SEXP environment
 * @param lambda pointer to a SEXP numeric coefficient value
 * 
 * @return pointer to a SEXP numeric Jacobian value
 * 
 */
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


/**
 * Calculate the log likelihood term for spatial regression
 * using a structure to hold data
 *
 * @param alpha current value of coefficient being optimzed
 * @param INFO pointer to a Jac_SSE_info structure
 * 
 * @return double, value of the log likelihood for current coef
 *
 */
static double f_esar_ll(double alpha, Jac_SSE_info *INFO) {
  double loglik, s2, SSE, Jacobian;

  SSE = ml_sse_info(INFO, alpha);
  s2 = SSE/INFO->n;

  Jacobian = ml_Jac_info(INFO, alpha); 
  if (!R_FINITE(Jacobian)) {
	    warning("NA/Inf replaced by maximum positive value");
	    return(DBL_MAX);
  } 

  loglik = (Jacobian - ((INFO->n/2) * log(2 * PI)) -
    ((INFO->n/2) * log(s2)) - (1/(2 * (s2))) * SSE);

  if (INFO->verbose)
    Rprintf("val: %9.6f, SSE: %9.3f, Jacobian: %9.3f, LL: %10.3f\n",
    alpha, SSE, Jacobian, loglik);
  
  return -loglik;
}


/**
 * SEXP function to run Brent_fmin, first populating a Jac_SSE_info
 * structure from a SEXP environment to pass to the objective
 * function f_esar_ll
 *
 * @param env pointer to a SEXP environment
 * @param interval SEXP numeric vector with min and max interval bounds
 * @param tol SEXP numeric tolerance for stopping fmin
 * 
 * @return SEXP list of maximum and objective values
 *
 */
SEXP do_LLa(SEXP env, SEXP interval, SEXP tol) {
  double low, up, toler;
  double max, obj;
  SEXP res, resnames;
  int pc=0;
  Jac_SSE_info INFO, *infptr;

  INFO.verbose = LOGICAL_POINTER(findVarInFrame(env, install("verbose")))[0];

  INFO.n = INTEGER_POINTER(findVarInFrame(env, install("n")))[0];
  INFO.y = findVarInFrame(env, install("y"));
  INFO.x = findVarInFrame(env, install("x"));
  INFO.p = length(INFO.x) / INFO.n;
  INFO.wy = findVarInFrame(env, install("wy"));
  INFO.WX = findVarInFrame(env, install("WX"));

  INFO.yl = (double *) R_alloc(INFO.n, sizeof(double));
  INFO.xlq = (double *) R_alloc(INFO.n*INFO.p, sizeof(double));
  INFO.qy = (double *) R_alloc(INFO.n*INFO.p, sizeof(double));
  INFO.xlqyl = (double *) R_alloc(INFO.p, sizeof(double));
  INFO.jpvt = (int *) R_alloc(INFO.p, sizeof(int));
  INFO.work = (double *) R_alloc(INFO.p*2, sizeof(double));
  INFO.qraux = (double *) R_alloc(INFO.p, sizeof(double));

  INFO.a = NUMERIC_POINTER(findVarInFrame(env, install("a")))[0];
  INFO.b = NUMERIC_POINTER(findVarInFrame(env, install("b")))[0];
  INFO.dWd = AS_CHM_SP(findVarInFrame(env, install("dWd")));
  INFO.ndWd = AS_CHM_SP(findVarInFrame(env, install("ndWd")));
  INFO.C1p = AS_CHM_FR(findVarInFrame(env, install("C1p")));
  INFO.C1n = AS_CHM_FR(findVarInFrame(env, install("C1n")));

  infptr = &INFO;

  low = NUMERIC_POINTER(interval)[0];
  up = NUMERIC_POINTER(interval)[1];
  toler = NUMERIC_POINTER(tol)[0];
  
 
  max = Brent_fmin(low, up, (double (*)(double, void*)) f_esar_ll,
    infptr, toler);
  obj = f_esar_ll(max, infptr);

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



