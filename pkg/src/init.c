#include "spdep2.h"

static R_CallMethodDef CallEntries[] = {
    {"R_ml_sse_env", (DL_FUNC) &R_ml_sse_env, 2},
    {"R_ml_Jac_env", (DL_FUNC) &R_ml_Jac_env, 2},
    {"do_LLa", (DL_FUNC) &do_LLa, 3},
    {NULL, NULL, 0}
};

/* static R_CMethodDef CEntries[] = {
    {"f_esar_ll", (DL_FUNC) &f_esar_ll, 2}, 
    {"ml_sse_env", (DL_FUNC) &ml_sse_env, 2},
    {"ml_Jac_env", (DL_FUNC) &ml_Jac_env, 2},
    {NULL, NULL, 0}
};*/

#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif

void R_init_spdep2(DllInfo *dll)
{
    R_useDynamicSymbols(dll, FALSE);
/*    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL); */
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
}

