#include "spdep2.h"

#include <R_ext/Rdynload.h>
#include "Matrix.h"

static R_CallMethodDef CallEntries[] = {
    {"R_ml_sse_env", (DL_FUNC) &R_ml_sse_env, 2},
    {"R_ml_Jac_env", (DL_FUNC) &R_ml_Jac_env, 2},
    {NULL, NULL, 0}
};


#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif

void R_init_spdep2(DllInfo *dll)
{
    R_useDynamicSymbols(dll, FALSE);
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
}

