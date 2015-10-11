#include <R_ext/Rdynload.h>
#include "fitter.h"
#include "random.h"

static const R_CMethodDef CEntries[]  = {
    {"EM_fit",                  (DL_FUNC) &EM_fit,                  9},
    {"equicorrelation_fit",     (DL_FUNC) &equicorrelation_fit,     11},
    {"rand_norm",               (DL_FUNC) &rand_norm,               4},
    {"rand_student",            (DL_FUNC) &rand_student,            5},
    {"restricted_center",       (DL_FUNC) &restricted_center,       9},
    {"variance_homogeneity",    (DL_FUNC) &variance_homogeneity,    11},
    {NULL, NULL, 0}
};

void R_init_heavy(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
