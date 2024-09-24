/* ID: R_init_MVT.c, last updated 2024-09-23, F.Osorio */

#include "base.h"
#include "MVT.h"
#include <R_ext/Rdynload.h>

#define CALLDEF(name, nargs)  {#name, (DL_FUNC) &name, nargs}
#define F77DEF(name, nargs)   {#name, (DL_FUNC) &F77_NAME(name), nargs}

/* estimation for the multivariate t-distribution */
extern void fitter_dispatcher(double *, int *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

/* density and RNG for the multivariate Student-t distribution */
extern void pdf_mstudent(double *, double *, int *, int *, double *, double *, double *);
extern void RNG_mstudent(double *, int *, double *, double *, double *);

/* density for the multivariate normal distribution */
extern void pdf_mnorm(double *, double *, int *, int *, double *, double *);

/* Wilson-Hilferty transformation */
extern void Wilson_Hilferty_chisq(double *, int *, int *, double *);
extern void Wilson_Hilferty_F(double *, int *, int *, double *, double *);

/* registering C symbols */
static const R_CMethodDef CEntries[]  = {
  CALLDEF(fitter_dispatcher,       13),
  CALLDEF(pdf_mnorm,                6),
  CALLDEF(pdf_mstudent,             7),
  CALLDEF(RNG_mstudent,             5),
  CALLDEF(Wilson_Hilferty_chisq,    4),
  CALLDEF(Wilson_Hilferty_F,        5),
  {NULL, NULL, 0}
};

void R_init_MVT(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
