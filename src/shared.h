/* ID: shared.h, last updated 2021-02-06, F.Osorio */

#ifndef MVT_SHARED_H
#define MVT_SHARED_H

#include "MVT.h"

/* shared functions for fitters */
extern void compSymm_pd(double, double, int, double *);
extern void E_step(double *, int, int, double *, double *, FAMILY, double *, double *);
extern double log_Lik(FAMILY, DIMS, double *, double *);

#endif /* MVT_SHARED_H */
