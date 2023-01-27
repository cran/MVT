/* ID: wilson-hilferty.c, last updated 2021-02-06, F.Osorio */

#include "base.h"
#include "interface.h"

void
Wilson_Hilferty_chisq(double *distances, int *n, int *p, double *z)
{ /* Wilson-Hilferty transformation for chi-squared variables */
  WH_chisq(distances, *n, *p, z);
}

void
Wilson_Hilferty_F(double *distances, int *n, int *p, double *eta, double *z)
{ /* Wilson-Hilferty transformation for F variables */
  WH_F(distances, *n, *p, *eta, z);
}
