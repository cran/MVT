/* ID: shared.c, last updated 2024-09-23, F.Osorio */

#include "base.h"
#include "interface.h"
#include "MVT.h"
#include "shared.h"

void
E_step(double *x, int n, int p, double *center, double *Scatter, FAMILY family, double *distances, double *weights)
{ /* E-step: computation of Mahalanobis distances and weights for the t-distribution */
  double *Root, *z;
  int errcode = 0, job = 0;

  Root = (double *) R_Calloc(p * p, double);
  z    = (double *) R_Calloc(p, double);

  copy_lower(Root, p, Scatter, p, p);
  chol_decomp(Root, p, p, job, &errcode);
  if (errcode)
    error("Cholesky decomposition in E-step gave code %d", errcode);

  for (int i = 0; i < n; i++) {
    copy_vec(z, 1, x + i, n, p);
    distances[i] = mahalanobis(z, p, center, Root);
    weights[i] = do_weight(family, (double) p, distances[i]);
  }

  R_Free(Root); R_Free(z);
}

void
compSymm_pd(double sigma2, double rho, int p, double *mat)
{ /* constructs the equicorrelation matrix (compound symmetry) */
  for (int i = 0; i < p; i++) {
    mat[i * (p + 1)] = sigma2;
    for (int j = i + 1; j < p; j++)
      *(mat + i + j * p) = *(mat + j + i * p) = sigma2 * rho;
  }
}

double
log_Lik(FAMILY family, DIMS dm, double *distances, double *Scatter)
{ /* evaluate the log-likelihood function */
  double *Root, val;
  int errcode = 0, job = 0;

  Root = (double *) R_Calloc(dm->p * dm->p, double);

  copy_lower(Root, dm->p, Scatter, dm->p, dm->p);
  chol_decomp(Root, dm->p, dm->p, job, &errcode);
  if (errcode)
    error("Cholesky decomposition in log-likelihood gave code %d", errcode);

  val  = logLik_kernel(family, dm, distances);
  val -= (double) dm->n * logAbsDet(Root, dm->p, dm->p);
  R_Free(Root);

  return val;
}
