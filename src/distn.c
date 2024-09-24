/* ID: distn.c, last updated 2024-09-23, F.Osorio */

#include "MVT.h"
#include "interface.h"

/* ========================================================================== *
 * pdf for the multivariate Student-t distribution
 * ========================================================================== */

void
pdf_mstudent(double *y, double *x, int *nobs, int *vars, double *center, double *Scatter, double *df)
{ /* evaluation of the multivariate Student-t log-density function */
  int errcode = 0, job = 0, n = *nobs, p = *vars;
  double c_eta, eta = *df, half_eta_inv, *Root, *z, D2, log_pdf;

  Root = (double *) R_Calloc(p * p, double);
  z    = (double *) R_Calloc(p, double);

  copy_lower(Root, p, Scatter, p, p);
  chol_decomp(Root, p, p, job, &errcode);
  if (errcode)
    error("Cholesky decomposition in pdf_mstudent gave code %d", errcode);
  
  c_eta = eta / (1. - 2. * eta);
  half_eta_inv = .5 / eta;
  log_pdf  = .5 * p * log(c_eta * M_1_PI) + lgammafn(half_eta_inv + .5 * p) - lgammafn(half_eta_inv);
  log_pdf -= logAbsDet(Root, p, p);
  
  for (int i = 0; i < n; i++) {
    copy_vec(z, 1, x + i, n, p);
    D2 = mahalanobis(z, p, center, Root);
    y[i] = log_pdf - (half_eta_inv + .5 * p) * log1p(c_eta * D2);
  }

  R_Free(Root); R_Free(z);
}

/* ========================================================================== *
 * pdf for the multivariate normal distribution
 * ========================================================================== */

void
pdf_mnorm(double *y, double *x, int *nobs, int *vars, double *center, double *Scatter)
{ /* evaluation of the multivariate normal log-density function */
  int errcode = 0, job = 0, n = *nobs, p = *vars;
  double *Root, *z, D2, log_pdf = 0.0;

  Root = (double *) R_Calloc(p * p, double);
  z    = (double *) R_Calloc(p, double);

  copy_lower(Root, p, Scatter, p, p);
  chol_decomp(Root, p, p, job, &errcode);
  if (errcode)
    error("Cholesky decomposition in pdf_mlaplace gave code %d", errcode);
  
  log_pdf -= (double) p * M_LN_SQRT_2PI;
  log_pdf -= logAbsDet(Root, p, p);
  
  for (int i = 0; i < n; i++) {
    copy_vec(z, 1, x + i, n, p);
    D2 = mahalanobis(z, p, center, Root);
    y[i] = log_pdf - 0.5 * D2;
  }

  R_Free(Root); R_Free(z);
}
