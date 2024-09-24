/* ID: fitter_DIAG.c, last updated 2024-09-23, F.Osorio */

#include "base.h"
#include "interface.h"
#include "shared.h"
#include "MVT.h"

/* static functions.. */
static void update_diagonal(double *, int, int, double *, double *, double *);
/* ..end declarations */

int
fitter_DIAG(MODEL model)
{ /* fits the multivariate t-model considering a diagonal (DIAG) covariance matrix */
  DIMS dm = model->dm;
  int iter = 0;
  double conv, tol = R_pow(model->tolerance, 2. / 3.), RSS = (double) dm->n * dm->p, newRSS;

  /* main loop */
  repeat {
    /* E-step */
    E_step(model->y, dm->n, dm->p, model->center, model->Scatter, model->family, model->distances, model->weights);

    /* M-step */
    center_online(model->y, dm->n, dm->p, model->weights, model->center);
    update_diagonal(model->y, dm->n, dm->p, model->weights, model->center, model->Scatter);

    /* updating 'eta' parameter */
    if (!(model->fixShape)) {
      E_step(model->y, dm->n, dm->p, model->center, model->Scatter, model->family, model->distances, model->weights);
      update_mixture(model->family, model->dm, model->distances, model->weights, tol);
    }

    iter++;

    /* eval convergence */
    newRSS = dot_product(model->weights, 1, model->distances, 1, dm->n);
    conv = fabs((newRSS - RSS) / (newRSS + ABSTOL));
    if (conv < model->tolerance)
      break; /* successful completion */
    if (iter >= model->maxIter)
      break;  /* maximum number of iterations exceeded */
    RSS = newRSS;
  }
  return iter;
}

static void
update_diagonal(double *x, int n, int p, double *weights, double *center, double *Scatter)
{ /* update diagonal elements of Scatter matrix */
  double res, wts, *z, *sigma;

  /* initialization */
  sigma = (double *) R_Calloc(p, double);
  z     = (double *) R_Calloc(p, double);
  setzero(Scatter, p, p, p);

  /* updating stage */
  for (int i = 1; i < n; i++) {
    wts = weights[i];
    copy_vec(z, 1, x + i, n, p);
    for (int j = 0; j < p; j++) {
      res = z[j] - center[j];
      sigma[j] += wts * SQR(res);
    }
  }

  /* saving results, it is assumed that 'Scatter' matrix is zeroed */
  for (int j = 0; j < p; j++)
    Scatter[j * (p + 1)] = sigma[j] / (double) n;

  R_Free(sigma); R_Free(z);
}
