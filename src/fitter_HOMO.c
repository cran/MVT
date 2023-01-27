/* ID: fitter_HOMO.c, last updated 2021-03-12, F.Osorio */

#include "base.h"
#include "interface.h"
#include "shared.h"
#include "MVT.h"

/* static functions.. */
static void update_common_center(double *, int, double *);
static void update_cor(double *, int, double *, double *);
static void update_sigma2(double *, double *, int, double *);
/* ..end declarations */

int
fitter_HOMO(MODEL model)
{ /* fits the multivariate t-model considering a diagonal (DIAG) covariance matrix */
  DIMS dm = model->dm;
  int iter = 0, info = 0;
  double conv, tol = R_pow(model->tolerance, 2. / 3.), RSS = (double) dm->n * dm->p, newRSS;

  /* main loop */
  repeat {
    /* E-step */
    E_step(model->y, dm->n, dm->p, model->center, model->Scatter, model->family, model->distances, model->weights);

    /* CM-steps */
    invert_mat(model->Phi, dm->p, dm->p, &info);
    if (info)
      error("matrix inversion in fitter_HOMO gave code %d", info);
    center_and_Scatter(model->y, dm->n, dm->p, model->weights, model->center, model->Scatter);
    if (model->both)
      update_common_center(model->Phi, dm->p, model->center);
    update_sigma2(model->Scatter, model->Phi, dm->p, model->sigma2);
    update_cor(model->Scatter, dm->p, model->sigma2, model->Phi);
    cov2cor(model->Phi, dm->p);
    scale_mat(model->Scatter, dm->p, *(model->sigma2), model->Phi, dm->p, dm->p, dm->p);

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
update_common_center(double *R, int p, double *center)
{ /* compute the 'common' center estimate */
  double accum, delta, prod = 0.0, total = 0.0;

  for (int j = 0; j < p; j++) {
    accum = 0.0;
    for (int i = 0; i < p; i++)
      accum += R[i + j * p];
    prod  += accum * center[j];
    total += accum;
  }
  delta = prod / total;

  for (int j = 0; j < p; j++)
    center[j] = delta;
}

static void
update_sigma2(double *Scatter, double *R, int p, double *sigma2)
{ /* compute the sigma2 estimate under variance homogeneity */
  double tr;

  tr = dot_product(R, 1, Scatter, 1, p * p);
  *sigma2 = tr / (double) p;
}

static void
update_cor(double *Scatter, int p, double *sigma2, double *R)
{ /* compute the Phi estimate under variance homogeneity */
  double scale = 1.;

  scale /= *sigma2;
  scale_mat(R, p, scale, Scatter, p, p, p);
}
