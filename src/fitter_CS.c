/* ID: fitter_CS.c, last updated 2021-03-12, F.Osorio */

#include "base.h"
#include "interface.h"
#include "shared.h"
#include "MVT.h"

/* static functions.. */
static void update_rho(double *, int, double *, double *);
static void update_sigma2(double *, int, double *);
/* ..end declarations */

int
fitter_CS(MODEL model)
{ /* fits the multivariate t-model considering a diagonal (DIAG) covariance matrix */
  DIMS dm = model->dm;
  int iter = 0;
  double conv, tol = R_pow(model->tolerance, 2. / 3.), RSS = (double) dm->n * dm->p, newRSS;

  /* main loop */
  repeat {
    /* E-step */
    E_step(model->y, dm->n, dm->p, model->center, model->Scatter, model->family, model->distances, model->weights);

    /* CM-steps */
    center_and_Scatter(model->y, dm->n, dm->p, model->weights, model->center, model->Scatter);
    update_sigma2(model->Scatter, dm->p, model->sigma2);
    update_rho(model->Scatter, dm->p, model->sigma2, model->rho);
    compSymm_pd(*(model->sigma2), *(model->rho), dm->p, model->Scatter);

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
update_sigma2(double *Scatter, int p, double *sigma2)
{ /* compute the sigma2 estimate under equicorrelation */
  *sigma2 = trace_mat(Scatter, p, p) / (double) p;
}

static void
update_rho(double *Scatter, int p, double *sigma2, double *rho)
{ /* compute the rho estimate under equicorrelation */
  int job = 0;
  double cor;

  cor = sum_lower_tri(Scatter, p, p, job); /* strictly lower part */
  *rho = 2. * cor / (*sigma2 * p * (p - 1.));
}
