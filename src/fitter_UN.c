/* ID: fitter_UN.c, last updated 2021-03-03, F.Osorio */

#include "base.h"
#include "interface.h"
#include "shared.h"
#include "MVT.h"

int
fitter_UN(MODEL model)
{ /* fits the multivariate t-model considering an unstructured (UN) covariance matrix */
  DIMS dm = model->dm;
  int iter = 0;
  double conv, tol = R_pow(model->tolerance, 2. / 3.), RSS = (double) dm->n * dm->p, newRSS;

  /* main loop */
  repeat {
    /* E-step */
    E_step(model->y, dm->n, dm->p, model->center, model->Scatter, model->family, model->distances, model->weights);

    /* M-step */
    center_and_Scatter(model->y, dm->n, dm->p, model->weights, model->center, model->Scatter);

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
