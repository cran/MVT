/* ID: fitter_UN.c, last updated 2024-09-23, F.Osorio */

#include "base.h"
#include "interface.h"
#include "MVT.h"
#include "shared.h"

/* static functions.. */
static MODEL model_init(double *, int *, double *, int *, double *, double *, double *, double *, double *);
static void CS_init(MODEL, double *, double *);
static void HOMO_init(MODEL, double *, double *, double *);
static void model_free(MODEL);
/* ..end declarations */

void
fitter_dispatcher(double *y, int *pdims, double *settings, int *covType, double *center,
  double *Scatter, double *sigma2, double *rho, double *Phi, double *distances,
  double *weights, double *logLik, double *control)
{ /* estimation for the multivariate t-distribution using the EM algorithm */
  MODEL model;

  model = model_init(y, pdims, settings, covType, center, Scatter, distances, weights, control);
  switch (model->covType) {
    case UN: /* Unstructured */
      control[3] = (double) fitter_UN(model);
      break;
    case CS: /* Compound Symmetry (Equicorrelation) */
      CS_init(model, sigma2, rho);
      control[3] = (double) fitter_CS(model);
      break;
    case DIAGONAL: /* Diagonal */
      control[3] = (double) fitter_DIAG(model);
      break;
    case HOMOGENEITY: /* Variance homogeneity */
      HOMO_init(model, sigma2, Phi, control);
      control[3] = (double) fitter_HOMO(model);
      break;
    default:
      control[3] = (double) fitter_UN(model);
      break;
  }
  *logLik = log_Lik(model->family, model->dm, model->distances, model->Scatter);
  model_free(model);
}

static MODEL
model_init(double *y, int *pdims, double *settings, int *covType, double *center,
  double *Scatter, double *distances, double *weights, double *control)
{ /* constructor for a multivariate object */
  MODEL model;

  model = (MODEL) R_Calloc(1, MODEL_struct);
  model->dm = dims(pdims);
  model->settings = settings;
  model->family = family_init(settings);
  model->covType = *covType;
  model->y = y;
  model->center = center;
  model->Scatter = Scatter;
  model->distances = distances;
  model->weights = weights;
  model->control = control;
  model->maxIter = (int) control[0];
  model->tolerance = control[1];
  model->fixShape = (int) control[2];
  return model;
}

static void
CS_init(MODEL model, double *sigma2, double *rho)
{ /* initialization of equicorrelation fitter */
  model->sigma2 = sigma2;
  model->rho = rho;
}

static void
HOMO_init(MODEL model, double *sigma2, double *Phi, double *control)
{ /* initialization of variance homogeneity fitter */
  model->sigma2 = sigma2;
  model->Phi = Phi;
  model->both = (int) control[3];
}

static void
model_free(MODEL this)
{ /* destructor for a model object */
  dims_free(this->dm);
  family_free(this->family);
  R_Free(this);
}
