/* ID: family.c, last updated 2024-09-23, F.Osorio */

#include "MVT.h"

/* static functions.. */
static double weight_normal(void);
static double weight_student(double, double, double);
static double negQfnc_student(double, void *);
static void update_eta_student(DIMS, double *, double *, double);
static double logLik_normal(DIMS, double *);
static double logLik_student(DIMS, double, double *);
/* ..end declarations */

/* functions for dealing with 'family' objects */

FAMILY
family_init(double *settings) { 
  /* constructor for a family object */
  FAMILY ans;

  ans = (FAMILY) R_Calloc(1, FAMILY_struct);
  ans->kind = (int) settings[0];
  ans->eta  = settings + 1;
  return ans;
}

void
family_free(FAMILY this) { 
  /* destructor for a family object */
  R_Free(this);
}

/* functions for computation of weights */

static double
weight_normal(void) { 
  /* normal weight */
  return 1.;
}

static double
weight_student(double eta, double p, double distance) { 
  /* Student-t weight */
  return (1. + eta * p) / (1. + eta * (distance - 2.));
}

double
do_weight(FAMILY family, double p, double distance)
{ /* weights dispatcher */
  double wts;

  switch (family->kind) {
    case NORMAL:
      wts = weight_normal();
      break;
    case STUDENT:
      wts = weight_student((family->eta)[0], p, distance);
      break;
    default:
      wts = weight_normal();
      break;
  }
  return wts;
}

/* functions for the estimation of the 'shape' parameter */

static double
negQfnc_student(double eta, void *pars)
{ /* for Brent's procedure */
  QTpars st = (QTpars) pars;
  DIMS dm = st->dm;
  double accum = 0.0, df, half_eta_inv, half_c_inv, val;

  /* definitions */
  half_eta_inv = .5 / eta;
  half_c_inv   = half_eta_inv - 1.;

  /* sum of differences between log-weights and weights */
  for (int i = 0; i < dm->n; i++)
    accum += log((st->weights)[i]) - (st->weights)[i];
  accum /= (double) dm->n;

  /* compute Q-function for Student-t */
  df = .5 * (1. / st->eta + (double) dm->p);
  val  = half_c_inv * (digamma(df) - log(df) + accum);
  val += half_eta_inv * log(half_c_inv) - lgammafn(half_eta_inv);
  val *= (double) dm->n;
  st->Qfnc = val;

  return -val;
}

static void
update_eta_student(DIMS dm, double *eta, double *weights, double tol)
{
  QTpars pars;
  pars = (QTpars) R_Calloc(1, QT_pars);

  /* constructs a Q-function object */
  pars->dm = dm;
  pars->weights = weights;
  pars->eta = *eta;

  /* call optimizer */
  *eta = brent(.0, .5, negQfnc_student, pars, tol);

  R_Free(pars);
}

void
update_mixture(FAMILY family, DIMS dm, double *distances, double *weights, double tol)
{ /* update dispatcher */
  switch (family->kind) {
    case NORMAL:
      break;
    case STUDENT:
      update_eta_student(dm, family->eta, weights, tol);
      break;
    default:
      break;
  }
}

/*  functions for evaluation of the log-likelihood */

static double
logLik_normal(DIMS dm, double *distances)
{ /* gaussian log-likelihood */
  double accum = .0, val;

  /* sum of Mahalanobis distances */
  for (int i = 0; i < dm->n; i++)
    accum += *distances++;

  /* computation of the log-likelihood */
  val  = (double) (dm->n * dm->p) * M_LN_SQRT_2PI;
  val += .5 * accum;
  return -val;
}

static double
logLik_student(DIMS dm, double eta, double *distances)
{ /* Student-t log-likelihood */
  double accum = .0, c_eta, half_eta_inv, val;

  /* definitions */
  c_eta = eta / (1. - 2. * eta);
  half_eta_inv = .5 / eta;

  /* sum the kernel of the log-density */
  for (int i = 0; i < dm->n; i++)
    accum += log1p(c_eta * *distances++);

  /* computation of the log-likelihood */
  val  = .5 * dm->p * log(c_eta * M_1_PI);
  val += lgammafn(half_eta_inv + .5 * dm->p) - lgammafn(half_eta_inv);
  val *= (double) dm->n;
  val -= (half_eta_inv + .5 * dm->p) * accum;
  return val;
}

double
logLik_kernel(FAMILY family, DIMS dm, double *distances)
{ /* logLik dispatcher */
  double ans;

  switch (family->kind) {
    case NORMAL:
      ans = logLik_normal(dm, distances);
      break;
    case STUDENT:
      ans = logLik_student(dm, (family->eta)[0], distances);
      break;
    default:
      ans = logLik_normal(dm, distances);
      break;
  }
  return ans;
}
