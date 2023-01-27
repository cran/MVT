/* ID: MVT.h, last updated 2021-03-12, F.Osorio */

#ifndef MVT_H
#define MVT_H

#include "base.h"

/* available covariance structures */
typedef enum {
  UN,
  CS,
  DIAGONAL,
  HOMOGENEITY
} covStruc;

/* available families */
typedef enum {
  NORMAL,
  STUDENT
} classes;

/* Student-t family structure */
typedef struct FAMILY_struct {
  classes kind;   /* family kind */
  double *eta;    /* shape parameter */
} FAMILY_struct, *FAMILY;

/* Q-function info required for shape parameter estimation */
typedef struct QT_pars {
  DIMS dm;
  double eta, Qfnc;
  double *lengths, *weights;
} QT_pars, *QTpars;

/* structure to hold model results */
typedef struct MODEL_struct {
  DIMS      dm;       /* dimension data info */
  FAMILY    family;   /* family data and info */
  covStruc  covType;  /* covariance structure */
  int
    *pdims;           /* dimensions */
  double
    *y,               /* data matrix */
    *settings,        /* settings */
    *center,          /* position parameter estimates */
    *delta,           /* common center estimate */
    *Scatter,         /* scatter matrix estimate */
    *Phi,             /* correlation matrix */
    *sigma2,          /* scale estimate */
    *rho,             /* equicorrelation parameter estimate */
    *distances,       /* Mahalanobis distances */
    *weights,         /* weights for heavy tailed distributions */
    *control;         /* control settings for estimation algorithm */
  int
    maxIter,          /* maximun number of iterations */
    fixShape,         /* must estimate shape parameters? */
    both;             /* info about the homogeneity test */
  double
    tolerance;        /* convergence tolerance */
} MODEL_struct, *MODEL;

/* fitter for different covariance structures */
extern int fitter_CS(MODEL);
extern int fitter_DIAG(MODEL);
extern int fitter_HOMO(MODEL);
extern int fitter_UN(MODEL);

/* functions for dealing with 'family' objects */
extern FAMILY family_init(double *);
extern void family_free(FAMILY);

/* routines for computation of weights */
extern double do_weight(FAMILY, double, double);
extern void update_mixture(FAMILY, DIMS, double *, double *, double);

/* evaluation of the log-likelihood */
extern double logLik_kernel(FAMILY, DIMS, double *);

/* Brent's method for unidimensional minimization */
extern double brent(double, double, double (*f)(double, void *), void *, double);

#endif /* MVT_H */
