#ifndef FITTER_H
#define FITTER_H

#include "base.h"
#include "matrix.h"
#include "family.h"
#include "random.h"

/* structure to hold model results */
typedef struct MODEL_struct {
    DIMS dm;        /* dimension data info */
    FAMILY family;  /* family data and info */
    int
      *pdims;       /* dimensions */
    double
      *y,           /* data matrix */
      *settings,    /* settings */
      *center,      /* position parameter estimates */
      *delta,       /* common center estimate */
      *Scatter,     /* scatter matrix estimate */
      *Phi,         /* correlation matrix */
      *sigma2,      /* scale estimate */
      *rho,         /* equicorrelation parameter estimate */
      *distances,   /* Mahalanobis distances */
      *weights,     /* weights for heavy tailed distributions */
      *control;     /* control settings for estimation algorithm */
    int
      maxIter,      /* maximun number of iterations */
      fixShape,     /* must estimate shape parameters? */
      both;         /* info about the homogeneity test */
    double
      tolerance;    /* convergence tolerance */
} MODEL_struct, *MODEL;

/* estimation in linear models under heavy tailed distributions */
extern void EM_fit(double *, int *, double *, double *, double *, double *, double *, double *, double *);
extern void equicorrelation_fit(double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void restricted_center(double *, int *, double *, double *, double *, double *, double *, double *, double *);
extern void variance_homogeneity(double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

#endif /* FITTER_H */
