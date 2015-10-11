#ifndef FAMILY_H
#define FAMILY_H

#include "base.h"
#include "optim.h"
#include "random.h"

/* available families */
typedef enum {
    NORMAL,
    STUDENT
} classes;

/* heavy tailed family structure */
typedef struct FAMILY_struct {
    classes kind;   /* family kind */
    double *eta;    /* shape parameter */
} FAMILY_struct, *FAMILY;

/* Q-function info required for the degrees of freedom estimation */
typedef struct QT_pars {
    DIMS dm;
    double eta, Qfnc;
    double *lengths, *weights;
} QT_pars, *QTpars;

/* functions for dealing with 'family' objects */
extern FAMILY family_init(double *);
extern void family_free(FAMILY);

/* routines for computation of weights */
extern double do_weight(FAMILY, double, double);
extern void update_mixture(FAMILY, DIMS, double *, double *, double);

/*  functions for evaluation of the log-likelihood */
extern double logLik_kernel(FAMILY, DIMS, double *);

#endif /* FAMILY_H */
