#include "family.h"

/* declaration of static functions */

/* functions for computation of weights */
static double weight_normal();
static double weight_student(double, double, double);

/* functions for update the parameters of mixture variables */
static double negQfnc_student(double, void *);
static void update_eta_student(DIMS, double *, double *, double);

/*  functions for evaluation of the log-likelihood */
static double logLik_normal(DIMS, double *);
static double logLik_student(DIMS, double, double *);

/* functions for dealing with 'family' objects */

FAMILY
family_init(double *settings)
{   /* constructor for a family object */
    FAMILY ans;

    ans = (FAMILY) Calloc(1, FAMILY_struct);
    ans->kind = (int) settings[0];
    ans->eta  = settings + 1;
    return ans;
}

void
family_free(FAMILY this)
{   /* destructor for a family object*/
    Free(this);
}

/* functions for computation of weights */

static double
weight_normal()
{   /* normal weight */
    return 1.;
}

static double
weight_student(double eta, double p, double distance)
{   /* Student-t weight */
    double c_inv, eta_inv = 1., wts;

    eta_inv /= eta;
    c_inv = eta_inv - 2.;
    wts = (eta_inv + p) / (c_inv + distance);
    return wts;
}

double
do_weight(FAMILY family, double p, double distance)
{   /* weights dispatcher */
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
{   /* for brent procedure */
    QTpars st = (QTpars) pars;
    DIMS dm = st->dm;
    double accum = .0, df, c_inv, eta_inv = 1., val;
    int i;

    eta_inv /= eta;
    c_inv = eta_inv - 2.;
    for (i = 0; i < dm->n; i++) 
        accum += log((st->weights)[i]) - (st->weights)[i];
    accum /= dm->n;

    /* compute Q-function for Student-t */
    df = 1. / st->eta + (double) dm->p;
    val  = .5 * c_inv * (digamma(.5 * df) - log(.5 * df));
    val += .5 * eta_inv * log(.5 * c_inv) - lgammafn(.5 * eta_inv);
    val += .5 * c_inv * accum;
    val *= dm->n;
    st->Qfnc = val;

    return -val;
}

static void
update_eta_student(DIMS dm, double *eta, double *weights, double tol)
{
    QTpars pars;
    
    pars = (QTpars) Calloc(1, QT_pars);

    /* constructs a Q-function object */
    pars->dm = dm;
    pars->weights = weights;
    pars->eta = *eta;
    
    /* call optimizer */
    *eta = brent(.0, .5, negQfnc_student, pars, tol);
    
    Free(pars);
}

void
update_mixture(FAMILY family, DIMS dm, double *distances, double *weights, double tol)
{   /* update dispatcher */
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
{   /* gaussian log-likelihood */
    int i;
    double accum = .0, val;

    val  = M_LN_SQRT_2PI;
    val *= (double) dm->n * dm->p;
    for (i = 0; i < dm->n; i++)
        accum += *distances++;
    val += .5 * accum;
    return -val;
}

static double
logLik_student(DIMS dm, double eta, double *distances)
{   /* Student-t log-likelihood */
    int i;
    double accum = .0, c_eta, eta_inv = 1., val;

    eta_inv /= eta;
    c_eta = eta / (1. - 2. * eta);
    val  = .5 * log(c_eta) - M_LN_SQRT_PI;
    val *= (double) dm->p;
    val += lgammafn(.5 * (eta_inv + dm->p)) - lgammafn(.5 * eta_inv);
    val *= (double) dm->n;
    for (i = 0; i < dm->n; i++)
        accum += log1p(c_eta * *distances++);
    val -= .5 * (eta_inv + dm->p) * accum;
    return val;
}

double
logLik_kernel(FAMILY family, DIMS dm, double *distances)
{   /* logLik dispatcher */
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
