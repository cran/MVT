#include "fitter.h"

/* declaration of static functions */

/* functions to deal with dims objects */
static DIMS dims(int *);
static void dims_free(DIMS);

/* routines for estimation in multivariate Student-t distributions */
static MODEL model_init(double *, int *, double *, double *, double *, double *, double *, double *);
static void model_free(MODEL);
static int EM_iterate(MODEL);
static void E_step(MODEL);
static double mahalanobis(double *, int, double *, double *);
static void update_center(MODEL);
static void update_Scatter(MODEL);

/* routine for estimation under restricted center */
static int restricted_center_iterate(MODEL);

/* routines for estimation under variance homogeneity */
static MODEL variance_homogeneity_init(double *, int *, double *, double *, double *, double *, double *, double *, double *, double *);
static void update_homogeneity_center(MODEL);
static void update_homogeneity_sigma2(MODEL);
static void update_homogeneity_Phi(MODEL);
static void cov2cor(double *, int);
static int variance_homogeneity_iterate(MODEL);

/* routines for estimation under equicorrelation */
static MODEL equicorrelation_init(double *, int *, double *, double *, double *, double *, double *, double *, double *, double *);
static void update_equicorrelation_sigma2(MODEL);
static void update_equicorrelation_rho(MODEL);
static void compSymm_pd(double, double, int, double *);
static int equicorrelation_iterate(MODEL);

/* routine for evaluation of marginal log-likelihood function */
static double log_Lik(FAMILY, DIMS, double *, double *);

void
EM_fit(double *y, int *pdims, double *settings, double *center, double *Scatter,
    double *distances, double *weights, double *logLik, double *control)
{   /* estimation for multivariate heavy-tailed distributions */
    MODEL model;

    model = model_init(y, pdims, settings, center, Scatter, distances, weights, control);
    control[3] = (double) EM_iterate(model);
    *logLik = log_Lik(model->family, model->dm, model->distances, model->Scatter);
    model_free(model);
}

static DIMS
dims(int *pdims)
{   /* dims object for multivariate models */
    DIMS ans;

    ans = (DIMS) Calloc(1, DIMS_struct);
    ans->N = (int) pdims[0];
    ans->n = ans->N;
    ans->p = (int) pdims[1];
    return ans;
}

static void
dims_free(DIMS this)
{   /* destructor for a dims object */
    Free(this);
}

static MODEL
model_init(double *y, int *pdims, double *settings, double *center, double *Scatter,
    double *distances, double *weights, double *control)
{   /* constructor for a multivariate object */
    MODEL model;

    model = (MODEL) Calloc(1, MODEL_struct);
    model->dm = dims(pdims);
    model->settings = settings;
    model->family = family_init(settings);
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
model_free(MODEL this)
{   /* destructor for a model object */
    dims_free(this->dm);
    family_free(this->family);
    Free(this);
}

static int
EM_iterate(MODEL model)
{
    int iter = 0;
    double conv, tol = model->tolerance, logLik, newlogLik;

    /* initialization */
    logLik = log_Lik(model->family, model->dm, model->distances, model->Scatter);

    /* main loop */
    repeat {
        /* E-step */
        E_step(model);
        
        /* M-step */
        update_center(model);
        update_Scatter(model);
        if (!(model->fixShape))
            update_mixture(model->family, model->dm, model->distances, model->weights, tol);

        /* eval convergence */
        iter++;
        newlogLik = log_Lik(model->family, model->dm, model->distances, model->Scatter);
        conv = fabs((newlogLik - logLik) / (newlogLik + EPS_CONV));
        if (conv < model->tolerance)
            break;  /* successful completion */
        if (iter >= model->maxIter)
            break;  /* maximum number of iterations exceeded */
        logLik = newlogLik;
    }
    return iter;
}

static void
E_step(MODEL model)
{
    DIMS dm = model->dm;
    double *Root;
    int i, info = 0, job = 0;

    Root = (double *) Calloc(SQR(dm->p), double);
    copy_mat(Root, dm->p, model->Scatter, dm->p, dm->p, dm->p);
    chol_decomp(Root, dm->p, dm->p, job, &info);
    if (info)
        error("chol_decomp in E_step gave code %d", info);
    for (i = 0; i < dm->n; i++) {
        (model->distances)[i] = mahalanobis(model->y + i * dm->p, dm->p, model->center, Root);
        (model->weights)[i] = do_weight(model->family, (double) dm->p, (model->distances)[i]);
    }
    Free(Root);
}

static double
mahalanobis(double *y, int p, double *center, double *Root)
{   /* Mahalanobis distances */
    double ans, minus = -1., *z;
    int inc = 1, info = 0, job = 0;

    z = (double *) Calloc(p, double);
    Memcpy(z, y, p);
    F77_CALL(daxpy)(&p, &minus, center, &inc, z, &inc);
    backsolve(Root, p, p, z, p, 1, job, &info);
    if (info)
        error("backsolve in mahalanobis gave code %d", info);
    ans = norm_sqr(z, p, 1);
    Free(z);
    return ans;
}

static void
update_center(MODEL model)
{   /* compute the center estimate */
    DIMS dm = model->dm;
    double factor = 1., wts, *center;
    int i, inc = 1;

    center = (double *) Calloc(dm->p, double);
    for (i = 0; i < dm->n; i++) {
        wts = (model->weights)[i];
        F77_CALL(daxpy)(&(dm->p), &wts, model->y + i * dm->p, &inc, center, &inc);
    }
    factor /= F77_CALL(dasum)(&(dm->n), model->weights, &inc);
    F77_CALL(dscal)(&(dm->p), &factor, center, &inc);
    Memcpy(model->center, center, dm->p);
    Free(center);
}

static void
update_Scatter(MODEL model)
{   /* update the Scatter matrix estimate */
    DIMS dm = model->dm;
    double minus = -1., wts, *Scatter, *z;
    int i, inc = 1;

    Scatter = (double *) Calloc(SQR(dm->p), double);
    z = (double *) Calloc(dm->p, double);
    for (i = 0; i < dm->n; i++) {
        wts = (model->weights)[i];
        Memcpy(z, model->y + i * dm->p, dm->p);
        F77_CALL(daxpy)(&(dm->p), &minus, model->center, &inc, z, &inc);
        rank1_update(Scatter, dm->p, dm->p, dm->p, wts, z, z);
    }
    scale_mat(model->Scatter, dm->p, Scatter, dm->p, dm->p, dm->p, 1. / dm->n);
    Free(z); Free(Scatter);
}

void
restricted_center(double *y, int *pdims, double *settings, double *center,
    double *Scatter, double *distances, double *weights, double *logLik,
    double *control)
{   /* estimation for multivariate heavy-tailed distributions */
    MODEL model;

    model = model_init(y, pdims, settings, center, Scatter, distances, weights, control);
    control[3] = (double) restricted_center_iterate(model);
    *logLik = log_Lik(model->family, model->dm, model->distances, model->Scatter);
    model_free(model);
}

static int
restricted_center_iterate(MODEL model)
{
    int iter = 0;
    double conv, tol = model->tolerance, logLik, newlogLik;

    /* initialization */
    logLik = log_Lik(model->family, model->dm, model->distances, model->Scatter);

    /* main loop */
    repeat {
        /* E-step */
        E_step(model);
        
        /* M-step */
        update_Scatter(model);
        if (!(model->fixShape))
            update_mixture(model->family, model->dm, model->distances, model->weights, tol);

        /* eval convergence */
        iter++;
        newlogLik = log_Lik(model->family, model->dm, model->distances, model->Scatter);
        conv = fabs((newlogLik - logLik) / (newlogLik + EPS_CONV));
        if (conv < model->tolerance)
            break; /* successful completion */
        if (iter >= model->maxIter)
            break; /* maximum number of iterations exceeded */
        logLik = newlogLik;
    }
    return iter;
}

void
variance_homogeneity(double *y, int *pdims, double *settings, double *center, double *Scatter,
    double *sigma2, double *Phi, double *distances, double *weights, double *logLik, double *control)
{   /* ML estimation under variance homogeneity */
    MODEL model;

    /* initialization */
    model = variance_homogeneity_init(y, pdims, settings, center, Scatter, sigma2, Phi, distances, weights, control);
    
    /* ML estimation and log-likelihood evaluation */
    control[4] = (double) variance_homogeneity_iterate(model);
    *logLik = log_Lik(model->family, model->dm, model->distances, model->Scatter);
    
    model_free(model);
}

static MODEL
variance_homogeneity_init(double *y, int *pdims, double *settings, double *center, double *Scatter,
    double *sigma2, double *Phi, double *distances, double *weights, double *control)
{   /* constructor for a multivariate object */
    MODEL model;

    model = (MODEL) Calloc(1, MODEL_struct);
    model->dm = dims(pdims);
    model->settings = settings;
    model->family = family_init(settings);
    model->y = y;
    model->center = center;
    model->Scatter = Scatter;
    model->sigma2 = sigma2;
    model->Phi = Phi;
    model->distances = distances;
    model->weights = weights;
    model->control = control;
    model->maxIter = (int) control[0];
    model->tolerance = control[1];
    model->fixShape = (int) control[2];
    model->both = (int) control[3];
    return model;
}

static void
update_homogeneity_center(MODEL model)
{   /* compute the 'common' center estimate */
    DIMS dm = model->dm;
    int i, j;
    double accum, delta, prod = 0.0, total = 0.0;

    for (j = 0; j < dm->p; j++) {
        accum = 0.0;
        for (i = 0; i < dm->p; i++)
            accum += (model->Phi)[i + j * dm->p];
        prod  += accum * (model->center)[j];
        total += accum;
    }
    delta = prod / total;
    
    for (i = 0; i < dm->p; i++)
        (model->center)[i] = delta;
}

static void
update_homogeneity_sigma2(MODEL model)
{   /* compute the sigma2 estimate under variance homogeneity */
    DIMS dm = model->dm;
    int inc = 1;
    double pp = SQR(dm->p), tr;

    tr = dot_product(model->Phi, inc, model->Scatter, inc, pp);
    (model->sigma2)[0] = tr / dm->p;
}

static void
update_homogeneity_Phi(MODEL model)
{   /* compute the Phi estimate estimate under variance homogeneity */
    DIMS dm = model->dm;
    double scale = 1.;

    scale /= (model->sigma2)[0];
    scale_mat(model->Phi, dm->p, model->Scatter, dm->p, dm->p, dm->p, scale);
}

static void
cov2cor(double *cov, int p)
{   /* scales a 'covariance' matrix into the corresponding correlation matrix */
    int i, j;
    double *s;

    s = (double *) Calloc(p, double);
    for (i = 0; i < p; i++)
        s[i] = cov[i * (p + 1)];
    for (i = 0; i < p; i++) {
        cov[i * (p + 1)] = 1.0;
        for (j = i + 1; j < p; j++) {
            *(cov + i + j * p) = *(cov + i + j * p) / sqrt(s[i] * s[j]);
            *(cov + j + i * p) = *(cov + i + j * p);
        }
    }
    Free(s);
}

static int
variance_homogeneity_iterate(MODEL model)
{
    int iter = 0, info = 0;
    double conv, tol = model->tolerance, logLik, newlogLik;

    /* initialization */
    logLik  = log_Lik(model->family, model->dm, model->distances, model->Scatter);

    /* main loop */
    repeat {
        /* E-step */
        E_step(model);
        
        /* CM-steps */
        invert_mat(model->Phi, model->dm->p, model->dm->p, &info);
        if (info)
            error("invert_mat in variance_homogeneity_iterate gave code %d", info);
        update_center(model);
        if (model->both)
            update_homogeneity_center(model);
        update_Scatter(model);
        update_homogeneity_sigma2(model);
        update_homogeneity_Phi(model);
        cov2cor(model->Phi, model->dm->p);
        scale_mat(model->Scatter, model->dm->p, model->Phi, model->dm->p, model->dm->p, model->dm->p, (model->sigma2)[0]);
        if (!(model->fixShape))
            update_mixture(model->family, model->dm, model->distances, model->weights, tol);

        /* eval convergence */
        iter++;
        newlogLik = log_Lik(model->family, model->dm, model->distances, model->Scatter);
        conv = fabs((newlogLik - logLik) / (newlogLik + EPS_CONV));
        if (conv < model->tolerance)
            break; /* successful completion */
        if (iter >= model->maxIter)
            break; /* maximum number of iterations exceeded */
        logLik = newlogLik;
    }
    return iter;
}

void
equicorrelation_fit(double *y, int *pdims, double *settings, double *center, double *Scatter,
    double *sigma2, double *rho, double *distances, double *weights, double *logLik, double *control)
{   /* ML estimation under equicorrelation matrix */
    MODEL model;

    /* initialization */
    model = equicorrelation_init(y, pdims, settings, center, Scatter, sigma2, rho, distances, weights, control);
    
    /* ML estimation and log-likelihood evaluation */
    control[3] = (double) equicorrelation_iterate(model);
    *logLik = log_Lik(model->family, model->dm, model->distances, model->Scatter);

    model_free(model);
}

static MODEL
equicorrelation_init(double *y, int *pdims, double *settings, double *center, double *Scatter,
    double *sigma2, double *rho, double *distances, double *weights, double *control)
{   /* constructor for a multivariate object */
    MODEL model;

    model = (MODEL) Calloc(1, MODEL_struct);
    model->dm = dims(pdims);
    model->settings = settings;
    model->family = family_init(settings);
    model->y = y;
    model->center = center;
    model->Scatter = Scatter;
    model->sigma2 = sigma2;
    model->rho = rho;
    model->distances = distances;
    model->weights = weights;
    model->control = control;
    model->maxIter = (int) control[0];
    model->tolerance = control[1];
    model->fixShape = (int) control[2];
    return model;
}

static void
update_equicorrelation_sigma2(MODEL model)
{   /* compute the sigma2 estimate under equicorrelation */
    DIMS dm = model->dm;

    (model->sigma2)[0] = trace_mat(model->Scatter, dm->p, dm->p) / dm->p;
}

static void
update_equicorrelation_rho(MODEL model)
{   /* compute the rho estimate under equicorrelation */
    DIMS dm = model->dm;
    double rho;

    rho = sum_lower_tri(model->Scatter, dm->p, dm->p);
    (model->rho)[0] = 2. * rho / ((model->sigma2)[0] * dm->p * (dm->p - 1));
}

static void
compSymm_pd(double sigma2, double rho, int p, double *mat)
{   /* construcs the equicorrelation matrix (compound symmetry) */
    int i, j;
    
    for (i = 0; i < p; i++) {
        mat[i * (p + 1)] = sigma2;
        for (j = i + 1; j < p; j++)
            *(mat + i + j * p) = *(mat + j + i * p) = sigma2 * rho;
    }
}

static int
equicorrelation_iterate(MODEL model)
{
    int iter = 0;
    double conv, tol = model->tolerance, logLik, newlogLik;

    /* initialization */
    logLik = log_Lik(model->family, model->dm, model->distances, model->Scatter);

    /* main loop */
    repeat {
        /* E-step */
        E_step(model);
        
        /* CM-steps */
        update_center(model);
        update_Scatter(model);
        update_equicorrelation_sigma2(model);
        update_equicorrelation_rho(model);
        compSymm_pd((model->sigma2)[0], (model->rho)[0], model->dm->p, model->Scatter);
        if (!(model->fixShape))
            update_mixture(model->family, model->dm, model->distances, model->weights, tol);

        /* eval convergence */
        iter++;
        newlogLik = log_Lik(model->family, model->dm, model->distances, model->Scatter);
        conv = fabs((newlogLik - logLik) / (newlogLik + EPS_CONV));
        if (conv < model->tolerance)
            break; /* successful completion */
        if (iter >= model->maxIter)
            break; /* maximum number of iterations exceeded */
        logLik = newlogLik;
    }
    return iter;
}

static double
log_Lik(FAMILY family, DIMS dm, double *distances, double *Scatter)
{   /* evaluate the log-likelihood function */
    double ans = 0., *Root;
    int i, info = 0, job = 0;

    Root = (double *) Calloc(SQR(dm->p), double);
    copy_mat(Root, dm->p, Scatter, dm->p, dm->p, dm->p);
    chol_decomp(Root, dm->p, dm->p, job, &info);
    if (info)
        error("chol_decomp in log_Lik gave code %d", info);
    ans += logLik_kernel(family, dm, distances);
    ans -= dm->n * logAbsDet(Root, dm->p, dm->p);
    Free(Root);
    return ans;
}
