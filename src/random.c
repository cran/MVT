#include "random.h"

/* declaration of static functions */

/* functions to deal with dims objects */
static DIMS dims(int *);
static void dims_free(DIMS);

/* spherical random generation */
static void rand_spherical_student(double *, double, int, int);
static void rand_spherical_norm(double *, int, int);

/* 'dims' functions */

static DIMS
dims(int *pdims)
{   /* dims object */
    DIMS ans;

    ans = (DIMS) Calloc(1, DIMS_struct);
    ans->n = (int) pdims[0];
    ans->p = (int) pdims[1];
    return ans;
}

static void
dims_free(DIMS this)
{   /* destructor for a dims object */
    Free(this);
}

/* multivariate normal random generation */

void
rand_norm(double *y, int *pdims, double *center, double *Scatter)
{   /* multivariate normal random generation */
    DIMS dm;
    char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
    double one = 1.;
    int i, inc = 1, info = 0, job = 1;
    
    dm = dims(pdims);
    GetRNGstate();
    chol_decomp(Scatter, dm->p, dm->p, job, &info);
    if (info)
	error("DPOTRF in cholesky decomposition gave code %d", info);
    rand_spherical_norm(y, dm->n, dm->p);
    F77_CALL(dtrmm)(side, uplo, trans, diag, &(dm->p), &(dm->n), &one, Scatter,
		    &(dm->p), y, &(dm->p));
    for (i = 0; i < dm->n; i++) {
	F77_CALL(daxpy)(&(dm->p), &one, center, &inc, y, &inc);
	y += dm->p;
    }
    PutRNGstate();
    dims_free(dm);
}

static void
rand_spherical_norm(double *y, int n, int p)
{   /* independent standard normal variates */
    int i, j;

    for (i = 0; i < n; i++) {
	for (j = 0; j < p; j++)
	    y[j] = norm_rand();
	y += p;
    }
}

/* multivariate Student-t random generation */

void
rand_student(double *y, int *pdims, double *center, double *Scatter, double *eta)
{   /* multivariate Student-t random generation */
    DIMS dm;
    char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
    double one = 1.;
    int i, inc = 1, info = 0, job = 1;
    
    dm = dims(pdims);
    GetRNGstate();
    chol_decomp(Scatter, dm->p, dm->p, job, &info);
    if (info)
	error("DPOTRF in cholesky decomposition gave code %d", info);
    rand_spherical_student(y, *eta, dm->n, dm->p);
    F77_CALL(dtrmm)(side, uplo, trans, diag, &(dm->p), &(dm->n), &one, Scatter,
		    &(dm->p), y, &(dm->p));
    for (i = 0; i < dm->n; i++) {
	F77_CALL(daxpy)(&(dm->p), &one, center, &inc, y, &inc);
	y += dm->p;
    }
    PutRNGstate();
    dims_free(dm);
}

static void
rand_spherical_student(double *y, double eta, int n, int p)
{   /* standard Student-t variates */
    int i, j, inc = 1;
    double tau, radial, scale, shape;

    shape = .5 / eta;
    scale = 2. * eta / (1. - 2. * eta);

    for (i = 0; i < n; i++) {
	for (j = 0; j < p; j++)
	    y[j] = norm_rand();
	tau = rgamma(shape, scale);
	radial = R_pow(tau, -.5);
	F77_CALL(dscal)(&p, &radial, y, &inc);
	y += p;
    }
}
