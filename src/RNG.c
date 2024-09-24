/* ID: RNG.c, last updated 2024-09-23, F.Osorio */

#include "base.h"
#include "interface.h"

/* static functions.. */
static void std_student_rand(double *, double, int, int);
/* ..end declarations */

void
RNG_mstudent(double *y, int *pdims, double *center, double *Scatter, double *eta)
{ /* multivariate Student-t random generation */
  DIMS dm;
  char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
  int info = 0, job = 1;

  dm = dims(pdims);
  GetRNGstate();
  chol_decomp(Scatter, dm->p, dm->p, job, &info);
  if (info)
    error("DPOTRF in cholesky decomposition gave code %d", info);
  std_student_rand(y, *eta, dm->n, dm->p);
  mult_triangular_mat(1.0, Scatter, dm->p, dm->p, dm->n, side, uplo, trans, diag, y, dm->p);
  for (int i = 0; i < dm->n; i++) {
    ax_plus_y(1.0, center, 1, y, 1, dm->p);
    y += dm->p;
  }
  PutRNGstate();
  dims_free(dm);
}

static void
std_student_rand(double *y, double eta, int n, int p)
{ /* standard Student-t variates */
  double tau, radial, scale, shape;

  shape = .5 / eta;
  scale = 2. * eta / (1. - 2. * eta);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    tau = rgamma(shape, scale);
    radial = R_pow(tau, -.5);
    scale_vec(y, p, 1, radial);
    y += p;
  }
}
