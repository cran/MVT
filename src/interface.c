/* ID: interface.c, last updated 2024-09-23, F.Osorio */

#include "base.h"
#include "interface.h"
#include <fastmatrix_API.h>

/* ========================================================================== *
 * basic matrix manipulations
 * ========================================================================== */

void
ax_plus_y(double alpha, double *x, int incx, double *y, int incy, int n)
{ /* y <- alpha * x + y */
  BLAS1_axpy(alpha, x, incx, y, incy, n);
}

double
dot_product(double *x, int incx, double *y, int incy, int n)
{ /* sum(x * y) */
  return BLAS1_dot_product(x, incx, y, incy, n);
}

void
scale_vec(double *x, int n, int inc, double alpha)
{ /* x <- alpha * x */
  BLAS1_scale(alpha, x, inc, n);
}

void
copy_vec(double *y, int incy, double *x, int incx, int n)
{ /* y <- x (alternative to Memcpy with increments not equal to 1) */
  BLAS1_copy(y, incy, x, incx, n);
}

void
copy_mat(double *y, int ldy, double *x, int ldx, int nrow, int ncol)
{ /* y <- x[,] */
  FM_copy_mat(y, ldy, x, ldx, nrow, ncol);
}

void
copy_lower(double *y, int ldy, double *x, int ldx, int p)
{ /* lower.tri(y) <- lower.tri(x) */
  FM_cpy_lower(x, ldx, p, y, ldy);
}

void
scale_mat(double *y, int ldy, double alpha, double *x, int ldx, int nrow, int ncol)
{ /* y <- alpha * x[,] */
  FM_scale_mat(y, ldy, alpha, x, ldx, nrow, ncol);
}

void
setzero(double *y, int ldy, int nrow, int ncol)
{ /* y[,] <- 0 */
  FM_setzero(y, ldy, nrow, ncol);
}

void
GAXPY(double *y, double alpha, double *a, int lda, int nrow, int ncol, double *x, double beta)
{ /* y <- alpha * a %*% x + beta * y */
  FM_GAXPY(y, alpha, a, lda, nrow, ncol, x, beta, 0); /* job = 0 */
}

void
mult_triangular_vec(double *a, int lda, int n, char *uplo, char *trans, char *diag, double *x, int inc)
{ /* wrapper to BLAS-2 routine 'DTRMV' */
  BLAS2_trmv(a, lda, n, uplo, trans, diag, x, inc);
}

void
mult_triangular_mat(double alpha, double *a, int lda, int nrow, int ncol, char *side, char *uplo, char *trans, char *diag, double *y, int ldy)
{ /* wrapper to BLAS-3 routine 'DTRMM' */
  BLAS3_trmm(alpha, a, lda, nrow, ncol, side, uplo, trans, diag, y, ldy);
}

double
trace_mat(double *a, int lda, int n)
{ /* sum(diag(a)) */
  return FM_trace(a, lda, n);
}

double
sum_lower_tri(double *a, int lda, int p, int job)
{ /* sum(lower.tri(a)): diag = FALSE, job = 0,
   *                    diag = TRUE,  job = 1 */
  return FM_sum_lower_tri(a, lda, p, job);
}

double
logAbsDet(double *a, int lda, int n)
{ /* log(abs(det(a))), where 'a' is a triangular matrix */
  return FM_logAbsDet(a, lda, n);
}

/* ========================================================================== *
 * routines for matrix decompositions
 * ========================================================================== */

void
chol_decomp(double *a, int lda, int p, int job, int *info)
{ /* cholesky factorization of a real symmetric positive definite matrix a.
   * the factorization has the form:
   * a <- l %*% t(l) (job = 0), or a <- t(u) %*% u (job = 1),
   * where u is an upper triangular matrix and l is lower triangular */
  FM_chol_decomp(a, lda, p, job, info);
}

/* ========================================================================== *
 * matrix inversion and triangular solver
 * ========================================================================== */

void
invert_mat(double *a, int lda, int n, int *info)
{ /* performs matrix inversion */
  FM_invert_mat(a, lda, n, info);
}

void
backsolve(double *r, int ldr, int n, double *b, int ldb, int nrhs, int *info)
{ /* backsolve solve triangular systems of the form r %*% x = b, where r
   * is an upper triangular matrix and b is a matrix containing the right-hand
   * sides to equations */
  FM_backsolve(r, ldr, n, b, ldb, nrhs, info);
}

/* ========================================================================== *
 * (multivariate) descriptive statistics 
 * ========================================================================== */

void
center_and_Scatter(double *x, int n, int p, double *weights, double *center, double *Scatter)
{ /* compute center and Scatter estimates using an online algorithm */
  FM_center_and_Scatter(x, n, p, weights, center, Scatter);
}

void
center_online(double *x, int n, int p, double *weights, double *center)
{ /* compute center estimate using an online algorithm */
  FM_online_center(x, n, p, weights, center);
}

/* ========================================================================== *
 * Mahalanobis distance
 * ========================================================================== */

double
mahalanobis(double *x, int p, double *center, double *Root)
{ /* Mahalanobis distance, argument 'Root' is the lower triangular factor of the
   * covariance matrix */
  return FM_mahalanobis(x, p, center, Root);
}

/* ========================================================================== *
 * Wilson-Hilferty transformation
 * ========================================================================== */

void
WH_chisq(double *distances, int n, int p, double *z)
{ /* Wilson-Hilferty transformation for chi-squared variables */
  FM_WH_chisq(distances, n, p, z);
}

void
WH_F(double *distances, int n, int p, double eta, double *z)
{ /* Wilson-Hilferty transformation for F deviates */
  FM_WH_F(distances, n, p, eta, z);
}

/* ========================================================================== *
 * Utilities
 * ========================================================================== */

void
cov2cor(double *cov, int p)
{ /* scales a 'covariance' matrix into the corresponding correlation matrix */
  FM_cov2cor(cov, p);
}
