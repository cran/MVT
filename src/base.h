#ifndef BASE_H
#define BASE_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>

/* some definitions */
#define NULLP    (void *) 0
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define SQR(x)   R_pow_di(x, 2)
#define EPS_CONV 1.0e-2
#define GOLDEN   0.3819660112501051
#define repeat   for(;;)

/* dims structure */
typedef struct DIMS_struct {
    int
      N,        /* total number of observations */
      n,        /* number of observations */
      p;        /* number of variables */
} DIMS_struct, *DIMS;

/* QR structure */
typedef struct QR_struct {
    double *mat, *qraux;
    int ldmat, nrow, ncol;
} QR_struct, *QRStruct;

#endif /* BASE_H */
