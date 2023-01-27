/* ID: base.h, last updated 2023-01-24, F.Osorio */

#ifndef BASE_H
#define BASE_H

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <R.h>
#include <Rconfig.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Print.h>

/* some definitions */
#define ABSTOL    1.0e-2
#define GOLDEN    0.3819660112501051
#define DNULLP    (double *) 0
#define NULLP     (void *) 0
#define REPORT    5
#define CUBE(x)   R_pow_di(x, 3)
#define MAX(a,b)  (((a)>(b)) ? (a) : (b))
#define MIN(a,b)  (((a)<(b)) ? (a) : (b))
#define SQR(x)    R_pow_di(x, 2)
#define repeat    for(;;)

/* dims structure */
typedef struct DIMS_struct {
  int
    n,      /* number of observations */
    p;      /* number of variables */
} DIMS_struct, *DIMS;

/* dims functions */
extern DIMS dims(int *);
extern void dims_free(DIMS);

#endif /* BASE_H */
