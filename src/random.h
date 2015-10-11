#ifndef RANDOM_H
#define RANDOM_H

#include "matrix.h"

/* multivariate Student-t random generation (to be called by R) */
extern void rand_student(double *, int *, double *, double *, double *);
extern void rand_norm(double *, int *, double *, double *);

#endif /* RANDOM_H */
