#ifndef OPTIM_H
#define OPTIM_H

#include "base.h"

/* Brent's method for unidimensional optimization */
extern double brent(double, double, double (*f)(double, void *), void *, double);

#endif /* OPTIM_H */
