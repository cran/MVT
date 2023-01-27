/* ID: interface.h, last updated 2022-08-21, F.Osorio */

#ifndef MVT_INTERFACE_H
#define MVT_INTERFACE_H

/* basic matrix manipulations */
extern void ax_plus_y(double, double *, int, double *, int, int);
extern void copy_lower(double *, int, double *, int, int);
extern void copy_mat(double *, int, double *, int, int, int);
extern void copy_vec(double *, int, double *, int, int);
extern double dot_product(double *, int, double *, int, int);
extern void GAXPY(double *, double, double *, int, int, int, double *, double);
extern double logAbsDet(double *, int, int);
extern void mult_triangular_vec(double *, int, int, char *, char *, char *, double *, int);
extern void mult_triangular_mat(double, double *, int, int, int, char *, char *, char *, char *, double *, int);
extern void scale_mat(double *, int, double, double *, int, int, int);
extern void scale_vec(double *, int, int, double);
extern void setzero(double *, int, int, int);
extern double sum_lower_tri(double *, int, int, int);
extern double trace_mat(double *, int, int);

/* routines for matrix decompositions */
extern void chol_decomp(double *, int, int, int, int *);

/* triangular solver */
extern void backsolve(double *, int, int, double *, int, int, int *);
extern void invert_mat(double *, int, int, int *);

/* descriptive statistics */
extern void center_and_Scatter(double *, int, int, double *, double *, double *);
extern void center_online(double *, int, int, double *, double *);

/* Mahalanobis distance */
extern double mahalanobis(double *, int, double *, double *);

/* Wilson-Hilferty transformation */
extern void WH_chisq(double *, int, int, double *);
extern void WH_F(double *, int, int, double, double *);

/* Utilities */
extern void cov2cor(double *, int);

#endif /* MVT_INTERFACE_H */
