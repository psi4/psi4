#ifndef __ERD_INTEGRAL_H__
#define __ERD_INTEGRAL_H__


#define ERD_SCREEN  0
#define ERD_SPHERIC 1


// Main fortran function of the ERD package.

	void erd__gener_eri_batch_ (int *imax, int *zmax,
                                   int *nalpha, int *ncoeff, int *ncsum,
                                   int *ncgto1, int *ncgto2, int *ncgto3,
                                   int *ncgto4, int *npgto1, int *npgto2,
                                   int *npgto3, int *npgto4, int *shell1,
                                   int *shell2, int *shell3, int *shell4,
                                   double *x1, double *y1, double *z1,
                                   double *x2, double *y2, double *z2,
                                   double *x3, double *y3, double *z3,
                                   double *x4, double *y4, double *z4,
                                   double *alpha, double *cc, int *ccbeg,
                                   int *ccend, int *spheric, int *screen,
                                   int *icore, int *nbatch, int *nfirst,
                                   double *zcore);
void erd__memory_eri_batch_ (int *nalpha, int *ncoeff,
                                    int *ncgto1, int *ncgto2, int *ncgto3,
                                    int *ncgto4, int *npgto1, int *npgto2,
                                    int *npgto3, int *npgto4, int *shell1,
                                    int *shell2, int *shell3, int *shell4,
                                    double *x1, double *y1, double *z1,
                                    double *x2, double *y2, double *z2,
                                    double *x3, double *y3, double *z3,
                                    double *x4, double *y4, double *z4,
                                    double *alpha, double *cc, int *spheric,
                                    int *imin, int *iopt, int *zmin,
                                    int *zopt);

#endif /* __ERD_INTEGRAL_H__ */
