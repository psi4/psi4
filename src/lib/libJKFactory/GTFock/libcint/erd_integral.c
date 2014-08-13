#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "erd_integral.h"
#include "basisset.h"
#include "config.h"
#include "cint_def.h"

#include <sys/param.h>
#include <sys/times.h>


static void config_erd (ERD_t erd, int A, int B, int C, int D, BasisSet_t basis)
{
    int cc_offset_A;
    int cc_offset_B;
    int cc_offset_C;
    int cc_offset_D;
    int alpha_offset_A;
    int alpha_offset_B;
    int alpha_offset_C;
    int alpha_offset_D;
    
    alpha_offset_A = 0;
    cc_offset_A = 0;
    alpha_offset_B = basis->nexp[A];
    cc_offset_B = basis->nexp[A];
    alpha_offset_C = alpha_offset_B + basis->nexp[B];
    cc_offset_C = cc_offset_B + basis->nexp[B];
    alpha_offset_D = alpha_offset_C + basis->nexp[C];
    cc_offset_D = cc_offset_C + basis->nexp[C];
    erd->nalpha = alpha_offset_D + basis->nexp[D]; 
    erd->ncoef =  cc_offset_D + basis->nexp[D];
    erd->spheric = basis->basistype;
         
    memcpy (&(erd->alpha[alpha_offset_A]), basis->exp[A], sizeof(double) * basis->nexp[A]);
    memcpy (&(erd->alpha[alpha_offset_B]), basis->exp[B], sizeof(double) * basis->nexp[B]);
    memcpy (&(erd->alpha[alpha_offset_C]), basis->exp[C], sizeof(double) * basis->nexp[C]);
    memcpy (&(erd->alpha[alpha_offset_D]), basis->exp[D], sizeof(double) * basis->nexp[D]); 
    memcpy (&(erd->cc[cc_offset_A]), basis->cc[A], sizeof(double) * basis->nexp[A]);
    memcpy (&(erd->cc[cc_offset_B]), basis->cc[B], sizeof(double) * basis->nexp[B]);
    memcpy (&(erd->cc[cc_offset_C]), basis->cc[C], sizeof(double) * basis->nexp[C]);
    memcpy (&(erd->cc[cc_offset_D]), basis->cc[D], sizeof(double) * basis->nexp[D]);

    erd->npgto1 = basis->nexp[A];
    erd->npgto2 = basis->nexp[B];
    erd->npgto3 = basis->nexp[C];
    erd->npgto4 = basis->nexp[D];
    erd->cc_end[0] = erd->npgto1;
    erd->cc_end[1] = erd->npgto2;
    erd->cc_end[2] = erd->npgto3;
    erd->cc_end[3] = erd->npgto4;
    erd->x1 = basis->x[A];
    erd->y1 = basis->y[A];
    erd->z1 = basis->z[A];
    erd->x2 = basis->x[B];
    erd->y2 = basis->y[B];
    erd->z2 = basis->z[B];
    erd->x3 = basis->x[C];
    erd->y3 = basis->y[C];
    erd->z3 = basis->z[C];
    erd->x4 = basis->x[D];
    erd->y4 = basis->y[D];
    erd->z4 = basis->z[D];
    erd->shell1 = basis->momentum[A];
    erd->shell2 = basis->momentum[B];
    erd->shell3 = basis->momentum[C];
    erd->shell4 = basis->momentum[D];
}


static void erd_max_scratch (BasisSet_t basis, ERD_t erd)
{
    int max_momentum;
    int max_primid;
    int int_memory_min = 0;
    int int_memory_opt = 0;
    int fp_memory_min = 0;
    int fp_memory_opt = 0;
    
    _maxMomentum (basis, &max_momentum);
    _maxPrimid (basis, &max_primid);
    
    config_erd (erd,
                max_primid, max_primid,
                max_primid, max_primid,
                basis);
    
    erd->shell1 = max_momentum;
    erd->shell2 = max_momentum;
    erd->shell3 = max_momentum;
    erd->shell4 = max_momentum;
    erd->x1 = 1.0;
    erd->x2 = 2.0;
    erd->x3 = 3.0;
    erd->x4 = 4.0;
    erd->y1 = 1.0;
    erd->y2 = 2.0;
    erd->y3 = 3.0;
    erd->y4 = 4.0;
    erd->z1 = 1.0;
    erd->z2 = 2.0;
    erd->z3 = 3.0;
    erd->z4 = 4.0;
  
    erd__memory_eri_batch_ (&(erd->nalpha), &(erd->ncoef),
                            &(erd->ncgto1), &(erd->ncgto2),
                            &(erd->ncgto3), &(erd->ncgto4),
                            &(erd->npgto1), &(erd->npgto2),
                            &(erd->npgto3), &(erd->npgto4),
                            &(erd->shell1), &(erd->shell2),
                            &(erd->shell3), &(erd->shell4),
                            &(erd->x1), &(erd->y1), &(erd->z1),
                            &(erd->x2), &(erd->y2), &(erd->z2),
                            &(erd->x3), &(erd->y3), &(erd->z3),
                            &(erd->x4), &(erd->y4), &(erd->z4),
                            erd->alpha, erd->cc, &(erd->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);
    erd->int_memory_opt = erd->int_memory_opt > int_memory_opt ?
            erd->int_memory_opt : int_memory_opt;
    erd->fp_memory_opt = erd->fp_memory_opt > fp_memory_opt ?
            erd->fp_memory_opt : fp_memory_opt;
}


CIntStatus_t CInt_createERD (BasisSet_t basis, ERD_t *erd)
{
    ERD_t e;
    int max_ncoef;
    int max_nexp;

    e = (ERD_t)calloc (1, sizeof(struct ERD));
    if (NULL == e)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;
    }

    e->ncgto1 = 1;
    e->ncgto2 = 1;
    e->ncgto3 = 1;
    e->ncgto4 = 1;    
    e->cc_beg[0] = 1;
    e->cc_beg[1] = 1;
    e->cc_beg[2] = 1;
    e->cc_beg[3] = 1;
    e->ncsum = 4;
    e->spheric = ERD_SPHERIC;
    e->screen = ERD_SCREEN;

    _maxnumExp (basis, &max_nexp);
    max_ncoef = max_nexp;
    e->cc = (double *)malloc (4 * max_ncoef * sizeof(double));
    e->alpha = (double *)malloc (4 * max_nexp * sizeof(double));
    if (NULL == e->cc || NULL == e->alpha)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;
    }

    erd_max_scratch (basis, e);
    e->zcore = (double *)malloc (e->fp_memory_opt * sizeof(double));
    e->icore = (int *)malloc (e->int_memory_opt * sizeof(int));   
    if (NULL == e->zcore || NULL == e->icore)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;
    }
    e->zmax = e->fp_memory_opt;
    e->imax = e->int_memory_opt;
    
    *erd = e;

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_destroyERD (ERD_t erd)
{
    free (erd->zcore);
    free (erd->icore);
    free (erd->alpha);
    free (erd->cc);
    free (erd->coef_offset);
    free (erd->exp_offset);
    free (erd);
    
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_computeShellQuartet ( BasisSet_t basis, ERD_t erd,
                                        int A, int B, int C, int D,
                                        double **integrals, int *nints)
{
    int nfirst;

#if ( _DEBUG_LEVEL_ == 3 )
    if (A < 0 || A >= basis->nshells ||
        B < 0 || B >= basis->nshells ||
        C < 0 || C >= basis->nshells ||
        D < 0 || D >= basis->nshells)
    {
        CINT_PRINTF (1, "invalid shell indices\n");
        *nints = 0;
        return CINT_STATUS_INVALID_VALUE;
    }
#endif
    
    config_erd (erd, A, B, C, D, basis);

#if ( _DEBUG_LEVEL_ == 3 )
    int int_memory_min;
    int int_memory_opt;
    int fp_memory_min;
    int fp_memory_opt;
    erd__memory_eri_batch_ (&(erd->nalpha), &(erd->ncoef),
                            &(erd->ncgto1), &(erd->ncgto2),
                            &(erd->ncgto3), &(erd->ncgto4),
                            &(erd->npgto1), &(erd->npgto2),
                            &(erd->npgto3), &(erd->npgto4),
                            &(erd->shell1), &(erd->shell2),
                            &(erd->shell3), &(erd->shell4),
                            &(erd->x1), &(erd->y1), &(erd->z1),
                            &(erd->x2), &(erd->y2), &(erd->z2),
                            &(erd->x3), &(erd->y3), &(erd->z3),
                            &(erd->x4), &(erd->y4), &(erd->z4),
                            erd->alpha, erd->cc, &(erd->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);
    assert (fp_memory_opt <= erd->fp_memory_opt);
    assert (int_memory_opt <= erd->int_memory_opt);   
#endif

    erd__gener_eri_batch_ (&(erd->imax), &(erd->zmax),
                           &(erd->nalpha), &(erd->ncoef),
                           &(erd->ncsum), &(erd->ncgto1),
                           &(erd->ncgto2), &(erd->ncgto3),
                           &(erd->ncgto4), &(erd->npgto1),
                           &(erd->npgto2), &(erd->npgto3),
                           &(erd->npgto4), &(erd->shell1),
                           &(erd->shell2), &(erd->shell3),
                           &(erd->shell4),
                           &(erd->x1), &(erd->y1), &(erd->z1),
                           &(erd->x2), &(erd->y2), &(erd->z2),
                           &(erd->x3), &(erd->y3), &(erd->z3),
                           &(erd->x4), &(erd->y4), &(erd->z4),
                           erd->alpha, erd->cc,
                           erd->cc_beg, erd->cc_end,
                           &(erd->spheric), &(erd->screen),
                           erd->icore, nints,
                           &nfirst, erd->zcore);

    *integrals = &(erd->zcore[nfirst - 1]);

    return CINT_STATUS_SUCCESS;
}
