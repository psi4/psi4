#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "oed_integral.h"
#include "basisset.h"
#include "cint_def.h"
#include "config.h"


static void config_oed (OED_t  oed, int A, int B, BasisSet_t basis)
{
    int cc_offset_A;
    int cc_offset_B;
    int alpha_offset_A;
    int alpha_offset_B;

    alpha_offset_A = 0;
    cc_offset_A = 0;
    alpha_offset_B = basis->nexp[A];
    cc_offset_B = basis->nexp[A];
    oed->ncoeff = alpha_offset_B + basis->nexp[B];
    oed->nalpha = cc_offset_B + basis->nexp[B];
    
    memcpy (&(oed->alpha[alpha_offset_A]), basis->exp[A], sizeof(double) * basis->nexp[A]);
    memcpy (&(oed->alpha[alpha_offset_B]), basis->exp[B], sizeof(double) * basis->nexp[B]);
    memcpy (&(oed->cc[cc_offset_A]), basis->cc[A], sizeof(double) * basis->nexp[A]);
    memcpy (&(oed->cc[cc_offset_B]), basis->cc[B], sizeof(double) * basis->nexp[B]);

    oed->npgto1 = basis->nexp[A];
    oed->npgto2 = basis->nexp[B];
    oed->cc_end[0] = oed->npgto1;
    oed->cc_end[1] = oed->npgto2;    
    oed->x1 = basis->x[A];
    oed->y1 = basis->y[A];
    oed->z1 = basis->z[A];
    oed->x2 = basis->x[B];
    oed->y2 = basis->y[B];
    oed->z2 = basis->z[B];
    oed->shell1 = basis->momentum[A];
    oed->shell2 = basis->momentum[B];
}


static void oed_max_scratch (BasisSet_t basis, OED_t oed)
{
    int max_momentum;
    int max_primid;
    int int_memory_min = 0;
    int int_memory_opt = 0;
    int fp_memory_min = 0;
    int fp_memory_opt = 0;
    
    _maxMomentum (basis, &max_momentum);
    _maxPrimid (basis, &max_primid);
    
    config_oed (oed,
                max_primid, max_primid,
                basis); 

    oed->shell1 = max_momentum;
    oed->shell2 = max_momentum;
    oed->x1 = 1.0;
    oed->x2 = 2.0;
    oed->y1 = 1.0;
    oed->y2 = 2.0;
    oed->z1 = 1.0;
    oed->z2 = 2.0;

    oed__memory_kin_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            oed->alpha, oed->cc, &(oed->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);
    oed->int_memory_opt = oed->int_memory_opt > int_memory_opt ?
            oed->int_memory_opt : int_memory_opt;
    oed->fp_memory_opt = oed->fp_memory_opt > fp_memory_opt ?
            oed->fp_memory_opt : fp_memory_opt;

    oed__memory_ovl_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            oed->alpha, oed->cc, &(oed->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);
    oed->int_memory_opt = oed->int_memory_opt > int_memory_opt ?
            oed->int_memory_opt : int_memory_opt;
    oed->fp_memory_opt = oed->fp_memory_opt > fp_memory_opt ?
            oed->fp_memory_opt : fp_memory_opt;
    
    oed__memory_nai_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            &(oed->natoms), oed->alpha, oed->cc, &(oed->spheric),                          
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);     
    oed->int_memory_opt = oed->int_memory_opt > int_memory_opt ?
            oed->int_memory_opt : int_memory_opt;
    oed->fp_memory_opt = oed->fp_memory_opt > fp_memory_opt ?
            oed->fp_memory_opt : fp_memory_opt;
}


CIntStatus_t CInt_createOED (BasisSet_t basis, OED_t *oed)
{
    OED_t o;
    int max_nexp;
    
    o = (OED_t )calloc (1, sizeof(struct OED));
    if (NULL == o)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;
    }
    o->ncsum = 2;
    o->ncgto1 = 1;
    o->ncgto2 = 1;
    o->cc_beg[0] = 1;
    o->cc_beg[1] = 1;
    o->natoms = basis->natoms;
    o->xn = basis->xn;
    o->yn = basis->yn;
    o->zn = basis->zn;
    o->ncharge = basis->ncharge;
    o->spheric = OED_SPHERIC;
    o->screen = OED_SCREEN;
    
    _maxnumExp (basis, &max_nexp);
    o->cc = (double *)malloc (2 * max_nexp * sizeof(double));
    o->alpha = (double *)malloc (2 * max_nexp * sizeof(double));
    if (NULL == o->cc || NULL == o->alpha)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;
    }
    o->coef_offset = (int *)malloc (basis->nshells * sizeof(int));
    o->exp_offset = (int *)malloc (basis->nshells * sizeof(int));
    if (NULL == o->coef_offset || NULL == o->exp_offset)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;
    }
        
    oed_max_scratch (basis, o);
    o->zcore = (double *)malloc (o->fp_memory_opt * sizeof(double));
    o->zcore2 = (double *)malloc (o->fp_memory_opt * sizeof(double));
    o->icore = (int *)malloc (o->int_memory_opt * sizeof(int));
    if (NULL == o->zcore ||
        NULL == o->zcore2 ||
        NULL == o->icore)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;
    }
    
    o->zmax = o->fp_memory_opt;
    o->imax = o->int_memory_opt;
    
    *oed = o;

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_destroyOED (OED_t  oed)
{
    free (oed->zcore);
    free (oed->icore);
    free (oed->alpha);
    free (oed->cc);
    free (oed->coef_offset);
    free (oed->exp_offset);
    free (oed);

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_computePairKin (BasisSet_t basis, OED_t oed,
                                  int A, int B,
                                  double **integrals, int *nints)
{
    int nfirst;
    
    if (A < 0 || A >= basis->nshells ||
        B < 0 || B >= basis->nshells)
    {
        CINT_PRINTF (1, "invalid shell indices\n");
        return CINT_STATUS_INVALID_VALUE;
    }

    config_oed (oed, A, B, basis);

#if ( _DEBUG_LEVEL_ == 3 )
    int int_memory_min;
    int int_memory_opt;
    int fp_memory_min;
    int fp_memory_opt;
    oed__memory_kin_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            oed->alpha, oed->cc, &(oed->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);
    
    assert (fp_memory_opt <= oed->fp_memory_opt);
    assert (int_memory_opt <= oed->int_memory_opt);   
#endif

    oed__gener_kin_batch_ (&(oed->imax), &(oed->zmax),
                           &(oed->nalpha), &(oed->ncoeff), &(oed->ncsum),
                           &(oed->ncgto1), &(oed->ncgto2),
                           &(oed->npgto1), &(oed->npgto2),
                           &(oed->shell1), &(oed->shell2),
                           &(oed->x1), &(oed->y1), &(oed->z1),
                           &(oed->x2), &(oed->y2), &(oed->z2),
                           oed->alpha, oed->cc,
                           oed->cc_beg, oed->cc_end, &(oed->spheric), &(oed->screen),
                           oed->icore, nints, &nfirst, oed->zcore);

    *integrals = &(oed->zcore[nfirst - 1]);
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_computePairOvl (BasisSet_t basis, OED_t oed,
                                  int A, int B,
                                  double **integrals, int *nints)
{
    int nfirst;
    
    if (A < 0 || A >= basis->nshells ||
        B < 0 || B >= basis->nshells)
    {
        CINT_PRINTF (1, "invalid shell indices\n");
        return CINT_STATUS_INVALID_VALUE;
    }
    
    config_oed (oed, A, B, basis);

#if ( _DEBUG_LEVEL_ == 3 )
    int int_memory_min;
    int int_memory_opt;
    int fp_memory_min;
    int fp_memory_opt;
    oed__memory_ovl_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            oed->alpha, oed->cc, &(oed->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);
    
    assert (fp_memory_opt <= oed->fp_memory_opt);
    assert (int_memory_opt <= oed->int_memory_opt);   
#endif

    oed__gener_ovl_batch_ (&(oed->imax), &(oed->zmax),
                           &(oed->nalpha), &(oed->ncoeff), &(oed->ncsum),
                           &(oed->ncgto1), &(oed->ncgto2),
                           &(oed->npgto1), &(oed->npgto2),
                           &(oed->shell1), &(oed->shell2),
                           &(oed->x1), &(oed->y1), &(oed->z1),
                           &(oed->x2), &(oed->y2), &(oed->z2),
                           oed->alpha, oed->cc, oed->cc_beg, oed->cc_end,
                           &(oed->spheric), &(oed->screen),
                           oed->icore, nints, &nfirst, oed->zcore);

    *integrals = &(oed->zcore[nfirst - 1]);  
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_computePairPot (BasisSet_t basis, OED_t oed,
                                  int A, int B,
                                  double **integrals, int *nints)
{
    int nfirst;
    if (A < 0 || A >= basis->nshells ||
        B < 0 || B >= basis->nshells)
    {
        CINT_PRINTF (1, "invalid shell indices\n");
        return CINT_STATUS_INVALID_VALUE;
    }
    
    config_oed (oed, A, B, basis);

#if ( _DEBUG_LEVEL_ == 3 )
    int int_memory_min;
    int int_memory_opt;
    int fp_memory_min;
    int fp_memory_opt;
    oed__memory_nai_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            &(oed->natoms), oed->alpha, oed->cc, &(oed->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);
    
    assert (fp_memory_opt <= oed->fp_memory_opt);
    assert (int_memory_opt <= oed->int_memory_opt);   
#endif

    oed__gener_nai_batch_ (&(oed->imax), &(oed->zmax),
                           &(oed->nalpha), &(oed->ncoeff), &(oed->ncsum),
                           &(oed->ncgto1), &(oed->ncgto2),
                           &(oed->npgto1), &(oed->npgto2),
                           &(oed->shell1), &(oed->shell2),
                           &(oed->x1), &(oed->y1), &(oed->z1),
                           &(oed->x2), &(oed->y2), &(oed->z2),
                           &(oed->natoms),
                           oed->xn, oed->yn, oed->zn,
                           oed->ncharge, oed->alpha, oed->cc,
                           oed->cc_beg, oed->cc_end,
                           &(oed->spheric), &(oed->screen),
                           oed->icore, nints, &nfirst, oed->zcore);

    *integrals = &(oed->zcore[nfirst - 1]);
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_computePairCoreH (BasisSet_t basis, OED_t oed,
                                    int A, int B,
                                    double **integrals, int *nints)
{
    int nfirst;
    int nfirst2;
    int ni;
    int ni2;
    int i;

    if (A < 0 || A >= basis->nshells ||
        B < 0 || B >= basis->nshells)
    {
        CINT_PRINTF (1, "invalid shell indices\n");
        return CINT_STATUS_INVALID_VALUE;
    }
    
    config_oed (oed, A, B, basis);

#if ( _DEBUG_LEVEL_ == 3 )
    int int_memory_min;
    int int_memory_opt;
    int fp_memory_min;
    int fp_memory_opt;
    oed__memory_kin_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            oed->alpha, oed->cc, &(oed->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);

    assert (fp_memory_opt <= oed->fp_memory_opt);
    assert (int_memory_opt <= oed->int_memory_opt);
    oed__memory_nai_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            &(oed->natoms), oed->alpha, oed->cc, &(oed->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);

    assert (fp_memory_opt <= oed->fp_memory_opt);
    assert (int_memory_opt <= oed->int_memory_opt);   
#endif

    oed__gener_kin_batch_ (&(oed->imax), &(oed->zmax),
                           &(oed->nalpha), &(oed->ncoeff), &(oed->ncsum),
                           &(oed->ncgto1), &(oed->ncgto2),
                           &(oed->npgto1), &(oed->npgto2),
                           &(oed->shell1), &(oed->shell2),
                           &(oed->x1), &(oed->y1), &(oed->z1),
                           &(oed->x2), &(oed->y2), &(oed->z2),
                           oed->alpha, oed->cc, oed->cc_beg, oed->cc_end,
                           &(oed->spheric), &(oed->screen),
                           oed->icore, &ni2, &nfirst2, oed->zcore2);
    
    oed__gener_nai_batch_ (&(oed->imax), &(oed->zmax),
                           &(oed->nalpha), &(oed->ncoeff), &(oed->ncsum),
                           &(oed->ncgto1), &(oed->ncgto2),
                           &(oed->npgto1), &(oed->npgto2),
                           &(oed->shell1), &(oed->shell2),
                           &(oed->x1), &(oed->y1), &(oed->z1),
                           &(oed->x2), &(oed->y2), &(oed->z2),
                           &(oed->natoms),
                           oed->xn, oed->yn, oed->zn,
                           oed->ncharge, oed->alpha, oed->cc,
                           oed->cc_beg, oed->cc_end,
                           &(oed->spheric), &(oed->screen),
                           oed->icore, &ni, &nfirst, oed->zcore);

    *integrals = &(oed->zcore[nfirst - 1]);
    if (ni != 0 && ni2 != 0)
    {
        for (i = 0; i < ni; i++)
        {
            oed->zcore[nfirst - 1 + i] += oed->zcore2[nfirst2 - 1 + i];           
        }
        *nints = ni;
    }
    else if (ni != 0)
    {
        *nints = ni;
    }
    else if (ni2 != 0)
    {
        *integrals = &(oed->zcore2[nfirst2 - 1]);
        *nints = ni2;
    }
    else
    {
        *nints = 0;
    }
    
    return CINT_STATUS_SUCCESS;
}
