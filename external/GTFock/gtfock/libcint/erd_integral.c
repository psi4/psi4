#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "erd_integral.h"
#include "basisset.h"
#include "config.h"
#include "cint_def.h"

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif


static void erd_max_scratch(BasisSet_t basis, ERD_t erd) {
    const int max_momentum = basis->max_momentum;
    const int max_primid = basis->max_nexp_id;
    const int maxnpgto = basis->nexp[max_primid];
        
    if (max_momentum < 2) {
        erd->capacity = 81;
    } else {
        erd->capacity = erd__memory_csgto(
            maxnpgto, maxnpgto, maxnpgto, maxnpgto,
            max_momentum, max_momentum,
            max_momentum, max_momentum,
            1.0, 1.0, 1.0, 2.0, 2.0, 2.0,
            3.0, 3.0, 3.0, 4.0, 4.0, 4.0,
            ERD_SPHERIC);
    }
}


static CIntStatus_t create_vrrtable(BasisSet_t basis, ERD_t erd) {
    const int max_shella = basis->max_momentum + 1;
    erd->max_shella = max_shella;
    const int max_shellp = 2 * (max_shella - 1);
    const int tablesize = max_shellp + 1;
    const int total_combinations = (max_shellp + 1) * (max_shellp + 2) * (max_shellp + 3) / 6;
    int **vrrtable = (int **)malloc(sizeof(int *) * tablesize);
    CINT_ASSERT(vrrtable != NULL);

    int *vrrtable__ = (int *)malloc(sizeof(int) * 4 * total_combinations);
    CINT_ASSERT(vrrtable__ != NULL);

    int n = 0;
    for(int shella = 0; shella <= max_shellp; shella++)
    {
        vrrtable[shella] = &vrrtable__[n];
        int count = 0;
        for(int xf = shella; xf >= 0; xf--)
        {
            for(int yf = shella - xf; yf >= 0; yf--)
            {
                int zf = shella - xf - yf;
                    vrrtable[shella][count + 0] = xf;
                    vrrtable[shella][count + 1] = yf;
                    vrrtable[shella][count + 2] = zf;
                    vrrtable[shella][count + 3] = count;
                    count += 4;
            }
        }
        n += count;
    }
    erd->vrrtable = vrrtable;
    return CINT_STATUS_SUCCESS;
}

static CIntStatus_t destroy_vrrtable(ERD_t erd) {
    free(erd->vrrtable[0]);
    free(erd->vrrtable);
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_createERD(BasisSet_t basis, ERD_t *erd, int nthreads) {      
    CINT_ASSERT(nthreads > 0);

    // malloc erd
    ERD_t e = (ERD_t)calloc(1, sizeof(struct ERD));
    CINT_ASSERT(e != NULL);
    erd_max_scratch(basis, e);

    // memory scratch memory
    e->nthreads = nthreads;
    e->buffer = (double **)malloc(nthreads * sizeof(double *));
    CINT_ASSERT(e->buffer != NULL);
    for (int i = 0; i < nthreads; i++) {
        e->buffer[i] = (double *)ALIGNED_MALLOC(e->capacity * sizeof(double));
        CINT_ASSERT(e->buffer[i] != NULL);
    }

    // create vrr table
    const CIntStatus_t status = create_vrrtable(basis, e);
    CINT_ASSERT(status == CINT_STATUS_SUCCESS);
    CINT_INFO("totally use %.3lf MB (%.3lf MB per thread)",
        (e->fp_memory_opt * sizeof(double)
        + e->int_memory_opt * sizeof(int)) * nthreads/1024.0/1024.0,
        (e->fp_memory_opt * sizeof(double)
        + e->int_memory_opt * sizeof(int))/1024.0/1024.0);
    *erd = e;
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_destroyERD(ERD_t erd) {
    for (uint32_t i = 0; i < erd->nthreads; i++) {
        ALIGNED_FREE(erd->buffer[i]);
    }
    free(erd->buffer);

    destroy_vrrtable(erd);
    free(erd);

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_computeShellQuartet( BasisSet_t basis, ERD_t erd, int tid,
                                        int A, int B, int C, int D,
                                        double **integrals, int *nints)
{
#if ( _DEBUG_LEVEL_ == 3 )
    if (A < 0 || A >= basis->nshells ||
        B < 0 || B >= basis->nshells ||
        C < 0 || C >= basis->nshells ||
        D < 0 || D >= basis->nshells)
    {
#ifndef __INTEL_OFFLOAD
        CINT_PRINTF(1, "invalid shell indices\n");
#endif
        *nints = 0;
        return CINT_STATUS_INVALID_VALUE;
    }
    if (tid <= 0 ||
        tid >= erd->nthreads)
    {
#ifndef __INTEL_OFFLOAD
        CINT_PRINTF(1, "invalid thread id\n");
#endif
        *nints = 0;
        return CINT_STATUS_INVALID_VALUE;    
    }
#endif

    const uint32_t shell1 = basis->momentum[A];
    const uint32_t shell2 = basis->momentum[B];
    const uint32_t shell3 = basis->momentum[C];
    const uint32_t shell4 = basis->momentum[D];
    const uint32_t orshell = shell1 | shell2 | shell3 | shell4;
    if (orshell < 2) {
        uint32_t integrals_count = 0;
        erd__1111_csgto(
            A, B, C, D,
            basis->nexp, basis->momentum, basis->xyz0,
            (const double**)basis->exp, basis->minexp, (const double**)basis->cc, (const double**)basis->norm,
            erd->capacity, &integrals_count, erd->buffer[tid]);
        *nints = integrals_count;
    } else {
        uint32_t integrals_count = 0;
        erd__csgto(
            A, B, C, D,
            basis->nexp, basis->momentum, basis->xyz0,
            (const double**)basis->exp, basis->minexp, (const double**)basis->cc, (const double**)basis->norm,
            erd->vrrtable,
            basis->basistype,
            erd->capacity, &integrals_count, erd->buffer[tid]);
        *nints = integrals_count;
    }

    *integrals = erd->buffer[tid];

    if (*nints != 0 && isnan((*integrals)[0]))
    {
        printf ("NAN %d %d %d %d\n", A, B, C, D);
    }
    return CINT_STATUS_SUCCESS;
}

CIntStatus_t CInt_computeShellQuartets(BasisSet_t basis,
                                       ERD_t erd,
                                       uint32_t threadId,
                                       uint32_t shellIndicexA,
                                       const uint32_t*restrict shellIndicesB,
                                       uint32_t shellIndicexC,
                                       const uint32_t*restrict shellIndicesD,
                                       uint32_t shellIndicesCount,
                                       double **integrals,
                                       int *integralsCountPtr)
{
    int totalIntegralsCount = 0;
    for (uint32_t shellIndicesIndex = 0; shellIndicesIndex != shellIndicesCount; shellIndicesIndex += 1) {
        int integralsCount = 0;
        CIntStatus_t status = CInt_computeShellQuartet(basis, erd, threadId,
            shellIndicexA, shellIndicesB[shellIndicesIndex], shellIndicexC, shellIndicesD[shellIndicesIndex],
            integrals, &integralsCount);
        if (status != CINT_STATUS_SUCCESS) {
            *integralsCountPtr = totalIntegralsCount;
            return status;
        }
        totalIntegralsCount += integralsCount;
    }
    *integralsCountPtr = totalIntegralsCount;
    return CINT_STATUS_SUCCESS;
}


void CInt_getMaxMemory(ERD_t erd, double *memsize)
{
    *memsize = erd->capacity * sizeof(double) * erd->nthreads;
}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
