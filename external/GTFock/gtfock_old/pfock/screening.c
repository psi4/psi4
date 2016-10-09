#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <ga.h>
#include <mpi.h>

#include "CInt.h"
#include "config.h"
#include "screening.h"
#include "taskq.h"

#include "../../CifiedFxns.h"

int schwartz_screening(PFock_t pfock, BasisSet_t basis)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // create shell pairs values
    //ERD_t erd;
    int nthreads = omp_get_max_threads();
    //CInt_createERD(basis, &erd, nthreads);
    int nshells = pfock->nshells;

    // create global arrays for screening
    int nprow = pfock->nprow;
    int npcol = pfock->npcol;
    int dims[2];
    int block[2];
    int map[nprow + npcol];
    for (int i = 0; i < nprow; i++) {
        map[i] = pfock->rowptr_sh[i];
    }
    for (int i = 0; i < npcol; i++) {
        map[i + nprow] = pfock->colptr_sh[i];
    }
    dims[0] = nshells;
    dims[1] = nshells;
    block[0] = nprow;
    block[1] = npcol;
    pfock->ga_screening =
        NGA_Create_irreg(C_DBL, 2, dims, "array Screening", block, map);
    if (0 == pfock->ga_screening) {
        return -1;
    }

    // compute the max shell value
    double *sq_values = (double *)PFOCK_MALLOC(sizeof(double) *
        pfock->nshells_row * pfock->nshells_col);
    if (NULL == sq_values) {
        return -1;
    }
    int startM = pfock->sshell_row;
    int startN = pfock->sshell_col;
    int endM = pfock->eshell_row;
    int endN = pfock->eshell_col;
    double maxtmp = 0.0;
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        #pragma omp for reduction(max:maxtmp)
        for (int M = startM; M <= endM; M++) {
            int dimM = CInt_getShellDim(basis, M);
            for (int N = startN; N <= endN; N++) {
                int dimN = CInt_getShellDim(basis, N);
                double *integrals;
                int nints=
                ComputeShellQuartet(basis,tid,M,N,M,N,&integrals);
                //CInt_computeShellQuartet(basis, erd, tid, M, N, M, N,
                //                         &integrals, &nints);
                double maxvalue = 0.0;
                if (nints != 0) {
                    for (int iM = 0; iM < dimM; iM++) {
                        for (int iN = 0; iN < dimN; iN++) {
                            int index =
                                iM * (dimN*dimM*dimN+dimN) + iN * (dimM*dimN+1);
                            if (maxvalue < fabs(integrals[index])) {
                                maxvalue = fabs(integrals[index]);
                            }
                        }
                    }
                }
                sq_values[(M - startM) * (endN - startN + 1)  + (N - startN)]
                    = maxvalue;
                if (maxvalue > maxtmp) {
                    maxtmp = maxvalue;
                }
            }
        }
    }
    int lo[2];
    int hi[2];
    lo[0] = startM;
    hi[0] = endM;
    lo[1] = startN;
    hi[1] = endN;
    int ld = endN - startN + 1;
    NGA_Put(pfock->ga_screening, lo, hi, sq_values, &ld);
    // max value
    MPI_Allreduce(&maxtmp, &(pfock->maxvalue), 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    //CInt_destroyERD(erd);
    PFOCK_FREE(sq_values);

    // init shellptr
    sq_values = (double *)PFOCK_MALLOC(sizeof(double) * nshells);
    if (NULL == sq_values) {
        return -1;
    }
    int nnz = 0;
    double eta = pfock->tolscr2 / pfock->maxvalue;
    pfock->shellptr = (int *)PFOCK_MALLOC(sizeof(int) * (nshells + 1));
    pfock->mem_cpu += 1.0 * sizeof(int) * (nshells + 1);
    if (NULL == pfock->shellptr) {
        return -1;
    }
    memset(pfock->shellptr, 0, sizeof(int) * (nshells + 1));
    for (int M = 0; M < nshells; M++) {
        pfock->shellptr[M] = nnz;
        lo[0] = M;
        hi[0] = M;
        lo[1] = 0;
        hi[1] = nshells - 1;
        ld = nshells;
        NGA_Get(pfock->ga_screening, lo, hi, sq_values, &ld);
        for (int N = 0; N < nshells; N++) {
            double maxvalue = sq_values[N];
            if (maxvalue > eta) {
                if (M > N && (M + N) % 2 == 1 || M < N && (M + N) % 2 == 0) {
                    continue;
                } else {
                    nnz++;
                }
            }
        }
        pfock->shellptr[M + 1] = nnz;
    }
    pfock->nnz = nnz;

    double maxvalue;
    pfock->shellvalue  = (double *)PFOCK_MALLOC(sizeof(double) * nnz);
    pfock->shellid  = (int *)PFOCK_MALLOC(sizeof(int) * nnz);
    pfock->shellrid  = (int *)PFOCK_MALLOC(sizeof(int) * nnz);
    pfock->mem_cpu += 1.0 * sizeof(double) * nnz + 2.0 * sizeof(int) * nnz;
    nshells = pfock->nshells;
    if (pfock->shellvalue == NULL ||
        pfock->shellid == NULL ||
        pfock->shellrid == NULL) {
        return -1;
    }
    nnz = 0;
    for (int A = 0; A < nshells; A++) {
        pfock->shellptr[A] = nnz;
        lo[0] = A;
        hi[0] = A;
        lo[1] = 0;
        hi[1] = nshells - 1;
        ld = nshells;
        NGA_Get(pfock->ga_screening, lo, hi, sq_values, &ld);
        for (int B = 0; B < nshells; B++) {
            maxvalue = sq_values[B];
            if (maxvalue > eta) {
                if (A > B && (A + B) % 2 == 1 || A < B && (A + B) % 2 == 0)
                    continue;
                if (A == B) {
                    pfock->shellvalue[nnz] = maxvalue;
                } else {
                    pfock->shellvalue[nnz] = -maxvalue;
                }
                pfock->shellid[nnz] = B;
                pfock->shellrid[nnz] = A;
                nnz++;
            }
        }
    }
    PFOCK_FREE(sq_values);
    GA_Destroy(pfock->ga_screening);

    return 0;
}


void clean_screening(PFock_t pfock)
{
    PFOCK_FREE(pfock->shellid);
    PFOCK_FREE(pfock->shellrid);
    PFOCK_FREE(pfock->shellptr);
    PFOCK_FREE(pfock->shellvalue);
}
