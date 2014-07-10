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


int init_screening (PFock_t pfock, BasisSet_t basis)
{
    ERD_t erd;
    int M;
    int N;
    int dimM;
    int dimN;
    int iM;
    int iN;
    int index;
    double maxvalue;
    double maxtmp;    
    int startM;
    int endM;
    int startN;
    int endN;
    int nshells;
    int myrank;
    double *sq_values;
    double *integrals;
    int nints;
    int lo[2];
    int hi[2];
    int ld;
    int nnz;
    double eta;
    int nprow;
    int npcol;
    int *map;
    int i;
    int dims[2];
    int block[2];

    MPI_Comm_rank (MPI_COMM_WORLD, &myrank); 

    // create shell pairs values
    CInt_createERD (basis, &erd);
    maxtmp = 0.0;   
    nshells = pfock->nshells;
    
    // create global arrays for screening 
    nprow = pfock->nprow;
    npcol = pfock->npcol;
    map = (int *)PFOCK_MALLOC (sizeof(int) * (nprow + npcol));
    if (NULL == map)
    {
        return -1;
    }

    for (i = 0; i < nprow; i++)
    {
        map[i] = pfock->rowptr_sh[i];
    }   
    for (i = 0; i < npcol; i++)
    {
        map[i + nprow] = pfock->colptr_sh[i];
    }
    dims[0] = nshells;
    dims[1] = nshells;
    block[0] = nprow;
    block[1] = npcol;             
    pfock->ga_screening = NGA_Create_irreg (C_DBL, 2, dims, "array Screening", block, map);
    if (0 == pfock->ga_screening)
    {
        return -1;
    }
    
    sq_values = (double *)PFOCK_MALLOC (sizeof(double) * pfock->nshells_row * pfock->nshells_col);
    if (NULL == sq_values)
    {
        return -1;
    }
    startM = pfock->sshell_row;
    startN = pfock->sshell_col;
    endM = pfock->eshell_row;
    endN = pfock->eshell_col;   
    for (M = startM; M <= endM; M++)
    {
        dimM = CInt_getShellDim (basis, M);
        for (N = startN; N <= endN; N++)
        {
            dimN = CInt_getShellDim (basis, N);
            CInt_computeShellQuartet (basis, erd, M, N, M, N, &integrals, &nints);            
            maxvalue = 0.0;
            if (nints != 0)
            {
                for (iM = 0; iM < dimM; iM++)
                {
                    for (iN = 0; iN < dimN; iN++)
                    {
                        index = iM * (dimN*dimM*dimN+dimN) + iN * (dimM*dimN+1);
                        if (maxvalue < fabs (integrals[index]))
                        {
                            maxvalue = fabs (integrals[index]);                    
                        }
                    }
                }
            }
            sq_values[(M - startM) * (endN - startN + 1)  + (N - startN)] = maxvalue;
            if (maxvalue > maxtmp)
            {
                maxtmp = maxvalue;
            }
        }
        lo[0] = startM;
        hi[0] = endM;
        lo[1] = startN;
        hi[1] = endN;
        ld = endN - startN + 1;
        NGA_Put (pfock->ga_screening, lo, hi, sq_values, &ld);
    }
    // max value  
    MPI_Allreduce (&maxtmp, &maxvalue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    pfock->maxvalue = maxvalue;
    CInt_destroyERD (erd);
    free (sq_values);

    // init shellptr
    sq_values = (double *)PFOCK_MALLOC (sizeof(double) * nshells);
    if (NULL == sq_values)
    {
        return -1;
    }
    nnz = 0;
    eta = pfock->tolscr2 / maxvalue;
    pfock->shellptr = (int *)PFOCK_MALLOC (sizeof(int) * (nshells + 1));
    if (NULL == pfock->shellptr)
    {
        return -1;
    }
    memset (pfock->shellptr, 0, sizeof(int) * (nshells + 1));
    for (M = 0; M < nshells; M++)
    {
        pfock->shellptr[M] = nnz;
        lo[0] = M;
        hi[0] = M;
        lo[1] = 0;
        hi[1] = nshells - 1;
        ld = nshells;
        NGA_Get (pfock->ga_screening, lo, hi, sq_values, &ld);
        for (N = 0; N < nshells; N++)
        {
            maxvalue = sq_values[N];
            if (maxvalue > eta)
            {
                if (M > N && (M + N) % 2 == 1 || M < N && (M + N) % 2 == 0)
                {
                    continue;
                }
                else
                {
                    nnz++;
                }
            }
        }
        pfock->shellptr[M + 1] = nnz;
    }
    pfock->nnz = nnz;
    
    free (sq_values);

    return 0;
}


int schwartz_screening (PFock_t pfock, int startshell, int endshell)
{
    int A;
    int B;
    double maxvalue;    
    int nshells;   
    double eta;
    long long nnz;
    int myrank;
    double *sq_values;
    int lo[2];
    int hi[2];
    int ld;
    int size;
    
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

    maxvalue = pfock->maxvalue;    
    eta = pfock->tolscr2 / maxvalue;
    
    size = pfock->shellptr[endshell + 1] - pfock->shellptr[startshell];
    pfock->shellvalue  = (double *)PFOCK_MALLOC (sizeof(double) * size);
    pfock->shellid  = (int *)PFOCK_MALLOC (sizeof(int) * size);
    pfock->shellrid  = (int *)PFOCK_MALLOC (sizeof(int) * size);
    nshells = pfock->nshells;
    sq_values = (double *)PFOCK_MALLOC (sizeof(double) * nshells);
    if (pfock->shellvalue == NULL ||
        pfock->shellid == NULL ||
        pfock->shellrid == NULL ||
        sq_values == NULL)
    {
        return -1;
    }
    
    nnz = 0;
    for (A = startshell; A <= endshell; A++)
    {
        pfock->shellptr[A] = nnz;
        lo[0] = A;
        hi[0] = A;
        lo[1] = 0;
        hi[1] = nshells - 1;
        ld = nshells;
        NGA_Get (pfock->ga_screening, lo, hi, sq_values, &ld);
        for (B = 0; B < nshells; B++)
        {
            maxvalue = sq_values[B];
            if (maxvalue > eta)
            {
                if (A > B && (A + B) % 2 == 1 || A < B && (A + B) % 2 == 0)
                    continue;
                if (A == B)
                {
                    pfock->shellvalue[nnz] = maxvalue;                       
                }
                else
                {
                    pfock->shellvalue[nnz] = -maxvalue;
                }
                pfock->shellid[nnz] = B;
                pfock->shellrid[nnz] = A;
                nnz++;
            }
        }
    }
    free (sq_values);

    return 0;
}


void clean_screening (PFock_t pfock)
{
    free (pfock->shellid);
    free (pfock->shellrid);
    free (pfock->shellptr);
    free (pfock->shellvalue);
    GA_Destroy (pfock->ga_screening);
}
