#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <ga.h>
#include <macdecls.h>
#include <time.h>

#include "config.h"
#include "taskq.h"
#include "fock_task.h"



static void update_F (double *integrals, int numdmat2,
                      int dimM, int dimN, int dimP, int dimQ,
                      int flag1, int flag2, int flag3, int flagsymm,
                      int iMN, int iPQ, int iMP, int iNP, int iMQ, int iNQ,
                      double **D1, double **D2, double **D3,
                      double **J1, double **J2, double **K3,
                      int ldX1, int ldX2, int ldX3)
{
    int iM;
    int iN;
    int iP;
    int iQ;
    double I;
    int flag4;
    int flag5;
    int flag6;
    int flag7;
    double *D_MN;
    double *D_PQ;
    double *D_MP;
    double *D_NP;
    double *D_MQ;
    double *D_NQ;    
    double *J_MN;
    double *J_PQ;
    double *K_MP;
    double *K_NP;
    double *K_MQ;
    double *K_NQ;
    int i;
    int index;
    int indexi;
    int indexj;
    int indexk;
    int indexl;

    flag4 = (flag1 == 1 && flag2 == 1) ? 1 : 0;
    flag5 = (flag1 == 1 && flag3 == 1) ? 1 : 0;
    flag6 = (flag2 == 1 && flag3 == 1) ? 1 : 0;
    flag7 = (flag4 == 1 && flag3 == 1) ? 1 : 0;                    
    for (i = 0 ; i < numdmat2; i++)
    {
        D_MN = D1[i] + iMN;               
        D_PQ = D2[i] + iPQ;                
        D_MP = D3[i] + iMP;
        D_MQ = D3[i] + iMQ;
        D_NP = D3[i] + iNP;
        D_NQ = D3[i] + iNQ;        
        J_MN = J1[i] + iMN;
        J_PQ = J2[i] + iPQ;
        K_MP = K3[i] + iMP;
        K_MQ = K3[i] + iMQ;
        K_NP = K3[i] + iNP;
        K_NQ = K3[i] + iNQ;

        for (iM = 0; iM < dimM; iM++)
        {
            for (iN = 0; iN < dimN; iN++)
            {
                for (iP = 0; iP < dimP; iP++)
                {
                    for (iQ = 0; iQ < dimQ; iQ++)
                    {
                        I = integrals[iM + dimM * (iN + dimN * (iP + dimP * iQ))];                    

                        // F(m, n) += D(p, q) * 2 * I(m, n, p, q)
                        // F(n, m) += D(p, q) * 2 * I(n, m, p, q)
                        // F(m, n) += D(q, p) * 2 * I(m, n, q, p)
                        // F(n, m) += D(q, p) * 2 * I(n, m, q, p)
                        J_MN[iM * ldX1 + iN] += 0.25 * (1 + flag1 + flag2 + flag4) *
                                D_PQ[iP * ldX2 + iQ] * I * (1 + flagsymm);

                        // F(p, q) += D(m, n) * 2 * I(p, q, m, n)
                        // F(p, q) += D(n, m) * 2 * I(p, q, n, m)
                        // F(q, p) += D(m, n) * 2 * I(q, p, m, n)
                        // F(q, p) += D(n, m) * 2 * I(q, p, n, m)
                        J_PQ[iP * ldX2 + iQ] += 0.25 * (flag3 + flag5 + flag6 + flag7) *
                                D_MN[iM * ldX1 + iN] * I * (1 + flagsymm);

                        // F(m, p) -= D(n, q) * I(m, n, p, q)
                        // F(p, m) -= D(q, n) * I(p, q, m, n)
                        K_MP[iM * ldX3 + iP] += (1 + flag3) *
                                0.5 * D_NQ[iN * ldX3 + iQ] * I;

                        // F(n, p) -= D(m, q) * I(n, m, p, q)
                        // F(p, n) -= D(q, m) * I(p, q, n, m)
                        K_NP[iN * ldX3 + iP] += (flag1 + flag5) *
                                0.5 * D_MQ[iM * ldX3 + iQ] * I;

                        // F(m, q) -= D(n, p) * I(m, n, q, p)
                        // F(q, m) -= D(p, n) * I(q, p, m, n)
                        K_MQ[iM * ldX3 + iQ] += (flag2 + flag6) *
                                0.5 * D_NP[iN * ldX3 + iP] * I;

                        // F(n, q) -= D(m, p) * I(n, m, q, p)
                        // F(q, n) -= D(p, m) * I(q, p, n, m)
                        K_NQ[iN * ldX3 + iQ] += (flag4 + flag7) *
                               0.5 * D_MP[iM * ldX3 + iP] * I;
                    }
                }
            }
        }
    }
}


void load_local_bufD (PFock_t pfock)
{
    int A;
    int B;
    int lo[2];
    int hi[2];
    int posrow;
    int poscol;
    double *VD;
    int i; 
    int ldD1;
    int ldD2;
    int ldD3;
    int *loadrow;
    int sizerow;
    int *loadcol;
    int sizecol;

    loadrow = pfock->loadrow;
    loadcol = pfock->loadcol;
    sizerow = pfock->sizeloadrow;
    sizecol = pfock->sizeloadcol;

    for (i = 0; i < pfock->numdmat2; i++)
    {
        NGA_Access (pfock->ga_bufD1[i],
                    pfock->lo_D1, pfock->hi_D1,
                    &(pfock->D1[i]), &ldD1);
        NGA_Access (pfock->ga_bufD2[i],
                    pfock->lo_D2, pfock->hi_D2,
                    &(pfock->D2[i]), &ldD2);
        NGA_Access (pfock->ga_bufD3[i],
                    pfock->lo_D3, pfock->hi_D3,
                    &(pfock->D3[i]), &ldD3);
    }
        
    /* D1 */
    lo[0] = pfock->sfunc_row;
    hi[0] = pfock->efunc_row;
    for (A = 0; A < sizerow; A++)
    {
        lo[1] = loadrow[PLEN * A + P_LO];
        hi[1] = loadrow[PLEN * A + P_HI];
        posrow = loadrow[PLEN * A + P_W];
        for (i = 0; i < pfock->numdmat2; i++)
        {
            VD = pfock->D1[i];
            MY_GET (pfock->ga_D[i], lo, hi, &(VD[posrow]), &ldD1);
        }
    }

    /* D2 */
    lo[0] = pfock->sfunc_col;
    hi[0] = pfock->efunc_col;
    for (B = 0; B < sizecol; B++)
    {
        lo[1] = loadcol[PLEN * B + P_LO];
        hi[1] = loadcol[PLEN * B + P_HI];
        poscol = loadcol[PLEN * B + P_W];
        for (i = 0; i < pfock->numdmat2; i++)
        {
            VD = pfock->D2[i];
            MY_GET (pfock->ga_D[i], lo, hi, &(VD[poscol]), &ldD2);
        }
    }

    /* D3 */
    for (A = 0; A < sizerow; A++)
    {
        lo[0] = loadrow[PLEN * A + P_LO];
        hi[0] = loadrow[PLEN * A + P_HI];
        posrow = loadrow[PLEN * A + P_W];
        for (B = 0; B < sizecol; B++)
        {
            lo[1] = loadcol[PLEN * B + P_LO];
            hi[1] = loadcol[PLEN * B + P_HI];
            poscol = loadcol[PLEN * B + P_W];
            for (i = 0; i < pfock->numdmat2; i++)
            {
                VD = pfock->D3[i];
                MY_GET (pfock->ga_D[i], lo, hi, &(VD[posrow * ldD3 + poscol]), &ldD3);
            }
        }
    }

    for (i = 0; i < pfock->numdmat2; i++)
    {
        NGA_Release_update (pfock->ga_bufD1[i], pfock->lo_D1, pfock->hi_D1);
        NGA_Release_update (pfock->ga_bufD2[i], pfock->lo_D2, pfock->hi_D2);
        NGA_Release_update (pfock->ga_bufD3[i], pfock->lo_D3, pfock->hi_D3);
    }

    NGA_Sync ();
}


void store_local_bufF (PFock_t pfock)
{
    int A;
    int B;
    int lo[2];
    int hi[2];
    int posrow;
    int poscol;
    double *VF;
    int i;
    double done = 1.0;
    int *ga_J;
    int *ga_K;
#ifdef _SCF_
    ga_J = pfock->ga_F;
    ga_K = pfock->ga_F;
#else
    ga_J = pfock->ga_J;
    ga_K = pfock->ga_K;
#endif
    int ldF1;
    int ldF2;
    int ldF3;
    int *loadrow;
    int sizerow;
    int *loadcol;
    int sizecol;

    ldF1 = pfock->maxrowsize;
    ldF2 = pfock->maxcolsize;
    ldF3 = ldF2;
    loadrow = pfock->loadrow;
    loadcol = pfock->loadcol;
    sizerow = pfock->sizeloadrow;
    sizecol = pfock->sizeloadcol;

    for (i = 0; i < pfock->numdmat2; i++)
    {
        if (i < pfock->numdmat)
        {
            NGA_Access (pfock->ga_bufF1[i],
                        pfock->lo_D1, pfock->hi_D1,
                        &(pfock->F1[i]), &ldF1);
            NGA_Access (pfock->ga_bufF2[i],
                        pfock->lo_D2, pfock->hi_D2,
                        &(pfock->F2[i]), &ldF2);
        }       
        NGA_Access (pfock->ga_bufF3[i],
                    pfock->lo_D3, pfock->hi_D3,
                    &(pfock->F3[i]), &ldF3);
    }
    
    /* F1 */
    lo[0] = pfock->sfunc_row;
    hi[0] = pfock->efunc_row;
    for (A = 0; A < sizerow; A++)
    {
        lo[1] = loadrow[PLEN * A + P_LO];
        hi[1] = loadrow[PLEN * A + P_HI];
        posrow = loadrow[PLEN * A + P_W];
        for (i = 0; i < pfock->numdmat2; i++)
        {
            VF = pfock->F1[i];
            MY_ACC (ga_J[i], lo, hi, &(VF[posrow]), &ldF1, &done);
        }
    }

    /* D2 */
    lo[0] = pfock->sfunc_col;
    hi[0] = pfock->efunc_col;
    for (B = 0; B < sizecol; B++)
    {
        lo[1] = loadcol[PLEN * B + P_LO];
        hi[1] = loadcol[PLEN * B + P_HI];
        poscol = loadcol[PLEN * B + P_W];
        for (i = 0; i < pfock->numdmat2; i++)
        {
            VF = pfock->F2[i];
            MY_ACC (ga_J[i], lo, hi, &(VF[poscol]), &ldF2, &done);
        }
    }

    /* D3 */
    for (A = 0; A < sizerow; A++)
    {
        lo[0] = loadrow[PLEN * A + P_LO];
        hi[0] = loadrow[PLEN * A + P_HI];
        posrow = loadrow[PLEN * A + P_W];
        for (B = 0; B < sizecol; B++)
        {
            lo[1] = loadcol[PLEN * B + P_LO];
            hi[1] = loadcol[PLEN * B + P_HI];
            poscol = loadcol[PLEN * B + P_W];
            for (i = 0; i < pfock->numdmat2; i++)
            {
                VF = pfock->F3[i];
                MY_ACC (ga_K[i], lo, hi, &(VF[posrow * ldF3 + poscol]), &ldF3, &done);
            }
        }
    }

    for (i = 0; i < pfock->numdmat2; i++)
    {
        if (i < pfock->numdmat)
        {
            NGA_Release (pfock->ga_bufF1[i], pfock->lo_D1, pfock->hi_D1);
            NGA_Release (pfock->ga_bufF2[i], pfock->lo_D2, pfock->hi_D2);
        }
        NGA_Release (pfock->ga_bufF3[i], pfock->lo_D3, pfock->hi_D3);
    }

    NGA_Sync ();
}


void compute_FD_ptr (PFock_t pfock, int startM, int endM,
                     int *ptrrow, int *rowsize)
{
    int A;
    int B;
    int i;
    int start;
    int end;
    int flag;

    for (A = 0; A < pfock->nshells; A++)
    {
        ptrrow[A] = -1;
    }    
    // init row pointers
    for (A = startM; A <= endM; A++)
    {
        start = pfock->shellptr[A];
        end = pfock->shellptr[A + 1]; 
        for (i = start; i < end; i++)
        {
            B = pfock->shellid[i];
            ptrrow[B] = 1;
        }
    }
    for (i = 0; i < pfock->natoms; i++)
    {
        start = pfock->s_startind[i];
        end = pfock->s_startind[i + 1];
        flag = -1;
        for (A = start; A < end; A++)
        {
            if (ptrrow[A] != -1)
                flag = 1;
        }
        for (A = start; A < end; A++)
        {
            ptrrow[A] = flag;
        }
    }
    *rowsize = 0;
    for (A = 0; A < pfock->nshells; A++)
    {
        if (ptrrow[A] == 1)
        {
            ptrrow[A] = *rowsize;           
            *rowsize += pfock->f_startind[A + 1] - pfock->f_startind[A];
        }
    }
}


void init_FD_load (PFock_t pfock, int *ptrrow,
                   int **loadrow, int *loadsize)
{
    int loadcount;
    int A;
    int idx;
    int lo;
    int hi;
    
    loadcount = 0;
    for (A = 0; A < pfock->nshells; A++)
    {
        if (ptrrow[A] != -1)
        {
            while (ptrrow[A] != -1 && A < pfock->nshells)
            {
                A++;
            }           
            loadcount++;
        }
    }
    *loadrow = (int *)PFOCK_MALLOC (sizeof(int) * PLEN * loadcount);
    assert (NULL != *loadrow);
    *loadsize = loadcount;
    
    loadcount = 0;
    for (A = 0; A < pfock->nshells; A++)
    {
        idx = ptrrow[A];
        if (idx != -1)
        {
            lo = pfock->f_startind[A];
            while (ptrrow[A] != -1 && A < pfock->nshells)
            {
                A++;
            }           
            hi = pfock->f_startind[A] - 1;
            (*loadrow)[loadcount * PLEN + P_LO] = lo;
            (*loadrow)[loadcount * PLEN + P_HI] = hi;
            (*loadrow)[loadcount * PLEN + P_W] = idx;
            loadcount++;
        }
    }
}


// for SCF, J = K
void fock_task (PFock_t pfock, BasisSet_t basis, int startrow, int startcol,
                int startM, int endM, int startP, int endP,
                double **D1, double **D2, double **D3,
                double ***VJ1, double ***VJ2, double ***VK3,
                int ldX1, int ldX2, int ldX3,
                double *nsq, double *nitl,double *inttime)
{
    int start;
    int end;    
    int flagsymm;    
    int ii;
    
    *nsq = 0.0;
    *nitl = 0.0;
    flagsymm = (pfock->nosymm == 1 ? 0 : 1);
    start = pfock->shellptr[startM];
    end = pfock->shellptr[endM + 1];

    #pragma omp parallel
    {
        int jj;
        double value1;  
        int start2;
        int end2;       
        double value2;
        int dimM;
        int dimN;
        int dimP;
        int dimQ;
        int flag1;
        int flag2;
        int flag3;
        int iX1M;
        int iX2P;
        int iX3M;       
        int iX3N;
        int iX3P;
        int iX3Q;
        int iX1N;
        int iX2Q;
        int iMN;
        int iPQ;
        int iMP;
        int iNP;
        int iMQ;
        int iNQ;
        int M;
        int N;
        int P;
        int Q;
        ERD_t erd;
        double *integrals;
        int nints;
        int nt;
        double mynsq;
        double mynitl;
        double **J1;
        double **J2;
        double **K3;
        
        nt = omp_get_thread_num ();

        erd = pfock->erd[nt];
        J1 = VJ1[0];
        J2 = VJ2[nt];
        K3 = VK3[nt];
        mynsq = 0;
        mynitl = 0;

        #pragma omp for schedule(dynamic)

        for (ii = start; ii < end; ii++)
        {
            M = pfock->shellrid[ii];
            dimM = pfock->f_startind[M + 1] - pfock->f_startind[M];
            iX3M = pfock->rowpos[M]; 
            iX1M = pfock->f_startind[M] - pfock->f_startind[startrow];      
            N = pfock->shellid[ii];
            value1 = pfock->shellvalue[ii];
            dimN = pfock->f_startind[N + 1] - pfock->f_startind[N];
            flag1 = (value1 < 0.0) ? 1 : 0;
            iX1N = iX3N = pfock->rowptr[ii];
            iMN = iX1M * ldX1+ iX1N;
            for (P = startP; P <= endP; P++)
            {                               
                start2 = pfock->shellptr[P];
                end2 = pfock->shellptr[P + 1];
                dimP = pfock->f_startind[P + 1] - pfock->f_startind[P];
                if ((M > P && (M + P) % 2 == 1) || 
                    (M < P && (M + P) % 2 == 0))                
                    continue;
                iX3P = pfock->colpos[P];
                iX2P = pfock->f_startind[P] - pfock->f_startind[startcol];
                iMP = iX3M * ldX3 + iX3P;
                iNP = iX3N * ldX3 + iX3P;
                for (jj = start2; jj < end2; jj++)
                {                
                    Q = pfock->shellid[jj];
                    if ((M == P) &&
                        ((N > Q && (N + Q) % 2 == 1) ||
                         (N < Q && (N + Q) % 2 == 0)))
                        continue;
                    flag3 = (M == P && Q == N) ? 0 : 1;
                    value2 = pfock->shellvalue[jj];
                    dimQ =  pfock->f_startind[Q + 1] - pfock->f_startind[Q];
                    flag2 = (value2 < 0.0) ? 1 : 0; 
                    iX3Q = iX2Q = pfock->colptr[jj];
                    iPQ = iX2P * ldX2+ iX2Q;
                    iMQ = iX3M * ldX3 + iX3Q;
                    iNQ = iX3N * ldX3 + iX3Q;
                    if (fabs(value1 * value2) >= pfock->tolscr2)
                    {
                        mynsq += 1.0;
                        mynitl += dimM*dimN*dimP*dimQ;
                        struct timeval tv3,tv4;
                        gettimeofday (&tv3,NULL);
                        CInt_computeShellQuartet (basis, erd, M, N, P, Q, &integrals, &nints);
                        gettimeofday (&tv4, NULL);
                        *inttime += (tv4.tv_sec - tv3.tv_sec) +
                                    (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
                        if (nints != 0)
                        {
                            update_F (integrals, pfock->numdmat2,
                                      dimM, dimN, dimP, dimQ,
                                      flag1, flag2, flag3, flagsymm,
                                      iMN, iPQ, iMP, iNP, iMQ, iNQ,
                                      D1, D2, D3, J1, J2, K3,
                                      ldX1, ldX2, ldX3);
                        }
                    }
                }
            }
        }


        omp_set_lock (&(pfock->lock));
        *nsq += mynsq;
        *nitl += mynitl;
        omp_unset_lock (&(pfock->lock));        
    } /* #pragma omp parallel */
}


void correct_F (PFock_t pfock)
{
    int lo[2];
    int hi[2];
    int lo2[2];
    int hi2[2];
    int ldF;
    int ld2;   
    int i;
    int j;
    int k;    
    int myrank;
    
    int nfuncs_row = pfock->nfuncs_row;
    int nfuncs_col = pfock->nfuncs_col;
    double *F_block;
    double *FT_block = pfock->FT_block;
#ifndef _SCF_
    double *K_block;
    double *KT_block = pfock->FT_block2;    
    int ldK;
#endif

    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    NGA_Distribution (pfock->ga_D[0], myrank, lo, hi);
    lo2[0] = lo[1];
    lo2[1] = lo[0];
    hi2[0] = hi[1];
    hi2[1] = hi[0];
    ld2 = hi2[1] - lo2[1] + 1;
    for (k = 0; k < pfock->numdmat2; k++)
    {
        // ga_F = ga_J when _SCF_ is defined
        NGA_Access (pfock->ga_F[k], lo, hi, &F_block, &ldF);
        NGA_Get (pfock->ga_F[k], lo2, hi2, FT_block, &ld2);
    #ifndef _SCF_
        NGA_Access (pfock->ga_K[k], lo, hi, &K_block, &ldK);
        NGA_Get (pfock->ga_K[k], lo2, hi2, KT_block, &ld2);
    #endif
        NGA_Sync ();
        for (i = 0; i < nfuncs_row; i++)
        {
            for (j = 0; j < nfuncs_col; j++)
            {
                F_block[i * ldF + j] += FT_block[j * nfuncs_row + i];
            #ifndef _SCF_    
                K_block[i * ldK + j] += KT_block[j * nfuncs_row + i];
            #endif
            }
        }
        NGA_Release_update (pfock->ga_F[k], lo, hi);
    #ifndef _SCF_    
        NGA_Release_update (pfock->ga_K[k], lo, hi);
    #endif
    }
    NGA_Sync ();
}

//RMR removed the inline because it was giving me trouble because the declaration was not in the header file
void access_bufD_GArrays (PFock_t pfock)
{
    int i;
    int ldD;

    for (i = 0; i < pfock->numdmat2; i++)
    {
        NGA_Access (pfock->ga_bufD1[i],
                    pfock->lo_D1, pfock->hi_D1,
                    &(pfock->D1[i]), &ldD);
        NGA_Access (pfock->ga_bufD2[i],
                    pfock->lo_D2, pfock->hi_D2,
                    &(pfock->D2[i]), &ldD);
        NGA_Access (pfock->ga_bufD3[i],
                    pfock->lo_D3, pfock->hi_D3,
                    &(pfock->D3[i]), &ldD);
    }
}

//RMR removed inline because it was giving problems because the declaration was not in the header file
void release_bufD_GArrays (PFock_t pfock)
{
    int i;
    for (i = 0; i < pfock->numdmat2; i++)
    {
        NGA_Release (pfock->ga_bufD1[i], pfock->lo_D1, pfock->hi_D1);
        NGA_Release (pfock->ga_bufD2[i], pfock->lo_D2, pfock->hi_D2);
        NGA_Release (pfock->ga_bufD3[i], pfock->lo_D3, pfock->hi_D3);
    }
}
