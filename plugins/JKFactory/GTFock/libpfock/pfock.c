 #include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <ga.h>
#include <macdecls.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#include <mkl.h>

#include "pfock_def.h"
#include "config.h"
#include "fock_task.h"
#include "taskq.h"
#include "screening.h"
#include "one_electron.h"


static void transform_D (PFock_t pfock)
{
    int lo[2];
    int hi[2];
    int lo2[2];
    int hi2[2];
    int ld;
    int ld2;
    double *D_block;
    double *DT_block;
    int i;
    int j;
    int k;
    int nfuncs_row;
    int nfuncs_col;
    int myrank;
    int numdmat;
    
    nfuncs_row = pfock->nfuncs_row;
    nfuncs_col = pfock->nfuncs_col;
    DT_block = pfock->FT_block;
    numdmat = pfock->numdmat;

    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    NGA_Distribution (pfock->ga_D[0], myrank, lo, hi);
    lo2[0] = lo[1];
    lo2[1] = lo[0];
    hi2[0] = hi[1];
    hi2[1] = hi[0];
    ld2 = hi2[1] - lo2[1] + 1;
    for (k = 0; k < numdmat; k++)
    {
        NGA_Access (pfock->ga_D[k + numdmat], lo, hi, &D_block, &ld);
        NGA_Get (pfock->ga_D[k], lo2, hi2, DT_block, &ld2);
        NGA_Sync ();
        for (i = 0; i < nfuncs_row; i++)
        {
            for (j = 0; j < nfuncs_col; j++)
            {
                D_block[i * nfuncs_col + j] = DT_block[j * nfuncs_row + i];
            }
        }
        NGA_Release_update (pfock->ga_D[k + numdmat], lo, hi);
    }

}


static void transform_F (PFock_t pfock)
{
    int lo[2];
    int hi[2];
    int lo2[2];
    int hi2[2];
    int ld;
    int ld2;
    double *F_block;
    double *FT_block;
    int i;
    int j;
    int k;
    int nfuncs_row;
    int nfuncs_col;
    int myrank;
    int numdmat;
    
    nfuncs_row = pfock->nfuncs_row;
    nfuncs_col = pfock->nfuncs_col;
    FT_block = pfock->FT_block;
    numdmat = pfock->numdmat;
#ifndef _SCF_
    double *K_block;
    double *KT_block = pfock->FT_block2;
#endif

    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    NGA_Distribution (pfock->ga_D[0], myrank, lo, hi);
    lo2[0] = lo[1];
    lo2[1] = lo[0];
    hi2[0] = hi[1];
    hi2[1] = hi[0];
    ld2 = hi2[1] - lo2[1] + 1;
    for (k = 0; k < numdmat; k++)
    {
        NGA_Access (pfock->ga_F[k], lo, hi, &F_block, &ld);
        NGA_Get (pfock->ga_F[k + numdmat], lo2, hi2, FT_block, &ld2);
    #ifndef _SCF_
        NGA_Access (pfock->ga_K[k], lo, hi, &K_block, &ld);
        NGA_Get (pfock->ga_K[k + numdmat], lo2, hi2, KT_block, &ld2);
    #endif
        NGA_Sync ();
        for (i = 0; i < nfuncs_row; i++)
        {
            for (j = 0; j < nfuncs_col; j++)
            {
                F_block[i * nfuncs_col + j] += FT_block[j * nfuncs_row + i];
            #ifndef _SCF_    
                K_block[i * nfuncs_col + j] += KT_block[j * nfuncs_row + i];
            #endif
            }
        }
        NGA_Release_update (pfock->ga_F[k], lo, hi);
    }
    GA_Sync ();
}


static PFockStatus_t init_fock (PFock_t pfock)
{
    int nshells;
    int nprow;
    int npcol;
    int i;
    int j;
    int n0;
    int n1;
    int t;
    int n2;
    int myrank;
    int nbp_row;
    int nbp_col;
    int nbp_p;
    int nshells_p;
        
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

    nbp_p = pfock->nbp_p;
    nbp_row = pfock->nprow * nbp_p;
    nbp_col = pfock->npcol *nbp_p;
    nshells = pfock->nshells;
    // partition task blocks
    nprow = pfock->nprow;
    npcol = pfock->npcol;
    pfock->rowptr_f = (int *)PFOCK_MALLOC (sizeof(int) * (nprow + 1));
    pfock->colptr_f = (int *)PFOCK_MALLOC (sizeof(int) * (npcol + 1));
    pfock->rowptr_sh = (int *)PFOCK_MALLOC (sizeof(int) * (nprow + 1));
    pfock->colptr_sh = (int *)PFOCK_MALLOC (sizeof(int) * (npcol + 1));
    pfock->rowptr_blk = (int *)PFOCK_MALLOC (sizeof(int) * (nprow + 1));
    pfock->colptr_blk = (int *)PFOCK_MALLOC (sizeof(int) * (npcol + 1));
    if (NULL == pfock->rowptr_f || NULL == pfock->colptr_f ||
        NULL == pfock->rowptr_sh || NULL == pfock->colptr_sh ||
        NULL == pfock->rowptr_blk || NULL == pfock->colptr_blk)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    // for row partition
    n0 = nshells/nprow;
    t = nshells%nprow;
    n1 = (nshells + nprow - 1)/nprow;    
    n2 = n1 * t;
    for (i = 0; i < nprow; i++)
    {
        pfock->rowptr_blk[i] = nbp_p * i;
        pfock->rowptr_sh[i] = i < t ? n1 * i : n2 + (i - t) * n0;
        pfock->rowptr_f[i] = pfock->f_startind[pfock->rowptr_sh[i]];
    }
    pfock->rowptr_blk[i] = nbp_row;
    pfock->rowptr_sh[i] = nshells;
    pfock->rowptr_f[i] = pfock->nbf;
    // set own
    pfock->sblk_row = pfock->rowptr_blk[myrank/npcol];
    pfock->eblk_row = pfock->rowptr_blk[myrank/npcol + 1] - 1;
    pfock->nblks_row = pfock->eblk_row - pfock->sblk_row + 1;    
    pfock->sshell_row = pfock->rowptr_sh[myrank/npcol];
    pfock->eshell_row = pfock->rowptr_sh[myrank/npcol + 1] - 1;
    pfock->nshells_row = pfock->eshell_row - pfock->sshell_row + 1;    
    pfock->sfunc_row = pfock->rowptr_f[myrank/npcol];
    pfock->efunc_row = pfock->rowptr_f[myrank/npcol + 1] - 1;
    pfock->nfuncs_row = pfock->efunc_row - pfock->sfunc_row + 1;   
    // for col partition
    n0 = nshells/npcol;
    t = nshells%npcol;
    n1 = (nshells + npcol - 1)/npcol;    
    n2 = n1 * t;
    for (i = 0; i < npcol; i++)
    {
        pfock->colptr_blk[i] = nbp_p * i;
        pfock->colptr_sh[i] = i < t ? n1 * i : n2 + (i - t) * n0;
        pfock->colptr_f[i] = pfock->f_startind[pfock->colptr_sh[i]];
    }
    pfock->colptr_blk[i] = nbp_col;
    pfock->colptr_sh[i] = nshells;
    pfock->colptr_f[i] = pfock->nbf;    
    // set own
    pfock->sblk_col = pfock->colptr_blk[myrank%npcol];
    pfock->eblk_col = pfock->colptr_blk[myrank%npcol + 1] - 1;
    pfock->nblks_col = pfock->eblk_col - pfock->sblk_col + 1;    
    pfock->sshell_col = pfock->colptr_sh[myrank%npcol];
    pfock->eshell_col = pfock->colptr_sh[myrank%npcol + 1] - 1;
    pfock->nshells_col = pfock->eshell_col - pfock->sshell_col + 1;    
    pfock->sfunc_col = pfock->colptr_f[myrank%npcol];
    pfock->efunc_col = pfock->colptr_f[myrank%npcol + 1] - 1;
    pfock->nfuncs_col = pfock->efunc_col - pfock->sfunc_col + 1;
     
    pfock->ntasks = nbp_p * nbp_p;
    pfock->blkrowptr_sh = (int *)PFOCK_MALLOC (sizeof(int) * (nbp_row + 1));
    pfock->blkcolptr_sh = (int *)PFOCK_MALLOC (sizeof(int) * (nbp_col + 1));
    if (NULL == pfock->blkrowptr_sh || NULL == pfock->blkcolptr_sh)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }

    // tasks 2D partitioning
    // row
    for (i = 0; i < nprow; i++)
    {
        nshells_p = pfock->rowptr_sh[i + 1] - pfock->rowptr_sh[i];
        n0 = nshells_p/nbp_p;
        t = nshells_p%nbp_p;
        n1 = (nshells_p + nbp_p - 1)/nbp_p;    
        n2 = n1 * t;
        for (j = 0; j < nbp_p; j++)
        {
            pfock->blkrowptr_sh[i *nbp_p + j] = pfock->rowptr_sh[i] +
                (j < t ? n1 * j : n2 + (j - t) * n0);
        }
    }
    pfock->blkrowptr_sh[i * nbp_p] = nshells;
    // col
    for (i = 0; i < npcol; i++)
    {
        nshells_p = pfock->colptr_sh[i + 1] - pfock->colptr_sh[i];
        n0 = nshells_p/nbp_p;
        t = nshells_p%nbp_p;
        n1 = (nshells_p + nbp_p - 1)/nbp_p;    
        n2 = n1 * t;
        for (j = 0; j < nbp_p; j++)
        {
            pfock->blkcolptr_sh[i *nbp_p + j] = pfock->colptr_sh[i] +
                (j < t ? n1 * j : n2 + (j - t) * n0);
        }
    }
    pfock->blkcolptr_sh[i * nbp_p] = nshells;
 
    // for correct_F
    pfock->FT_block = (double *)PFOCK_MALLOC (sizeof(double) *
        pfock->nfuncs_row * pfock->nfuncs_col);
    if (NULL == pfock->FT_block)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
#ifndef _SCF_    
    pfock->FT_block2 = (double *)PFOCK_MALLOC (sizeof(double) *
        pfock->nfuncs_row * pfock->nfuncs_col);
    if (NULL == pfock->FT_block2)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }    
#endif

    return PFOCK_STATUS_SUCCESS;
}


static PFockStatus_t init_GA (PFock_t pfock)
{
    int nbf;
    int nprow;
    int npcol;
    int *map;
    int i;
    int dims[2];
    int block[2];
    char str[8];
    int flag;

    // create global arrays
    nbf = pfock->nbf;
    nprow = pfock->nprow;
    npcol = pfock->npcol;
    map = (int *)PFOCK_MALLOC (sizeof(int) * (nprow + npcol));
    if (NULL == map)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    
    for (i = 0; i < nprow; i++)
    {       
        map[i] = pfock->rowptr_f[i];
    }   
    for (i = 0; i < npcol; i++)
    {
        map[i + nprow] = pfock->colptr_f[i];
    } 
    dims[0] = nbf;
    dims[1] = nbf;
    block[0] = nprow;
    block[1] = npcol; 

    pfock->ga_D = (int *)malloc (sizeof(int) * pfock->maxnumdmat2);
    pfock->ga_F = (int *)malloc (sizeof(int) * pfock->maxnumdmat2);  
    pfock->ga_J = (int *)malloc (sizeof(int) * pfock->maxnumdmat2);
    pfock->ga_K = (int *)malloc (sizeof(int) * pfock->maxnumdmat2);
    if (NULL == pfock->ga_D ||
        NULL == pfock->ga_F ||
        NULL == pfock->ga_J ||
        NULL == pfock->ga_K)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }

    flag = 1;
    sprintf (str, "D_0");
    pfock->ga_D[0] = NGA_Create_irreg (C_DBL, 2, dims, str, block, map);
    flag = flag * pfock->ga_D[0];
    for (i = 0; i < pfock->maxnumdmat2; i++)
    {
        if (i != 0)
        {
            flag = 1;
            sprintf (str, "D_%d", i);
            pfock->ga_D[i] = GA_Duplicate (pfock->ga_D[0], str);
            flag = flag * pfock->ga_D[i];
        }
    #ifdef _SCF_
        flag = 1;
        sprintf (str, "F_%d", i);
        pfock->ga_F[i] = GA_Duplicate (pfock->ga_D[0], str);
        pfock->ga_J[i] = pfock->ga_F[i];
        pfock->ga_K[i] = pfock->ga_F[i];
        flag = flag * pfock->ga_F[i];
    #else
        flag = 1;
        sprintf (str, "J_%d", i);
        pfock->ga_J[i] = GA_Duplicate (pfock->ga_D[0], str);
        flag = flag * pfock->ga_J[i];
        flag = 1;
        sprintf (str, "K_%d", i);
        pfock->ga_K[i] = GA_Duplicate (pfock->ga_D[0], str);
        flag = flag * pfock->ga_K[i];
        pfock->ga_F[i] = pfock->ga_J[i];            
    #endif
    }
    if (0 == flag)
    {
        PFOCK_PRINTF (1, "GA allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    
    pfock->gatable[PFOCK_MAT_TYPE_D] = pfock->ga_D;
    pfock->gatable[PFOCK_MAT_TYPE_F] = pfock->ga_F;
    pfock->gatable[PFOCK_MAT_TYPE_J] = pfock->ga_J;
    pfock->gatable[PFOCK_MAT_TYPE_K] = pfock->ga_K;

    return PFOCK_STATUS_SUCCESS;
}


static void clean_GA (PFock_t pfock)
{
    int i;

    for (i = 0; i < pfock->maxnumdmat2; i++)
    {
        GA_Destroy (pfock->ga_D[i]);        
    #ifdef _SCF_        
        GA_Destroy (pfock->ga_F[i]);
    #else
        GA_Destroy (pfock->ga_J[i]);
        GA_Destroy (pfock->ga_K[i]);
    #endif        
    }

    free (pfock->ga_D);
    free (pfock->ga_F);
    free (pfock->ga_J);
    free (pfock->ga_K);
}


static PFockStatus_t create_FD_GArrays (PFock_t pfock)
{
    int *map;
    int i;
    int dims[2];
    int block[2];
    char str[8];
    int maxrowfuncs;
    int maxcolfuncs;
    int maxrowsize;
    int maxcolsize;
    int myrank;

    maxrowsize = pfock->maxrowsize;
    maxcolsize = pfock->maxcolsize;
    maxrowfuncs = pfock->maxrowfuncs;
    maxcolfuncs = pfock->maxcolfuncs;

    pfock->ga_bufD1  = (int *)PFOCK_MALLOC (sizeof(int) * pfock->maxnumdmat2);
    pfock->ga_bufD2  = (int *)PFOCK_MALLOC (sizeof(int) * pfock->maxnumdmat2);
    pfock->ga_bufD3  = (int *)PFOCK_MALLOC (sizeof(int) * pfock->maxnumdmat2);
    pfock->ga_bufF1  = (int *)PFOCK_MALLOC (sizeof(int) * pfock->maxnumdmat2);
    pfock->ga_bufF2  = (int *)PFOCK_MALLOC (sizeof(int) * pfock->maxnumdmat2);
    pfock->ga_bufF3  = (int *)PFOCK_MALLOC (sizeof(int) * pfock->maxnumdmat2);
    pfock->ga_bufX1  = (int *)PFOCK_MALLOC (sizeof(int) * pfock->maxnumdmat);
    pfock->ga_bufX2  = (int *)PFOCK_MALLOC (sizeof(int) * pfock->maxnumdmat);
    pfock->ga_bufX3  = (int *)PFOCK_MALLOC (sizeof(int) * pfock->maxnumdmat);
    if (NULL == pfock->ga_bufD1 ||
        NULL == pfock->ga_bufD2 ||
        NULL == pfock->ga_bufD3 ||
        NULL == pfock->ga_bufF1 ||
        NULL == pfock->ga_bufF2 ||
        NULL == pfock->ga_bufF3 ||
        NULL == pfock->ga_bufX1 ||
        NULL == pfock->ga_bufX2 ||
        NULL == pfock->ga_bufX3)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }   
    map = (int *)PFOCK_MALLOC (sizeof(int) * (pfock->nprow + pfock->npcol));
    if (NULL == map)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    block[0] = pfock->nprow;
    block[1] = pfock->npcol;

    // for D1
    for (i = 0; i < pfock->nprow; i++)
    {
        map[i] = i * maxrowfuncs;
    }
    for (i = 0; i < pfock->npcol; i++)
    {       
        map[i + pfock->nprow] = i * maxrowsize;
    }
    dims[0] = pfock->nprow * maxrowfuncs;
    dims[1] = pfock->npcol * maxrowsize;
    sprintf (str, "bufD1_0");
    pfock->ga_bufD1[0] = NGA_Create_irreg (C_DBL, 2, dims, str, block, map);

    // for D2
    for (i = 0; i < pfock->nprow; i++)
    {       
        map[i] = maxcolfuncs * i;
    }
    for (i = 0; i < pfock->npcol; i++)
    {       
        map[i + pfock->nprow] = maxcolsize * i;
    }
    dims[0] = pfock->nprow * maxcolfuncs;
    dims[1] = pfock->npcol * maxcolsize;
    sprintf (str, "bufD2_0");
    pfock->ga_bufD2[0] = NGA_Create_irreg (C_DBL, 2, dims, str, block, map);

    // for D3
    for (i = 0; i < pfock->nprow; i++)
    {       
        map[i] =  maxrowsize * i;
    }
    for (i = 0; i < pfock->npcol; i++)
    {       
        map[i + pfock->nprow] = maxcolsize * i;
    }
    dims[0] = pfock->nprow * maxrowsize;
    dims[1] = pfock->npcol * maxcolsize;
    sprintf (str, "bufD3_0");
    pfock->ga_bufD3[0] = NGA_Create_irreg (C_DBL, 2, dims, str, block, map);

    if (0 == pfock->ga_bufD1[0] ||
        0 == pfock->ga_bufD2[0] ||
        0 == pfock->ga_bufD3[0])
    {
        PFOCK_PRINTF (1, "GA allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }

    for (i = 0; i < pfock->maxnumdmat2; i++)
    {
        if (i != 0)
        {
            sprintf (str, "bufD1_%d", i);
            pfock->ga_bufD1[i] = GA_Duplicate (pfock->ga_bufD1[0], str);
            sprintf (str, "bufD2_%d", i);
            pfock->ga_bufD2[i] = GA_Duplicate (pfock->ga_bufD2[0], str);
            sprintf (str, "bufD3_%d", i);
            pfock->ga_bufD3[i] = GA_Duplicate (pfock->ga_bufD3[0], str);
            if (0 == pfock->ga_bufD1[i] ||
                0 == pfock->ga_bufD2[i] ||
                0 == pfock->ga_bufD3[i])
            {
                PFOCK_PRINTF (1, "GA allocation failed\n");
                return PFOCK_STATUS_ALLOC_FAILED;
            }
        }
        if (i < pfock->maxnumdmat)
        {
            sprintf (str, "bufX1_%d", i);
            pfock->ga_bufX1[i] = GA_Duplicate (pfock->ga_bufD1[0], str);
            sprintf (str, "bufX2_%d", i);
            pfock->ga_bufX2[i] = GA_Duplicate (pfock->ga_bufD2[0], str);
            if (0 == pfock->ga_bufX1[i] ||
                0 == pfock->ga_bufX2[i])
            {
                PFOCK_PRINTF (1, "GA allocation failed\n");
                return PFOCK_STATUS_ALLOC_FAILED;
            }
        }
        sprintf (str, "bufX3_%d", i);
        pfock->ga_bufX3[i] = GA_Duplicate (pfock->ga_bufD3[0], str);
        if (0 == pfock->ga_bufX1[i])
        {
            PFOCK_PRINTF (1, "GA allocation failed\n");
            return PFOCK_STATUS_ALLOC_FAILED;
        }
    }

    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    NGA_Distribution (pfock->ga_bufD1[0], myrank, pfock->lo_D1, pfock->hi_D1);
    NGA_Distribution (pfock->ga_bufD2[0], myrank, pfock->lo_D2, pfock->hi_D2);
    NGA_Distribution (pfock->ga_bufD3[0], myrank, pfock->lo_D3, pfock->hi_D3);
    pfock->hi_F1[0] = pfock->lo_D1[0] + pfock->nfuncs_row - 1;
    pfock->hi_F1[1] = pfock->lo_D1[1] + pfock->sizemyrow - 1;
    pfock->hi_F2[0] = pfock->lo_D2[0] + pfock->nfuncs_col - 1;
    pfock->hi_F2[1] = pfock->lo_D2[1] + pfock->sizemycol - 1;
    pfock->hi_F3[0] = pfock->lo_D3[0] + pfock->sizemyrow - 1;
    pfock->hi_F3[1] = pfock->lo_D3[1] + pfock->sizemycol - 1;
    free (map);

    return PFOCK_STATUS_SUCCESS; 
}


static PFockStatus_t create_buffers (PFock_t pfock)
{
    int i;
    int j;
    int k;  
    int myrank;
    int maxcolsize;
    int maxrowsize;
    int maxrowfuncs;
    int maxcolfuncs;
    double totalFDsize = 0.0;
    int *ptrrow;
    int *ptrcol;
    int sh;
    int count;
    int myrow;
    int mycol;
    int nfuncs;
    
    
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    myrow = myrank/pfock->npcol;
    mycol = myrank%pfock->npcol;   
    ptrrow = (int *)malloc (sizeof(int) * pfock->nshells);
    ptrcol = (int *)malloc (sizeof(int) * pfock->nshells);
    if (NULL == ptrrow ||
        NULL == ptrcol)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;    
    }    

    // compute rowptr/pos and colptr/pos
    pfock->rowpos = (int *)PFOCK_MALLOC (sizeof(int) * pfock->nshells);
    pfock->colpos = (int *)PFOCK_MALLOC (sizeof(int) * pfock->nshells);
    pfock->rowptr = (int *)PFOCK_MALLOC (sizeof(int) * pfock->nnz);
    pfock->colptr = (int *)PFOCK_MALLOC (sizeof(int) * pfock->nnz);
    pfock->rowsize = (int *)PFOCK_MALLOC (sizeof(int) * pfock->nprow);
    pfock->colsize = (int *)PFOCK_MALLOC (sizeof(int) * pfock->npcol);
    if (NULL == pfock->rowpos  ||
        NULL == pfock->colpos  ||
        NULL == pfock->rowptr  || 
        NULL == pfock->colptr  ||
        NULL == pfock->rowsize ||
        NULL == pfock->colsize)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;    
    }   
    count = 0;
    maxrowsize = 0;
    maxrowfuncs = 0;
    for (i = 0; i < pfock->nprow; i++)
    {
        compute_FD_ptr (pfock,
                        pfock->rowptr_sh[i], pfock->rowptr_sh[i+1] - 1,
                        ptrrow, &(pfock->rowsize[i]));
        maxrowsize = pfock->rowsize[i] > maxrowsize ? pfock->rowsize[i] : maxrowsize;
        nfuncs = pfock->rowptr_f[i + 1] - pfock->rowptr_f[i];
        maxrowfuncs = nfuncs > maxrowfuncs ? nfuncs : maxrowfuncs;
        if (i == myrow)
        {
            pfock->sizemyrow = pfock->rowsize[i];
            init_FD_load (pfock, ptrrow, &(pfock->loadrow), &(pfock->sizeloadrow));  
        }
        for (j = pfock->rowptr_sh[i]; j < pfock->rowptr_sh[i+1]; j++)
        {
            pfock->rowpos[j] = ptrrow[j];
            for (k = pfock->shellptr[j]; k < pfock->shellptr[j+1]; k++)
            {
                sh = pfock->shellid[k];
                pfock->rowptr[count] = ptrrow[sh];
                count++;
            }
        }
    }
    count = 0;
    maxcolsize = 0;
    maxcolfuncs = 0;
    for (i = 0; i < pfock->npcol; i++)
    {
        compute_FD_ptr (pfock,
                        pfock->colptr_sh[i], pfock->colptr_sh[i+1] - 1,
                        ptrcol, &(pfock->colsize[i]));
        maxcolsize = pfock->colsize[i] > maxcolsize ? pfock->colsize[i] : maxcolsize;
        nfuncs = pfock->colptr_f[i + 1] - pfock->colptr_f[i];
        maxcolfuncs = nfuncs > maxcolfuncs ? nfuncs : maxcolfuncs;
        if (i == mycol)
        {
            pfock->sizemycol = pfock->colsize[i];
            init_FD_load (pfock, ptrcol, &(pfock->loadcol), &(pfock->sizeloadcol));  
        }
        for (j = pfock->colptr_sh[i]; j < pfock->colptr_sh[i+1]; j++)
        {
            pfock->colpos[j] = ptrcol[j];
            for (k = pfock->shellptr[j]; k < pfock->shellptr[j+1]; k++)
            {
                sh = pfock->shellid[k];
                pfock->colptr[count] = ptrcol[sh];
                count++;
            }
        }
    }
    free (ptrrow);
    free (ptrcol);
    pfock->maxrowsize = maxrowsize;
    pfock->maxcolsize = maxcolsize;
    pfock->maxrowfuncs = maxrowfuncs;
    pfock->maxcolfuncs = maxcolfuncs;
    if (myrank == 0)
    {
        PFOCK_INFO ("FD size (%d %d %d %d)\n",
            maxrowfuncs, maxcolfuncs, maxrowsize, maxcolsize);
    }

    // alloc memory for F and D
    pfock->VF1 = (double ***)PFOCK_MALLOC (sizeof(double **) * pfock->nthreads);
    pfock->VF2 = (double ***)PFOCK_MALLOC (sizeof(double **) * pfock->nthreads);
    pfock->VF3  = (double ***)PFOCK_MALLOC (sizeof(double **) * pfock->nthreads);
    if (NULL == pfock->VF1 ||
        NULL == pfock->VF2 ||
        NULL == pfock->VF3)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;    
    }
    pfock->D1  = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);
    pfock->D2  = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);
    pfock->D3  = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);  
    pfock->VD1 = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);
    pfock->VD2 = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);
    pfock->VD3 = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);
    pfock->F1  = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);
    pfock->F2  = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);
    pfock->F3  = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);
    if (NULL == pfock->D1  || NULL == pfock->VD1  ||
        NULL == pfock->D2  || NULL == pfock->VD2  ||
        NULL == pfock->D3  || NULL == pfock->VD3  ||
        NULL == pfock->F1  ||
        NULL == pfock->F2  ||
        NULL == pfock->F3)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    // F
    for (i = 0; i < pfock->nthreads; i++)
    {
        pfock->VF1[i] = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);
        pfock->VF2[i] = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);
        pfock->VF3[i] = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);

        if (NULL == pfock->VF1[i]  ||
            NULL == pfock->VF2[i]  ||
            NULL == pfock->VF3[i])
        {
            PFOCK_PRINTF (1, "memory allocation failed\n");
            return PFOCK_STATUS_ALLOC_FAILED;    
        }
    }

    // D
    for (i = 0; i < pfock->maxnumdmat2; i++)
    {
        pfock->VD1[i] = (double *)PFOCK_MALLOC 
            (sizeof(double) * maxrowfuncs * maxrowsize);
        pfock->VD2[i] = (double *)PFOCK_MALLOC 
            (sizeof(double) * maxcolfuncs * maxcolsize);
        pfock->VD3[i] = (double *)PFOCK_MALLOC 
            (sizeof(double) * maxrowsize  * maxcolsize);
        totalFDsize += 1.0 * sizeof(double) * maxrowsize * maxcolsize;
        totalFDsize += 1.0 * sizeof(double) * 
                (maxrowfuncs * maxrowsize + maxcolfuncs * maxcolsize);
        if (NULL == pfock->VD1[i] ||
            NULL == pfock->VD2[i] ||
            NULL == pfock->VD3[i])
        {
            PFOCK_PRINTF (1, "memory allocation failed\n");
            return PFOCK_STATUS_ALLOC_FAILED;
        }     
    }

    // X
    pfock->VX1 = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat);
    pfock->VX2 = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat);
    pfock->VX3 = (double **)PFOCK_MALLOC (sizeof(double *) * pfock->maxnumdmat2);
    for (j = 0; j < pfock->maxnumdmat2; j++)
    {
        pfock->VX3[j]  = (double *)PFOCK_MALLOC
                (sizeof(double) * maxrowsize  * maxcolsize * pfock->nthreads);
        totalFDsize += 1.0 * sizeof(double) * maxrowsize * maxcolsize * pfock->nthreads;
        if (NULL == pfock->VX3[j])
        {
            PFOCK_PRINTF (1, "memory allocation failed\n");
            return PFOCK_STATUS_ALLOC_FAILED;
        }
        if (j < pfock->maxnumdmat)
        {
            pfock->VX1[j] = (double *)PFOCK_MALLOC 
                (sizeof(double) * maxrowfuncs * maxrowsize * pfock->nthreads);       
            pfock->VX2[j] = (double *)PFOCK_MALLOC 
                (sizeof(double) * maxcolfuncs * maxcolsize * pfock->nthreads);            
            totalFDsize += 1.0 * sizeof(double) * 
                (maxrowfuncs * maxrowsize + maxcolfuncs * maxcolsize) * pfock->nthreads;
            if (NULL == pfock->VX1[j]  ||
                NULL == pfock->VX2[j])
            {
                PFOCK_PRINTF (1, "memory allocation failed\n");
                return PFOCK_STATUS_ALLOC_FAILED;
            } 
        }       
    }

    pfock->nbD1  = (ga_nbhdl_t *)PFOCK_MALLOC (sizeof(ga_nbhdl_t) * pfock->maxnumdmat2);
    pfock->nbD2  = (ga_nbhdl_t *)PFOCK_MALLOC (sizeof(ga_nbhdl_t) * pfock->maxnumdmat2);
    pfock->nbD3  = (ga_nbhdl_t *)PFOCK_MALLOC (sizeof(ga_nbhdl_t) * pfock->maxnumdmat2);
    pfock->nbF1  = (ga_nbhdl_t *)PFOCK_MALLOC (sizeof(ga_nbhdl_t) * pfock->maxnumdmat2);
    pfock->nbF2  = (ga_nbhdl_t *)PFOCK_MALLOC (sizeof(ga_nbhdl_t) * pfock->maxnumdmat2);
    pfock->nbF3  = (ga_nbhdl_t *)PFOCK_MALLOC (sizeof(ga_nbhdl_t) * pfock->maxnumdmat2);
    if (pfock->nbD1 == NULL &&
        pfock->nbD2 == NULL &&
        pfock->nbD3 == NULL &&
        pfock->nbF1 == NULL &&
        pfock->nbF2 == NULL &&
        pfock->nbF3 == NULL)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }         
    
    if (myrank == 0)
    {
        PFOCK_INFO ("fock uses %.3lf MB\n", totalFDsize/1024.0/1024.0);
    }

    return PFOCK_STATUS_SUCCESS;
}


static void destroy_buffers (PFock_t pfock)
{
    int i;
    int j;
    
    free (pfock->rowpos);
    free (pfock->colpos);
    free (pfock->rowptr);
    free (pfock->colptr);
    for (i = 0; i < pfock->maxnumdmat2; i++)
    {
        GA_Destroy (pfock->ga_bufD1[i]);
        GA_Destroy (pfock->ga_bufD2[i]);
        GA_Destroy (pfock->ga_bufD3[i]);     
        free (pfock->VD1[i]);
        free (pfock->VD2[i]);
        free (pfock->VD3[i]);
        free (pfock->VX3[i]);
        if (i < pfock->maxnumdmat)
        {
            free (pfock->VX1[i]);
            free (pfock->VX2[i]);                
        }
        GA_Destroy (pfock->ga_bufX3[i]);
        if (i < pfock->maxnumdmat)
        {
            GA_Destroy (pfock->ga_bufX1[i]);
            GA_Destroy (pfock->ga_bufX2[i]);   
        }
    }

    for (j = 0; j < pfock->nthreads; j++)
    {
        free (pfock->VF1[j]);
        free (pfock->VF2[j]);
        free (pfock->VF3[j]);
    }
    
    free (pfock->D1);
    free (pfock->D2);
    free (pfock->D3);
    free (pfock->VD1);
    free (pfock->VD2);
    free (pfock->VD3);
    free (pfock->F1);
    free (pfock->F2);
    free (pfock->F3);
    free (pfock->VF1);
    free (pfock->VF2);
    free (pfock->VF3);
    free (pfock->VX1);
    free (pfock->VX2);
    free (pfock->VX3);
    free (pfock->ga_bufD1);
    free (pfock->ga_bufD2);
    free (pfock->ga_bufD3);
    free (pfock->ga_bufF1);
    free (pfock->ga_bufF2);
    free (pfock->ga_bufF3);
    free (pfock->ga_bufX1);
    free (pfock->ga_bufX2);
    free (pfock->ga_bufX3);
    free (pfock->nbD1);
    free (pfock->nbD2);
    free (pfock->nbD3);
    free (pfock->nbF1);
    free (pfock->nbF2);
    free (pfock->nbF3);
}


PFockStatus_t PFock_create( BasisSet_t basis,
                            int nprow,
                            int npcol,
                            int ntasks,
                            int maxnumdmat,
                            int symm,
                            double tolscr,
                            PFock_t *_pfock )
{
    PFock_t pfock;
    int i;
    int nprocs;
    int myrank;
    PFockStatus_t ret;
    int minnshells;
    double t1;
    double t2;
    
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);         
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    pfock = (PFock_t)PFOCK_MALLOC (sizeof(struct PFock));    
    if (NULL == pfock)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    memset (pfock, 0, sizeof(PFock_t));

    // initialization
    if (myrank == 0)
    {
        PFOCK_INFO ("Initializing ...\n");
    }    
    pfock->maxnfuncs = CInt_getMaxShellDim (basis);
    pfock->nbf = CInt_getNumFuncs (basis);
    pfock->nshells = CInt_getNumShells (basis);
    pfock->natoms = CInt_getNumAtoms (basis);
    pfock->nosymm = (symm == 0 ? 1 : 0);
    pfock->nthreads = omp_get_max_threads ();
    omp_set_num_threads (pfock->nthreads);
    if (myrank == 0)
    {
        PFOCK_INFO ("  nthreads = %d\n", pfock->nthreads);
    }
    if (maxnumdmat <= 0)
    {
        PFOCK_PRINTF (1, "invalid number of density matrices\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }       
    pfock->maxnumdmat = maxnumdmat;
    pfock->maxnumdmat2 = (pfock->nosymm + 1) * maxnumdmat;

    // set nprow and nocol
    if (nprow <= 0 || nprow > pfock->nshells ||
        npcol <= 0 || npcol > pfock->nshells ||
        (nprow * npcol) > nprocs)
    {
        PFOCK_PRINTF (1, "invalid nprow or npcol\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
    else
    {
        pfock->nprow= nprow;
        pfock->npcol = npcol;
        pfock->nprocs = nprow * npcol;
    }
    
    // set tasks
    minnshells = (nprow > npcol ? nprow : npcol);
    minnshells = pfock->nshells/minnshells;
    if (ntasks >= minnshells)
    {
        pfock->nbp_p = minnshells;
    }
    else if (ntasks <= 0)
    {
        pfock->nbp_p = 3;
        pfock->nbp_p = MIN (pfock->nbp_p, minnshells);
    }
    else
    {
        pfock->nbp_p = ntasks;
    }
    // set screening threshold
    if (tolscr < 0.0)
    {
        PFOCK_PRINTF (1, "invalid screening threshold\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
    else
    {
        pfock->tolscr = tolscr;
        pfock->tolscr2 = tolscr * tolscr;
    }
    pfock->nbp_row = pfock->nbp_col = pfock->nbp_p;

    // functions starting positions
    pfock->f_startind = (int *)malloc (sizeof(int) * (pfock->nshells + 1));
    if (NULL == pfock->f_startind)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    for (i = 0; i < pfock->nshells; i++)
    {
        pfock->f_startind[i] = CInt_getFuncStartInd (basis, i);
    }
    pfock->f_startind[pfock->nshells] = pfock->nbf;

    // atoms starting positions
    pfock->s_startind = (int *)malloc (sizeof(int) * (pfock->natoms + 1));
    if (NULL == pfock->s_startind)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    for (i = 0; i < pfock->natoms; i++)
    {
        pfock->s_startind[i] = CInt_getAtomStartInd (basis, i);
    }
    pfock->s_startind[pfock->natoms] = pfock->nshells;

    // init comm
    if ((ret = init_fock (pfock)) != PFOCK_STATUS_SUCCESS)
    {
        return ret;
    }

    // create ERD
    pfock->erd = (ERD_t *)malloc (sizeof(ERD_t) * pfock->nthreads);
    for (i = 0; i < pfock->nthreads; i++)
    {
        CInt_createERD (basis, &(pfock->erd[i]));
    }
    
    // init scheduler
    if (init_taskq (pfock) != 0)
    {
        PFOCK_PRINTF (1, "task queue initialization failed\n");
        return PFOCK_STATUS_INIT_FAILED;
    }

    if (myrank == 0)
    {
        PFOCK_INFO ("Screening ...\n");
    }
    t1 = MPI_Wtime ();
    
    // screening
    if (init_screening (pfock, basis) != 0)
    {
        PFOCK_PRINTF (1, "screening initialization failed\n");
        return PFOCK_STATUS_INIT_FAILED;
    }

    if (schwartz_screening (pfock, 0, pfock->nshells - 1) != 0)
    {
        PFOCK_PRINTF (1, "schwartz screening failed\n");
        return PFOCK_STATUS_EXECUTION_FAILED;
    }

    t2 = MPI_Wtime ();
    if (myrank == 0)
    {
        PFOCK_INFO ("  takes %.3lf secs\n", t2 - t1);
        PFOCK_INFO ("Done\n");
    }
    
    // init global arrays
    if ((ret = init_GA (pfock)) != PFOCK_STATUS_SUCCESS)
    {
        return ret;
    }

    // create local buffers
    if ((ret = create_buffers (pfock)) != PFOCK_STATUS_SUCCESS)
    {
        return ret;
    }

    if ((ret = create_FD_GArrays (pfock)) != PFOCK_STATUS_SUCCESS)
    {
        return ret;
    }

    // statistics
    pfock->mpitime = (double *)PFOCK_MALLOC (sizeof(double) * pfock->nprocs);
    if (NULL == pfock->mpitime)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    pfock->computetime = (double *)PFOCK_MALLOC (sizeof(double) * pfock->nprocs);
    if (NULL == pfock->computetime)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    pfock->uitl = (double *)PFOCK_MALLOC (sizeof(double) * pfock->nprocs);
    if (NULL == pfock->uitl)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    pfock->usq = (double *)PFOCK_MALLOC (sizeof(double) * pfock->nprocs);
    if (NULL == pfock->usq)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    pfock->steals = (int *)PFOCK_MALLOC (sizeof(int) * pfock->nprocs);
    pfock->stealfrom = (int *)PFOCK_MALLOC (sizeof(int) * pfock->nprocs);
    if (NULL == pfock->steals ||
        NULL == pfock->stealfrom)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    // lock for statistics
    omp_init_lock (&(pfock->lock));

    *_pfock = pfock;
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_destroy (PFock_t pfock)
{
    int i;
    
    omp_destroy_lock (&(pfock->lock));
    free (pfock->blkrowptr_sh);
    free (pfock->blkcolptr_sh);
    free (pfock->rowptr_sh);
    free (pfock->colptr_sh);
    free (pfock->rowptr_f);
    free (pfock->colptr_f);
    free (pfock->rowptr_blk);
    free (pfock->colptr_blk);   
    free (pfock->usq);
    free (pfock->mpitime);
    free (pfock->computetime);
    free (pfock->steals);
    free (pfock->stealfrom);
    free (pfock->uitl);
    free (pfock->FT_block);
    free (pfock->f_startind);
    free (pfock->s_startind);

    for (i = 0; i < pfock->nthreads; i++)
    {
        CInt_destroyERD (pfock->erd[i]);
    }
    free (pfock->erd);
    
    clean_taskq (pfock);
    clean_screening (pfock);
    clean_GA (pfock);
    destroy_buffers (pfock);

    free (pfock);
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_setNumDenMat( PFock_t pfock,
                                  int numdmat )
{
    int i;
    int j;
    int maxrowfuncs;
    int maxcolfuncs;
    int maxrowsize;
    int maxcolsize;
    int numdmat2;

    maxrowfuncs = pfock->maxrowfuncs;
    maxcolfuncs = pfock->maxcolfuncs;
    maxrowsize = pfock->maxrowsize;
    maxcolsize = pfock->maxcolsize;
    
    if (pfock->committed == 1)
    {
        PFOCK_PRINTF (1, "can't change number of matrices after PFock_commitDenMats().\n");
        return PFOCK_STATUS_EXECUTION_FAILED;
    }
    if (numdmat <= 0 ||
        numdmat > pfock->maxnumdmat)
    {
        PFOCK_PRINTF (1, "invalid number of density matrices\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
    
    numdmat2 = numdmat * (pfock->nosymm + 1);
    pfock->numdmat = numdmat;
    pfock->numdmat2 = numdmat2;

    for (i = 0; i < numdmat2; i++)
    {
        if (i < numdmat)
        {
            pfock->ga_bufF1[i] = pfock->ga_bufX1[i];
            pfock->ga_bufF2[i] = pfock->ga_bufX2[i];
        }
        else
        {
            pfock->ga_bufF1[i] = pfock->ga_bufX1[i - numdmat];
            pfock->ga_bufF2[i] = pfock->ga_bufX2[i - numdmat];
        }
        pfock->ga_bufF3[i] = pfock->ga_bufX3[i];
    }

    for (j = 0; j < pfock->nthreads; j++)
    {
        for (i = 0; i < numdmat2; i++)
        {
            if (i < numdmat)
            {
                pfock->VF1[j][i] = &(pfock->VX1[i][j * maxrowfuncs * maxrowsize]);
                pfock->VF2[j][i] = &(pfock->VX2[i][j * maxcolfuncs * maxcolsize]);
            }
            else
            {
                pfock->VF1[j][i] = &(pfock->VX1[i - numdmat][j * maxrowfuncs * maxrowsize]);
                pfock->VF2[j][i] = &(pfock->VX2[i - numdmat][j * maxcolfuncs * maxcolsize]);
            }
            pfock->VF3[j][i] = &(pfock->VX3[i][j * maxrowsize  * maxcolsize]);
        }
    }
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_putDenMat( PFock_t pfock,
                               int index,
                               int rowstart,
                               int rowend,
                               int colstart,
                               int colend,
                               double *dmat,
                               int stride )
{
    int lo[2];
    int hi[2];

    if (pfock->committed == 1)
    {
        PFOCK_PRINTF (1, "can't change density matrix after PFock_commitDenMats().\n");
        return PFOCK_STATUS_EXECUTION_FAILED;
    }
    
    if (index < 0 ||
        index >= pfock->numdmat)
    {
        PFOCK_PRINTF (1, "invalid index\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
    lo[0] = rowstart;
    hi[0] = rowend;    
    lo[1] = colstart;
    hi[1] = colend;
    //Arguments are the handle of the array, {startrow,startcol},{endrow,endcol},local buffer (size of the block, not the entire matrix),
    //and the powerpoint I looked at makes it look like {nrows,ncols}, but the actual documentation says it is ndimensional-1, which means
    //they are ignoring the second arguement and the call here is correct
    NGA_Put (pfock->ga_D[index], lo, hi, dmat, &stride);
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_fillDenMat( PFock_t pfock,
                                int index,
                                double value )
{
    if (pfock->committed == 1)
    {
        PFOCK_PRINTF (1, "can't change density matrix after PFock_commitDenMats().\n");
        return PFOCK_STATUS_EXECUTION_FAILED;
    }
    
    if (index < 0 ||
        index >= pfock->numdmat)
    {
        PFOCK_PRINTF (1, "invalid index\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
    
    GA_Fill (pfock->ga_D[index], &value);
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_commitDenMats (PFock_t pfock)
{
    GA_Sync ();
    pfock->committed = 1;

    if (pfock->nosymm == 1)
    {
        transform_D (pfock);
    }
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_getMat( PFock_t pfock,
                            int index,
                            PFockMatType_t type,
                            int rowstart,
                            int rowend,
                            int colstart,
                            int colend,
                            double *mat,
                            int stride )
{
    int lo[2];
    int hi[2];
    int *ga;

    if (index < 0 ||
        index >= pfock->maxnumdmat)
    {
        PFOCK_PRINTF (1, "invalid index\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }    
#ifdef _SCF_
    if (PFOCK_MAT_TYPE_J == type || PFOCK_MAT_TYPE_K == type)
    {
        PFOCK_PRINTF (1, "invalid matrix type\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
#endif

    lo[0] = rowstart;
    hi[0] = rowend;    
    lo[1] = colstart;
    hi[1] = colend;
    // ga_F = ga_J when _SCF_ is not defined 
    ga = pfock->gatable[type];
    NGA_Get (ga[index], lo, hi, mat, &stride);
    
#ifndef _SCF_
    int ga_K;
    double *K;
    int i;
    int j;
    int sizerow;
    int sizecol;    
    if (PFOCK_MAT_TYPE_F == type)
    {
        sizerow = rowend - rowstart + 1;
        sizecol = colend - colstart + 1;
        K = (double *)malloc (sizerow*sizecol*sizeof(double));
        if (NULL == K)
        {
            PFOCK_PRINTF (1, "memory allocation failed\n");
            return PFOCK_STATUS_ALLOC_FAILED;
        }
        ga_K = pfock->ga_K[index];
        NGA_Get (ga_K, lo, hi, K, &stride);
        for (i = 0; i < sizerow; i++)
        {
            for (j = 0; j < sizecol; j++)                
            {
                mat[i * stride + j] += K[i * sizecol + j];
            }
        }
        free (K);
    }    
#endif

    return PFOCK_STATUS_SUCCESS;    
}


PFockStatus_t PFock_getLocalMatInds( PFock_t pfock,
                                     int *rowstart,
                                     int *rowend,
                                     int *colstart,
                                     int *colend )
{
    int lo[2];
    int hi[2];
    int myrank;

    
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    NGA_Distribution (pfock->ga_D[0], myrank, lo, hi);

    *rowstart = lo[0];
    *rowend = hi[0];
    *colstart = lo[1];
    *colend = hi[1];

    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_getLocalMatPtr ( PFock_t pfock,
                                     int index,
                                     PFockMatType_t type,
                                     int *rowstart,
                                     int *rowend,
                                     int *colstart,
                                     int *colend,
                                     double **mat,
                                     int *stride )
{
    int lo[2];
    int hi[2];
    int myrank;
    int *ga;

    if (index < 0 ||
        index >= pfock->maxnumdmat)
    {
        PFOCK_PRINTF (1, "invalid index\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
#ifdef _SCF_
    if (PFOCK_MAT_TYPE_J == type || PFOCK_MAT_TYPE_K == type)
#else
    if (PFOCK_MAT_TYPE_F == type)    
#endif 
    {
        PFOCK_PRINTF (1, "invalid matrix type\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
    ga = pfock->gatable[type];
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    NGA_Distribution (ga[index], myrank, lo, hi);
    NGA_Access (ga[index], lo, hi, mat, stride);
    *rowstart = lo[0];
    *rowend = hi[0];
    *colstart = lo[1];
    *colend = hi[1];
    
    return PFOCK_STATUS_SUCCESS;  
}


PFockStatus_t PFock_getMatGAHandle( PFock_t pfock,
                                    int index,
                                    PFockMatType_t type,
                                    int *ga )
{
    int *g;

    if (index < 0 ||
        index >= pfock->maxnumdmat)
    {
        PFOCK_PRINTF (1, "invalid index\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }

#ifdef _SCF_
    if (PFOCK_MAT_TYPE_J == type || PFOCK_MAT_TYPE_K == type)
#else
    if (PFOCK_MAT_TYPE_F == type)    
#endif 
    {
        PFOCK_PRINTF (1, "invalid matrix type\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }

    g = pfock->gatable[type];
    *ga = g[index];
    return PFOCK_STATUS_SUCCESS;
}


//#define __DETAILS__



PFockStatus_t PFock_computeFock (PFock_t pfock, BasisSet_t basis)
{
    double dzero = 0.0;
    int myrank;
    double usq;
    double uitl;
    double nsq;
    double nitl;
    int task;
    int my_sshellrow;
    int my_sshellcol;
    int ldX1;
    int ldX2;
    int ldX3;
    int startM;
    int endM;
    int startP;
    int endP;
    int rowid;
    int colid;
    int steals;
    int stealfrom;
    int myrow;
    int mycol;
    double **D1 = pfock->D1;
    double **D2 = pfock->D2;
    double **D3 = pfock->D3;
    double **J1 = pfock->VF1[0];
    double **J2 = pfock->VF2[0];
    double **K3 = pfock->VF3[0];
    double ***VJ1 = pfock->VF1;
    double ***VJ2 = pfock->VF2;
    double ***VK3 = pfock->VF3;
    ga_nbhdl_t *nbF1 = pfock->nbF1;
    ga_nbhdl_t *nbF2 = pfock->nbF2;
    ga_nbhdl_t *nbF3 = pfock->nbF3;
    int i;
    int k;
    int p;
    struct timeval tv1;
    struct timeval tv2;
    struct timeval tv3;
    struct timeval tv4;
    double timecomp;
    double t1;
    double t2;
    double timepass;
    int nthreads = pfock->nthreads;
    int maxrowfuncs = pfock->maxrowfuncs;
    int maxrowsize = pfock->maxrowsize;
    int maxcolfuncs = pfock->maxcolfuncs;
    int maxcolsize = pfock->maxcolsize;
    double done = 1.0;
    
#ifdef __DETAILS__
    double timetest[10];
    for (i = 0; i < 10; i++)
    {
        timetest[i] = 0.0;       
    }
#endif
    
    if (pfock->numdmat == 0)
    {
        return PFOCK_STATUS_SUCCESS;
    }   
    if (pfock->committed == 0)
    {
        PFOCK_PRINTF (1, "PFock_commitDenMats () must be called\n");
        return PFOCK_STATUS_INIT_FAILED;
    }
    pfock->committed = 0;
 
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    t1 = MPI_Wtime ();
    gettimeofday (&tv1, NULL);
    timecomp = 0.0;    

    if (myrank == 0)
    {
        PFOCK_INFO ("computing Fock matrices ...\n");
    }
    
    usq = 0.0;
    uitl = 0.0;
#ifdef __DETAILS__
    gettimeofday (&tv3, NULL);
#endif
    for (i = 0; i < pfock->numdmat2; i++)
    {
    #ifdef _SCF_
        GA_Fill (pfock->ga_F[i], &dzero);
    #else
        GA_Fill (pfock->ga_J[i], &dzero);
        GA_Fill (pfock->ga_K[i], &dzero);
    #endif
        if (i < pfock->numdmat)
        {
            GA_Fill (pfock->ga_bufF1[i], &dzero);
            GA_Fill (pfock->ga_bufF2[i], &dzero);
        }
        GA_Fill (pfock->ga_bufF3[i], &dzero);
    }   
    my_sshellrow = pfock->sshell_row;
    my_sshellcol = pfock->sshell_col;
    myrow = myrank/pfock->npcol;
    mycol = myrank%pfock->npcol;
#ifdef __DETAILS__
    gettimeofday (&tv4, NULL);
    timetest[0] += (tv4.tv_sec - tv3.tv_sec) +
                   (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
#endif

#ifdef __DETAILS__
    gettimeofday (&tv3, NULL);
#endif
    // local my D
    load_local_bufD (pfock);
    access_bufD_GArrays (pfock);
#ifdef __DETAILS__
    gettimeofday (&tv4, NULL);
    timetest[1] += (tv4.tv_sec - tv3.tv_sec) +
                   (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
#endif

#ifdef __DETAILS__
    gettimeofday (&tv3, NULL);
#endif
    // init my F
    for (i = 0; i < pfock->numdmat2; i++)
    {
        if (i < pfock->numdmat)
        {       
           #pragma omp parallel for
           for (k = 0; k < maxrowfuncs * maxrowsize; k++)
           {
               VJ1[i][0][k] = 0.0;    
           }
           #pragma omp parallel for
           for (k = 0; k < nthreads * maxcolfuncs * maxcolsize; k++)
           {
               VJ2[i][0][k] = 0.0;
           }
           #pragma omp parallel for
           for (k = 0; k < nthreads * maxrowsize * maxcolsize; k++)
           {
               VK3[i][0][k] = 0.0;
           }
        }
    }
#ifdef __DETAILS__
    gettimeofday (&tv4, NULL);
    timetest[2] += (tv4.tv_sec - tv3.tv_sec) +
                   (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
#endif
    
    ldX1 = pfock->maxrowsize;
    ldX2 = pfock->maxcolsize;
    ldX3 = pfock->maxcolsize;  
    /* own part */
    reset_taskq (pfock);
    double inttime=0;
    while ((task = taskq_next (pfock, myrow, mycol, 1)) < pfock->ntasks)
    {       
        rowid = task/pfock->nblks_col;
        colid = task%pfock->nblks_col;
        startM = pfock->blkrowptr_sh[pfock->sblk_row + rowid];
        endM = pfock->blkrowptr_sh[pfock->sblk_row + rowid + 1] - 1;
        startP = pfock->blkcolptr_sh[pfock->sblk_col + colid];
        endP = pfock->blkcolptr_sh[pfock->sblk_col + colid + 1] - 1;

        gettimeofday (&tv3, NULL);
        fock_task (pfock, basis, my_sshellrow, my_sshellcol,
                   startM, endM, startP, endP,
                   D1, D2, D3, VJ1, VJ2, VK3,
                   ldX1, ldX2, ldX3,
                   &nsq, &nitl,&inttime);
        gettimeofday (&tv4, NULL);
        timecomp += (tv4.tv_sec - tv3.tv_sec) +
                    (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
        usq += nsq;
        uitl += nitl;
    } /* own part */

    // reduction
#ifdef __DETAILS__
    gettimeofday (&tv3, NULL);
#endif
    for (i = 0; i < pfock->numdmat2; i++)
    {
        if (i < pfock->numdmat)
        {     
            #pragma omp parallel for
            for (k = 0; k < pfock->maxcolfuncs * pfock->maxcolsize; k++)
            {
                for (p = 1; p < pfock->nthreads; p++)
                {
                    VJ2[0][i][k] += VJ2[p][i][k];   
                }
                #if 0
                VJ2[0][i][k] += VJ2[1][i][k] + VJ2[2][i][k] +
                    VJ2[3][i][k] + VJ2[4][i][k] + 
                    VJ2[5][i][k] + VJ2[6][i][k] +
                    VJ2[7][i][k] + VJ2[8][i][k] +
                    VJ2[9][i][k] + VJ2[10][i][k] + VJ2[11][i][k];
                #endif
            }
        }
        #pragma omp parallel for
        for (k = 0; k < pfock->maxrowsize * pfock->maxcolsize; k++)
        {
            for (p = 1; p < pfock->nthreads; p++)
            {
                VK3[0][i][k] += VK3[p][i][k];   
            }
            #if 0
            VK3[0][i][k] += VK3[1][i][k] + VK3[2][i][k] +
                VK3[3][i][k] + VK3[4][i][k] + 
                VK3[5][i][k] + VK3[6][i][k] +
                VK3[7][i][k] + VK3[8][i][k] +
                VK3[9][i][k] + VK3[10][i][k] + VK3[11][i][k];
            #endif
        }
    }
#ifdef __DETAILS__
    gettimeofday (&tv4, NULL);
    timetest[3] += (tv4.tv_sec - tv3.tv_sec) +
                   (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
#endif

    // save results for local intergrals
    for (i = 0; i < pfock->numdmat2; i++)
    {
        if (i < pfock->numdmat)
        {
            NGA_NbAcc (pfock->ga_bufF1[i], pfock->lo_D1, pfock->hi_F1,
                       J1[i], &ldX1, &done, &(nbF1[i]));           
            NGA_NbAcc (pfock->ga_bufF2[i], pfock->lo_D2, pfock->hi_F2,
                       J2[i], &ldX2, &done, &(nbF2[i])); 
        }        
        NGA_NbAcc (pfock->ga_bufF3[i], pfock->lo_D3, pfock->hi_F3,
                   K3[i], &ldX3, &done, &(nbF3[i])); 
    }

    steals = 0;
    stealfrom = 0;
#ifdef _DYNAMIC_
    int idx;    
    // victim info
    int vpid;
    int vrow;
    int vcol;
    int vntasks;    
    int vsblk_row;
    int vsblk_col;
    int vnblks_col;
    int vnfuncs_row;
    int vnfuncs_col;
    int vsshellrow;
    int vsshellcol;
    int prevrow;
    int prevcol;
    int stealed;
    double **D1_task;
    double **D2_task;
    double **VD1 = pfock->VD1;
    double **VD2 = pfock->VD2;
    double **VD3 = pfock->VD3;
    int lo[2];
    int hi[2];
    ga_nbhdl_t *nbD1 = pfock->nbD1;
    ga_nbhdl_t *nbD2 = pfock->nbD2;
    ga_nbhdl_t *nbD3 = pfock->nbD3;
    
    prevrow = myrow;
    prevcol = mycol;
    /* steal tasks */
    for (idx = 0; idx < pfock->nprocs - 1; idx++)
    {
        vpid = (myrank + idx + 1)%pfock->nprocs;
        vrow = vpid/pfock->npcol;
        vcol = vpid%pfock->npcol;
        vsblk_row = pfock->rowptr_blk[vrow];
        vsblk_col = pfock->colptr_blk[vcol];
        vnblks_col = pfock->colptr_blk[vcol + 1] - vsblk_col;       
        vntasks = vnblks_col * (pfock->rowptr_blk[vrow + 1] - vsblk_row);
        vnfuncs_row = pfock->rowptr_f[vrow + 1] - pfock->rowptr_f[vrow];
        vnfuncs_col = pfock->colptr_f[vcol + 1] - pfock->colptr_f[vcol];
        vsshellrow = pfock->rowptr_sh[vrow];   
        vsshellcol = pfock->colptr_sh[vcol];
        stealed = 0;
        while ((task = taskq_next (pfock, vrow, vcol, 1)) < vntasks)
        {
            if (0 == stealed)
            {
            #ifdef __DETAILS__
                gettimeofday (&tv3, NULL);
            #endif
                if (vrow != prevrow && vrow != myrow)
                {
                    D1_task = VD1;
                    lo[0] = maxrowfuncs * vrow;
                    hi[0] = lo[0] + vnfuncs_row - 1;
                    lo[1] = maxrowsize * vcol;
                    hi[1] = lo[1] + pfock->rowsize[vrow] - 1;
                    for (i = 0; i < pfock->numdmat2; i++)
                    {
                        NGA_NbGet (pfock->ga_bufD1[i], lo, hi, VD1[i], &ldX1, &(nbD1[i]));
                    }
                }
                else if (vrow == myrow)
                {
                    D1_task = D1;
                }  
                if (vcol != prevcol && vcol != mycol)
                {
                    D2_task = VD2;
                    lo[0] = maxcolfuncs * vrow;
                    hi[0] = lo[0] + vnfuncs_col - 1;
                    lo[1] = maxcolsize * vcol;
                    hi[1] = lo[1] + pfock->colsize[vcol] - 1;
                    for (i = 0; i < pfock->numdmat2; i++)
                    {
                        NGA_NbGet (pfock->ga_bufD2[i], lo, hi, VD2[i], &ldX2, &(nbD2[i]));
                    }
                }
                else if (vcol == mycol)
                {
                    D2_task = D2;
                }
                lo[0] = pfock->maxrowsize * vrow;
                hi[0] = lo[0] + pfock->rowsize[vrow] - 1;
                lo[1] = pfock->maxcolsize * vcol;
                hi[1] = lo[1] + pfock->colsize[vcol] - 1;
                for (i = 0; i < pfock->numdmat2; i++)
                {
                    NGA_NbGet (pfock->ga_bufD3[i], lo, hi, VD3[i], &ldX3, &(nbD3[i]));
                }
                // wait for last NbAcc F
                for (i = 0; i < pfock->numdmat2; i++)
                {
                    if (i < pfock->numdmat)
                    {
                        NGA_NbWait (&(nbF1[i]));
                        NGA_NbWait (&(nbF2[i]));
                    }
                    NGA_NbWait (&(nbF3[i]));
                }
                // init F bufs
                for (i = 0; i < pfock->numdmat2; i++)
                {
                    if (i < pfock->numdmat)
                    {
                        #pragma omp parallel for
                        for (k = 0; k < maxrowfuncs * maxrowsize; k++)
                        {
                            VJ1[0][i][k] = 0.0;    
                        }
                        #pragma omp parallel for
                        for (k = 0; k < nthreads * maxcolfuncs * maxcolsize; k++)
                        {
                            VJ2[0][i][k] = 0.0;
                        }
                        #pragma omp parallel for
                        for (k = 0; k < nthreads * maxrowsize * maxcolsize; k++)
                        {
                            VK3[0][i][k] = 0.0;
                        }
                    }
                }
                // wait for NbGet
                if (vrow != prevrow && vrow != myrow)
                {
                    for (i = 0; i < pfock->numdmat2; i++)
                    {
                        NGA_NbWait (&(nbD1[i]));
                    }
                }
                if (vcol != prevcol && vcol != mycol)
                {
                    for (i = 0; i < pfock->numdmat2; i++)
                    {
                        NGA_NbWait (&(nbD2[i]));
                    }
                }
                for (i = 0; i < pfock->numdmat2; i++)
                {
                    NGA_NbWait (&(nbD3[i]));
                }
                stealfrom++;
            #ifdef __DETAILS__
                gettimeofday (&tv4, NULL);
                timetest[4] += (tv4.tv_sec - tv3.tv_sec) +
                   (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
            #endif
            }
            rowid = task/vnblks_col;
            colid = task%vnblks_col;
            // compute task
            startM = pfock->blkrowptr_sh[vsblk_row + rowid];
            endM = pfock->blkrowptr_sh[vsblk_row + rowid + 1] - 1;
            startP = pfock->blkcolptr_sh[vsblk_col + colid];
            endP = pfock->blkcolptr_sh[vsblk_col + colid + 1] - 1;
            gettimeofday (&tv3, NULL);
            fock_task (pfock, basis, vsshellrow, vsshellcol,
                       startM, endM, startP, endP,
                       D1_task, D2_task, VD3, VJ1, VJ2, VK3,
                       ldX1, ldX2, ldX3,
                       &nsq, &nitl,&inttime);
            gettimeofday (&tv4, NULL);
            timecomp += (tv4.tv_sec - tv3.tv_sec) +
                        (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
            usq += nsq;
            uitl += nitl;
            steals++;
            stealed = 1;
        }
        if (1 == stealed)
        {
        #ifdef __DETAILS__
            gettimeofday (&tv3, NULL);
        #endif
            for (i = 0; i < pfock->numdmat2; i++)
            {
                if (i < pfock->numdmat)
                {     
                    #pragma omp parallel for
                    for (k = 0; k < pfock->maxcolfuncs * pfock->maxcolsize; k++)
                    {
                        for (p = 1; p < pfock->nthreads; p++)
                        {
                            VJ2[0][i][k] += VJ2[p][i][k];   
                        }
                        #if 0
                        VJ2[0][i][k] += VJ2[1][i][k] + VJ2[2][i][k] +
                            VJ2[3][i][k] + VJ2[4][i][k] + 
                            VJ2[5][i][k] + VJ2[6][i][k] + 
                            VJ2[7][i][k] + VJ2[8][i][k] +
                            VJ2[9][i][k] + VJ2[10][i][k] + VJ2[11][i][k];
                        #endif
                    }
                }
                #pragma omp parallel for
                for (k = 0; k < pfock->maxrowsize * pfock->maxcolsize; k++)
                {
                    for (p = 1; p < pfock->nthreads; p++)
                    {
                        VK3[0][i][k] += VK3[p][i][k];   
                    }
                    #if 0
                    VK3[0][i][k] += VK3[1][i][k] + VK3[2][i][k] + VK3[3][i][k] + VK3[4][i][k] + 
                        VK3[5][i][k] + VK3[6][i][k] + VK3[7][i][k] + VK3[8][i][k] +
                        VK3[9][i][k] + VK3[10][i][k] + VK3[11][i][k];
                    #endif
                }
            }
        #ifdef __DETAILS__
            gettimeofday (&tv4, NULL);
            timetest[5] += (tv4.tv_sec - tv3.tv_sec) +
                           (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
        #endif

        #ifdef __DETAILS__
            gettimeofday (&tv3, NULL);
        #endif
            if (vrow != myrow)
            {
                lo[0] = maxrowfuncs * vrow;
                hi[0] = lo[0] + vnfuncs_row - 1;
                lo[1] = maxrowsize * vcol;
                hi[1] = lo[1] + pfock->rowsize[vrow] - 1;
                for (i = 0; i < pfock->numdmat; i++)
                {
                    NGA_NbAcc (pfock->ga_bufF1[i], lo, hi, J1[i], &ldX1,
                               &done, &(nbF1[i]));
                }
            }
            else
            {
                for (i = 0; i < pfock->numdmat; i++)
                {
                    NGA_NbAcc (pfock->ga_bufF1[i], pfock->lo_D1, pfock->hi_F1,
                               J1[i], &ldX1, &done, &(nbF1[i]));
                }
            }
            if (vcol != mycol)
            {
                lo[0] = maxcolfuncs * vrow;
                hi[0] = lo[0] + vnfuncs_col - 1;
                lo[1] = maxcolsize * vcol;
                hi[1] = lo[1] + pfock->colsize[vcol] - 1;
                for (i = 0; i < pfock->numdmat; i++)
                {
                    NGA_NbAcc (pfock->ga_bufF2[i], lo, hi, J2[i], &ldX2,
                               &done, &(nbF2[i]));
                }
            }
            else
            {
                for (i = 0; i < pfock->numdmat; i++)
                {
                    NGA_NbAcc (pfock->ga_bufF2[i], pfock->lo_D2, pfock->hi_F2,
                              J2[i], &ldX2, &done, &(nbF2[i]));
                }
            }
            lo[0] = pfock->maxrowsize * vrow;
            hi[0] = lo[0] + pfock->rowsize[vrow] - 1;
            lo[1] = pfock->maxcolsize * vcol;
            hi[1] = lo[1] + pfock->colsize[vcol] - 1;
            for (i = 0; i < pfock->numdmat2; i++)
            {
                NGA_NbAcc (pfock->ga_bufF3[i], lo, hi, K3[i], &ldX3,
                           &done, &(nbF3[i]));
            }
        #ifdef __DETAILS__
            gettimeofday (&tv4, NULL);
            timetest[6] += (tv4.tv_sec - tv3.tv_sec) +
                   (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
        #endif
            prevrow = vrow;
            prevcol = vcol;
        }      
    } /* steal tasks */
    
#endif /* #ifdef _DYNAMIC_ */

    // wait for last NbAcc F
    for (i = 0; i < pfock->numdmat2; i++)
    {
        if (i < pfock->numdmat)
        {
            NGA_NbWait (&(nbF1[i]));
            NGA_NbWait (&(nbF2[i]));
        }
        NGA_NbWait (&(nbF3[i]));
    }
    gettimeofday (&tv2, NULL);
    timepass = (tv2.tv_sec - tv1.tv_sec) +
               (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;
    t2 = MPI_Wtime ();
        
#ifdef __DETAILS__
    gettimeofday (&tv3, NULL);
#endif
    NGA_Sync ();
    release_bufD_GArrays (pfock);
    store_local_bufF (pfock);
#ifdef __DETAILS__    
    gettimeofday (&tv4, NULL);
    timetest[7] += (tv4.tv_sec - tv3.tv_sec) +
                   (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
#endif

#ifdef __DETAILS__
    gettimeofday (&tv3, NULL);
#endif
    // correct F
    if (pfock->nosymm == 1)
    {
        transform_F (pfock);
    }
    else
    {
        correct_F (pfock);    
    }
#ifdef __DETAILS__      
    gettimeofday (&tv4, NULL);
    timetest[8] += (tv4.tv_sec - tv3.tv_sec) +
                   (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
#endif

#ifdef __DETAILS__
    PFOCK_INFO ("  rank %7d: %9.3lf secs (%9.3lf secs), screening: %.4le (%.4le), %.3lf us/ints, %d from %d\n",
                 myrank, timepass, timecomp,
                 usq, uitl, timepass/uitl * 1000.0 * 1000.0, steals, stealfrom);
    PFOCK_INFO ("  rank %7d: %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf\n",
                 myrank, timetest[0], timetest[1], timetest[2], timetest[3], timetest[4],
                 timetest[5], timetest[6], timetest[7], timetest[8], timetest[9]);
#endif
    // statistics
    MPI_Gather (&steals, 1, MPI_INT, pfock->steals, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather (&stealfrom, 1, MPI_INT, pfock->stealfrom, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather (&usq, 1, MPI_DOUBLE, pfock->usq, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&uitl, 1, MPI_DOUBLE, pfock->uitl, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&timepass, 1, MPI_DOUBLE, pfock->mpitime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&timecomp, 1, MPI_DOUBLE, pfock->computetime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (myrank == 0)
    {
        double totaltime;
        double maxtime;
        double totalusq;
        double maxusq;
        double computetime;
        double totaluitl;
        double maxuitl;
        int totalsteals;
        int totalstealfrom;
        double tsq;      
        totaltime = 0.0;
        maxtime = 0.0;
        totalusq = 0.0;
        maxusq = 0.0;
        totaluitl = 0.0;
        maxuitl = 0.0;
        totalsteals = 0;
        totalstealfrom = 0;
        computetime = 0;
        for (i = 0; i < pfock->nprocs; i++)
        {
            totaltime += pfock->mpitime[i];
            maxtime = maxtime < pfock->mpitime[i] ? pfock->mpitime[i] : maxtime;
            totalusq += pfock->usq[i];
            maxusq = maxusq < pfock->usq[i] ? pfock->usq[i] : maxusq;
            totaluitl += pfock->uitl[i];
            maxuitl = maxuitl < pfock->uitl[i] ? pfock->uitl[i] : maxuitl;
            totalsteals += pfock->steals[i];
            totalstealfrom += pfock->stealfrom[i];
            computetime += pfock->computetime[i];
        }
        tsq = pfock->nshells;
        tsq = ((tsq + 1) * tsq/2.0 + 1) * tsq * (tsq + 1)/4.0;
        printf ("fock takes %.3lf secs (%.3lf, %.3lf: %.3lf)\n", t2 - t1,
            computetime/pfock->nprocs, totaltime/pfock->nprocs,
            computetime/totaltime);
        printf ("integrals took %.3lf secs\n",inttime/pfock->nprocs);
        printf ("total usq %.4le (%.3lf), nsq = %.4le (%.3lf), uitl %.4le (%.3lf)\n",
            totalusq, maxusq/(totalusq/pfock->nprocs),
            tsq, totalusq/tsq,
            totaluitl, maxuitl/(totaluitl/pfock->nprocs));
        PFOCK_INFO ("load balance %.3lf\n", maxtime/(totaltime/pfock->nprocs));
        PFOCK_INFO ("steals (%d %.3lf), from (%d %.3lf)\n",
            totalsteals, (double)totalsteals/pfock->nprocs,
            totalstealfrom, (double)totalstealfrom/pfock->nprocs);
    }

    NGA_Sync ();
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_createCoreHMat ( PFock_t pfock,
                                     BasisSet_t basis,
                                     CoreH_t *hmat )
{
    int startrow;
    int endrow;
    int startcol;
    int endcol;
    double *mtmp;
    int ga;
    int lo[2];
    int hi[2];
    int ld;
    CoreH_t h;

    h = (CoreH_t)malloc (sizeof(struct CoreH));
    startrow = pfock->sshell_row;
    endrow = pfock->eshell_row;
    startcol = pfock->sshell_col;
    endcol = pfock->eshell_col;
    mtmp = (double *)malloc (sizeof(double) * 
        pfock->nfuncs_row * pfock->nfuncs_col);
    if (NULL == mtmp || NULL == h)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    memset (mtmp, 0, sizeof(double) * 
        pfock->nfuncs_row * pfock->nfuncs_col);
    compute_H (pfock, basis, startrow, endrow, startcol, endcol, mtmp);

    ga = GA_Duplicate (pfock->ga_D[0], "tmp");
    if (0 == ga)
    {
        PFOCK_PRINTF (1, "GA allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }    
    lo[0] = pfock->sfunc_row;
    hi[0] = pfock->efunc_row;
    lo[1] = pfock->sfunc_col;
    hi[1] = pfock->efunc_col;
    ld = hi[1] - lo[1] + 1;
    NGA_Put (ga, lo, hi, mtmp, &ld);
    free (mtmp);
    NGA_Sync ();

    h->ga = ga;
    *hmat = h;
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_destroyCoreHMat ( CoreH_t hmat )
{
    GA_Destroy (hmat->ga);
    free (hmat);

    return PFOCK_STATUS_SUCCESS;    
}


PFockStatus_t PFock_getCoreHMatGAHandle( CoreH_t hmat,
                                         int *ga )
{
    *ga = hmat->ga;
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_getCoreHMat( CoreH_t hmat,
                                 int rowstart,
                                 int rowend,
                                 int colstart,
                                 int colend,
                                 double *mat,
                                 int stride )
{
    int lo[2];
    int hi[2];
    int ga;

    ga = hmat->ga;
    lo[0] = rowstart;
    hi[0] = rowend;    
    lo[1] = colstart;
    hi[1] = colend;    
    NGA_Get (ga, lo, hi, mat, &stride);
    return PFOCK_STATUS_SUCCESS;    
}


PFockStatus_t PFock_createOvlMat ( PFock_t pfock,
                                   BasisSet_t basis,
                                   Ovl_t *omat )
{
    int startrow;
    int endrow;
    int startcol;
    int endcol;
    double *mtmp;
    int ga;
    int lo[2];
    int hi[2];
    int ld;
    Ovl_t o;
    
    o = (Ovl_t)malloc (sizeof(struct Ovl));
    startrow = pfock->sshell_row;
    endrow = pfock->eshell_row;
    startcol = pfock->sshell_col;
    endcol = pfock->eshell_col;
    mtmp = (double *)malloc (sizeof(double) * 
        pfock->nfuncs_row * pfock->nfuncs_col);   
    if (NULL == mtmp || o == NULL)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    memset (mtmp, 0, sizeof(double) * 
        pfock->nfuncs_row * pfock->nfuncs_col);   
    compute_S (pfock, basis, startrow, endrow, startcol, endcol, mtmp);

    ga = GA_Duplicate (pfock->ga_D[0], "tmp");
    if (0 == ga)
    {
        PFOCK_PRINTF (1, "GA allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }    
    lo[0] = pfock->sfunc_row;
    hi[0] = pfock->efunc_row;
    lo[1] = pfock->sfunc_col;
    hi[1] = pfock->efunc_col;
    ld = hi[1] - lo[1] + 1;
    NGA_Put (ga, lo, hi, mtmp, &ld);
    free (mtmp);
    NGA_Sync ();

    o->ga = ga;
    *omat = o;
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_destroyOvlMat ( Ovl_t omat )
{
    GA_Destroy (omat->ga);
    free (omat);

    return PFOCK_STATUS_SUCCESS;    
}


PFockStatus_t PFock_getOvlMatGAHandle( Ovl_t omat,
                                       int *ga )
{
    *ga = omat->ga;
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_getOvlMat( Ovl_t omat,
                               int rowstart,
                               int rowend,
                               int colstart,
                               int colend,
                               double *mat,
                               int stride )
{
    int lo[2];
    int hi[2];
    int ga;

    ga = omat->ga;
    lo[0] = rowstart;
    hi[0] = rowend;    
    lo[1] = colstart;
    hi[1] = colend;    
    NGA_Get (ga, lo, hi, mat, &stride);
    return PFOCK_STATUS_SUCCESS;    
}


PFockStatus_t PFock_GAInit( int nbf,
                            int nprow,
                            int npcol,
                            int numdenmat,
                            int sizeheap,
                            int sizestack )
{
    int stack;
    int heap;
    int maxrowsize;
    int maxcolsize;
    //std::cout<<"GTFock::PFock_GAInit() with arguments: NBasis="<<nbf<<" NProccessors Per (Row,Column)=("<<nprow<<","<<npcol<<")";
    //std::cout<<" NDensMat="<<numdenmat<<" Size of (heap,stack)=("<<sizeheap<<","<<sizestack<<")"<<std::endl;
    maxrowsize = (nbf + nprow - 1)/nprow;
    maxcolsize = (nbf + npcol - 1)/npcol;    
    heap = numdenmat * 5 * maxrowsize * maxcolsize * sizeof(double);
    stack = heap;
    heap += sizeheap;
    stack += sizestack;
    //std::cout<<"The new stack and heap sizes: "<<stack<<" "<<heap<<std::endl;
#if ( _DEBUG_LEVEL_ > 2 )
    printf ("MA_init(): GA heap = %lf KB, stack = %lf KB\n",
            heap/1024.0, stack/1024.0);
#endif
    GA_Initialize ();
    //std::cout<<"Not failing in GA_Initialize...."<<std::endl;
    if (!MA_init (C_DBL, heap, stack))
    {
        printf ("MA_init() failed\n");
        return PFOCK_STATUS_INIT_FAILED;
    }
    
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_GAFinalize (void)
{    
    GA_Terminate ();
   
    return PFOCK_STATUS_SUCCESS;
}
