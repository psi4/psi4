#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <ga.h>
#include <macdecls.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#ifdef __INTEL_MKL__
  #include <mkl.h>
#endif
#include <assert.h>
#include <math.h>

#include "pfock.h"
#include "config.h"
#include "fock_task.h"
#include "fock_buf.h"
#include "taskq.h"
#include "screening.h"
#include "one_electron.h"


static PFockStatus_t init_fock(PFock_t pfock)
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
        
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    nbp_p = pfock->nbp_p;
    nbp_row = pfock->nprow * nbp_p;
    nbp_col = pfock->npcol *nbp_p;
    nshells = pfock->nshells;
    // partition task blocks
    nprow = pfock->nprow;
    npcol = pfock->npcol;
    pfock->rowptr_f = (int *)PFOCK_MALLOC(sizeof(int) * (nprow + 1));
    pfock->colptr_f = (int *)PFOCK_MALLOC(sizeof(int) * (npcol + 1));
    pfock->rowptr_sh = (int *)PFOCK_MALLOC(sizeof(int) * (nprow + 1));
    pfock->colptr_sh = (int *)PFOCK_MALLOC(sizeof(int) * (npcol + 1));
    pfock->rowptr_blk = (int *)PFOCK_MALLOC(sizeof(int) * (nprow + 1));
    pfock->colptr_blk = (int *)PFOCK_MALLOC(sizeof(int) * (npcol + 1));
    if (NULL == pfock->rowptr_f || NULL == pfock->colptr_f ||
        NULL == pfock->rowptr_sh || NULL == pfock->colptr_sh ||
        NULL == pfock->rowptr_blk || NULL == pfock->colptr_blk)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    pfock->mem_cpu += 3.0 * sizeof(int) * ((nprow + 1) + (npcol + 1));
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
    pfock->blkrowptr_sh = (int *)PFOCK_MALLOC(sizeof(int) * (nbp_row + 1));
    pfock->blkcolptr_sh = (int *)PFOCK_MALLOC(sizeof(int) * (nbp_col + 1));
    pfock->mem_cpu += sizeof(int) * ((nbp_row + 1) + (nbp_col + 1));
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
 
    return PFOCK_STATUS_SUCCESS;
}


static void recursive_bisection (int *rowptr, int first, int last,
                                 int npartitions, int *partition_ptr)
{
	int offset = rowptr[first];
	int nnz = rowptr[last] - rowptr[first];

	if(npartitions == 1)
	{
		partition_ptr[0] = first;
		return;
	}

	int left = npartitions/2;
	double ideal = ((double)nnz * (double)left)/npartitions;
	int i;
	for(i = first; i < last; i++)
	{
		double count = rowptr[i] - offset;
		double next_count = rowptr[i + 1] - offset;
		if(next_count > ideal)
		{
			if(next_count - ideal > ideal - count)
			{
				recursive_bisection(rowptr, first, i, left, partition_ptr);
				recursive_bisection(rowptr, i, last,
                                    npartitions - left, partition_ptr + left);
				return;
			}
			else
			{
				recursive_bisection(rowptr, first, i + 1, left, partition_ptr);
				recursive_bisection(rowptr, i + 1, last,
                                    npartitions - left, partition_ptr + left);
				return;
			}
		}
	}
}


static int nnz_partition (int m, int nnz, int min_nrows,
                          int *rowptr, int npartitions, int *partition_ptr)
{
	recursive_bisection(rowptr, 0, m, npartitions, partition_ptr);
	partition_ptr[npartitions] = m;

	for (int i = 0; i < npartitions; i++)
	{
		int nrows = partition_ptr[i + 1] - partition_ptr[i];
		if (nrows < min_nrows)
		{
			return -1;
		}
	}
    
	return 0;
}


static PFockStatus_t repartition_fock (PFock_t pfock)
{
    int nshells = pfock->nshells;
    int nnz = pfock->nnz;
    int nbp_p = pfock->nbp_p;
    int nprow = pfock->nprow;
    int npcol = pfock->npcol;
    int *shellptr = pfock->shellptr;
    int myrank;
    int ret;
    
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

    // for row partition
    int *newrowptr = (int *)malloc (sizeof(int) * (nprow + 1));
    int *newcolptr = (int *)malloc (sizeof(int) * (npcol + 1));
    ret = nnz_partition (nshells, nnz, nbp_p, shellptr, nprow, newrowptr);    
    if (ret != 0)
    {
        PFOCK_PRINTF (1, "nbp_p is too large\n");
        return PFOCK_STATUS_EXECUTION_FAILED;
    }
    ret = nnz_partition (nshells, nnz, nbp_p, shellptr, npcol, newcolptr);
    if (ret != 0)
    {
        PFOCK_PRINTF (1, "nbp_p is too large\n");
        return PFOCK_STATUS_EXECUTION_FAILED;
    }
    memcpy (pfock->rowptr_sh, newrowptr, sizeof(int) * (nprow + 1));    
    memcpy (pfock->colptr_sh, newcolptr, sizeof(int) * (npcol + 1));
    free (newrowptr);
    free (newcolptr);
    
    for (int i = 0; i < nprow; i++)
    {
        pfock->rowptr_f[i] = pfock->f_startind[pfock->rowptr_sh[i]];
    }
    pfock->rowptr_f[nprow] = pfock->nbf;
    // set own  
    pfock->sshell_row = pfock->rowptr_sh[myrank/npcol];
    pfock->eshell_row = pfock->rowptr_sh[myrank/npcol + 1] - 1;
    pfock->nshells_row = pfock->eshell_row - pfock->sshell_row + 1;    
    pfock->sfunc_row = pfock->rowptr_f[myrank/npcol];
    pfock->efunc_row = pfock->rowptr_f[myrank/npcol + 1] - 1;
    pfock->nfuncs_row = pfock->efunc_row - pfock->sfunc_row + 1;  
    
    // for col partition
    for (int i = 0; i < npcol; i++)
    {
        pfock->colptr_f[i] = pfock->f_startind[pfock->colptr_sh[i]];
    }
    pfock->colptr_f[npcol] = pfock->nbf;    
    // set own   
    pfock->sshell_col = pfock->colptr_sh[myrank%npcol];
    pfock->eshell_col = pfock->colptr_sh[myrank%npcol + 1] - 1;
    pfock->nshells_col = pfock->eshell_col - pfock->sshell_col + 1;    
    pfock->sfunc_col = pfock->colptr_f[myrank%npcol];
    pfock->efunc_col = pfock->colptr_f[myrank%npcol + 1] - 1;
    pfock->nfuncs_col = pfock->efunc_col - pfock->sfunc_col + 1;
     
    // tasks 2D partitioning
    // row
    int nshells_p;
    int n0;
    int t;
    int n1;
    int n2;
    for (int i = 0; i < nprow; i++)
    {
        nshells_p = pfock->rowptr_sh[i + 1] - pfock->rowptr_sh[i];
        n0 = nshells_p/nbp_p;
        t = nshells_p%nbp_p;
        n1 = (nshells_p + nbp_p - 1)/nbp_p;    
        n2 = n1 * t;
        for (int j = 0; j < nbp_p; j++)
        {
            pfock->blkrowptr_sh[i *nbp_p + j] = pfock->rowptr_sh[i] +
                (j < t ? n1 * j : n2 + (j - t) * n0);
        }
    }
    pfock->blkrowptr_sh[nprow * nbp_p] = nshells;
    // col
    for (int i = 0; i < npcol; i++)
    {
        nshells_p = pfock->colptr_sh[i + 1] - pfock->colptr_sh[i];
        n0 = nshells_p/nbp_p;
        t = nshells_p%nbp_p;
        n1 = (nshells_p + nbp_p - 1)/nbp_p;    
        n2 = n1 * t;
        for (int j = 0; j < nbp_p; j++)
        {
            pfock->blkcolptr_sh[i *nbp_p + j] = pfock->colptr_sh[i] +
                (j < t ? n1 * j : n2 + (j - t) * n0);
        }
    }
    pfock->blkcolptr_sh[npcol * nbp_p] = nshells;

    // for correct_F
    pfock->FT_block = (double *)PFOCK_MALLOC(sizeof(double) *
        pfock->nfuncs_row * pfock->nfuncs_col);
    pfock->mem_cpu += 1.0 * pfock->nfuncs_row * pfock->nfuncs_col * sizeof(double);
    if (NULL == pfock->FT_block)
    {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    
    return PFOCK_STATUS_SUCCESS;
}


static PFockStatus_t create_GA (PFock_t pfock)
{
    int nbf;
    int nprow;
    int npcol;
    int *map;
    int i;
    int dims[2];
    int block[2];
    char str[8];

    // create global arrays
    nbf = pfock->nbf;
    nprow = pfock->nprow;
    npcol = pfock->npcol;
    map = (int *)PFOCK_MALLOC(sizeof(int) * (nprow + npcol));
    if (NULL == map)
    {
        PFOCK_PRINTF(1, "memory allocation failed\n");
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

    pfock->ga_D = (int *)PFOCK_MALLOC(sizeof(int) * pfock->max_numdmat2);
    pfock->ga_F = (int *)PFOCK_MALLOC(sizeof(int) * pfock->max_numdmat2);
    pfock->ga_K = (int *)PFOCK_MALLOC(sizeof(int) * pfock->max_numdmat2);
    if (pfock->ga_D == NULL ||
        pfock->ga_F == NULL ||
        pfock->ga_K == NULL) {
        PFOCK_PRINTF(1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;        
    }

    sprintf(str, "D_0");
    pfock->ga_D[0] = NGA_Create_irreg(C_DBL, 2, dims, str, block, map);
    if (0 == pfock->ga_D[0]) {
        PFOCK_PRINTF(1, "GA allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    
    for (int i = 0; i < pfock->max_numdmat2; i++) {
        if (i != 0) {                
            sprintf(str, "D_%d", i);
            pfock->ga_D[i] = GA_Duplicate(pfock->ga_D[0], str);
            if (0 == pfock->ga_D[i]) {
                PFOCK_PRINTF(1, "GA allocation failed\n");
                return PFOCK_STATUS_ALLOC_FAILED;
            }
        }
    
        sprintf(str, "F_%d", i);
        pfock->ga_F[i] = GA_Duplicate(pfock->ga_D[0], str);
        if (0 == pfock->ga_F[i]) {
            PFOCK_PRINTF(1, "GA allocation failed\n");
            return PFOCK_STATUS_ALLOC_FAILED;
        }

        sprintf (str, "K_%d", i);
        pfock->ga_K[i] = GA_Duplicate(pfock->ga_D[0], str);
        if (0 == pfock->ga_K[i]) {
            PFOCK_PRINTF(1, "GA allocation failed\n");
            return PFOCK_STATUS_ALLOC_FAILED;
        }
    }
    pfock->gatable[PFOCK_MAT_TYPE_D] = pfock->ga_D;
    pfock->gatable[PFOCK_MAT_TYPE_F] = pfock->ga_F;
    pfock->gatable[PFOCK_MAT_TYPE_J] = pfock->ga_F;
    pfock->gatable[PFOCK_MAT_TYPE_K] = pfock->ga_K;

    PFOCK_FREE(map);
    
    return PFOCK_STATUS_SUCCESS;
}


static void destroy_GA(PFock_t pfock)
{
    for (int i = 0; i < pfock->max_numdmat2; i++) {
        GA_Destroy(pfock->ga_D[i]);        
        GA_Destroy(pfock->ga_F[i]);
        GA_Destroy(pfock->ga_K[i]);
    }
    PFOCK_FREE(pfock->ga_D);
    PFOCK_FREE(pfock->ga_F);
    PFOCK_FREE(pfock->ga_K);
    
}


static PFockStatus_t create_FD_GArrays (PFock_t pfock)
{
    int dims[2];
    int block[2];
    char str[8];
    
    int sizeD1 = pfock->sizeX1;
    int sizeD2 = pfock->sizeX2;
    int sizeD3 = pfock->sizeX3;  
    int *map = (int *)PFOCK_MALLOC(sizeof(int) * (1 + pfock->nprocs));
    if (NULL == map) {
        PFOCK_PRINTF(1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    block[0] = pfock->nprocs;
    block[1] = 1;

    // for D1
    for (int i = 0; i < pfock->nprocs; i++) {
        map[i] = i;
    }
    map[pfock->nprocs] = 0;
    dims[0] = pfock->nprocs;
    dims[1] = sizeD1;
    pfock->ga_D1 = (int *)PFOCK_MALLOC(sizeof(int) * pfock->max_numdmat2);
    if (pfock->ga_D1 == NULL) {
        PFOCK_PRINTF(1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    sprintf(str, "D1_0");
    pfock->ga_D1[0] = NGA_Create_irreg(C_DBL, 2, dims, str, block, map);
    for (int i = 0; i < pfock->max_numdmat2; i++) {
        if (i != 0) {
            sprintf(str, "D1_%d", i);
            pfock->ga_D1[i] = GA_Duplicate(pfock->ga_D1[0], str);
        }
        if (pfock->ga_D1[i] == 0) {
            PFOCK_PRINTF(1, "GA allocation failed\n");
            return PFOCK_STATUS_ALLOC_FAILED;
        }
    }
        
    // for D2
    for (int i = 0; i < pfock->nprocs; i++)
    {
        map[i] = i;
    }
    map[pfock->nprocs] = 0;
    dims[0] = pfock->nprocs;
    dims[1] = sizeD2;
    pfock->ga_D2 = (int *)PFOCK_MALLOC(sizeof(int) * pfock->max_numdmat2);
    if (pfock->ga_D2 == NULL) {
        PFOCK_PRINTF(1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    sprintf(str, "D2_0");
    pfock->ga_D2[0] = NGA_Create_irreg(C_DBL, 2, dims, str, block, map);
    for (int i = 0; i < pfock->max_numdmat2; i++) {
        if (i != 0) {
            sprintf(str, "D2_%d", i);
            pfock->ga_D2[i] = GA_Duplicate(pfock->ga_D2[0], str);
        }
        if (pfock->ga_D2[i] == 0) {
            PFOCK_PRINTF(1, "GA allocation failed\n");
            return PFOCK_STATUS_ALLOC_FAILED;
        }
    }
    
    // for D3
    for (int i = 0; i < pfock->nprocs; i++)
    {
        map[i] = i;
    }
    map[pfock->nprocs] = 0;
    dims[0] = pfock->nprocs;
    dims[1] = sizeD3;
    pfock->ga_D3 = (int *)PFOCK_MALLOC(sizeof(int) * pfock->max_numdmat2);
    if (pfock->ga_D3 == NULL) {
        PFOCK_PRINTF(1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    sprintf(str, "D3_0");
    pfock->ga_D3[0] = NGA_Create_irreg(C_DBL, 2, dims, str, block, map);
    for (int i = 0; i < pfock->max_numdmat2; i++) {
        if (i != 0) {
            sprintf(str, "D3_%d", i);
            pfock->ga_D3[i] = GA_Duplicate(pfock->ga_D3[0], str);
        }
        if (pfock->ga_D3[i] == 0) {
            PFOCK_PRINTF(1, "GA allocation failed\n");
            return PFOCK_STATUS_ALLOC_FAILED;
        }
    }

    // F1, F2, and F3
    pfock->ga_F1 = (int *)PFOCK_MALLOC(sizeof(int) * pfock->max_numdmat2);
    pfock->ga_F2 = (int *)PFOCK_MALLOC(sizeof(int) * pfock->max_numdmat2);
    pfock->ga_F3 = (int *)PFOCK_MALLOC(sizeof(int) * pfock->max_numdmat2);
    if (pfock->ga_F1 == NULL ||
        pfock->ga_F2 == NULL ||
        pfock->ga_F3 == NULL) {
        PFOCK_PRINTF(1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    for (int i = 0; i < pfock->max_numdmat2; i++) {
        sprintf(str, "F1_%d", i);
        pfock->ga_F1[i] = GA_Duplicate(pfock->ga_D1[0], str);
        sprintf(str, "F2_%d", i);
        pfock->ga_F2[i] = GA_Duplicate(pfock->ga_D2[0], str);
        sprintf(str, "F3_%d", i);
        pfock->ga_F3[i] = GA_Duplicate(pfock->ga_D3[0], str);
        if (pfock->ga_F1[i] == 0 ||
            pfock->ga_F2[i] == 0 ||
            pfock->ga_F3[i] == 0) {
            PFOCK_PRINTF(1, "GA allocation failed\n");
            return PFOCK_STATUS_ALLOC_FAILED;
        }
    }
    
    PFOCK_FREE(map);

    return PFOCK_STATUS_SUCCESS; 
}


static PFockStatus_t create_buffers (PFock_t pfock)
{
    int myrank;
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    int myrow = myrank/pfock->npcol;
    int mycol = myrank%pfock->npcol;   
    int *ptrrow = (int *)PFOCK_MALLOC(sizeof(int) * pfock->nshells);
    int *ptrcol = (int *)PFOCK_MALLOC(sizeof(int) * pfock->nshells);
    if (NULL == ptrrow ||
        NULL == ptrcol) {
        PFOCK_PRINTF(1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;    
    }    

    // compute rowptr/pos and colptr/pos
    pfock->rowpos = (int *)PFOCK_MALLOC(sizeof(int) * pfock->nshells);
    pfock->colpos = (int *)PFOCK_MALLOC(sizeof(int) * pfock->nshells);
    pfock->rowptr = (int *)PFOCK_MALLOC(sizeof(int) * pfock->nnz);
    pfock->colptr = (int *)PFOCK_MALLOC(sizeof(int) * pfock->nnz);
    pfock->rowsize = (int *)PFOCK_MALLOC(sizeof(int) * pfock->nprow);
    pfock->colsize = (int *)PFOCK_MALLOC(sizeof(int) * pfock->npcol);
    pfock->mem_cpu += 1.0 * sizeof(int) *
        (2.0 * pfock->nshells + 2.0 * pfock->nnz +
         pfock->nprow + pfock->npcol);
    if (NULL == pfock->rowpos  ||
        NULL == pfock->colpos  ||
        NULL == pfock->rowptr  || 
        NULL == pfock->colptr  ||
        NULL == pfock->rowsize ||
        NULL == pfock->colsize) {
        PFOCK_PRINTF(1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;    
    }   
    int count = 0;
    int maxrowsize = 0;
    int maxrowfuncs = 0;
    for (int i = 0; i < pfock->nprow; i++) {
        compute_FD_ptr (pfock,
                        pfock->rowptr_sh[i], pfock->rowptr_sh[i+1] - 1,
                        ptrrow, &(pfock->rowsize[i]));
        maxrowsize =
            pfock->rowsize[i] > maxrowsize ? pfock->rowsize[i] : maxrowsize;
        int nfuncs = pfock->rowptr_f[i + 1] - pfock->rowptr_f[i];
        maxrowfuncs = nfuncs > maxrowfuncs ? nfuncs : maxrowfuncs;
        if (i == myrow) {
            pfock->sizemyrow = pfock->rowsize[i];
            init_FD_load(pfock, ptrrow,
                         &(pfock->loadrow), &(pfock->sizeloadrow));  
        }
        for (int j = pfock->rowptr_sh[i]; j < pfock->rowptr_sh[i+1]; j++) {
            pfock->rowpos[j] = ptrrow[j];
            for (int k = pfock->shellptr[j]; k < pfock->shellptr[j+1]; k++)
            {
                int sh = pfock->shellid[k];
                pfock->rowptr[count] = ptrrow[sh];
                count++;
            }
        }
    }
    count = 0;
    int maxcolsize = 0;
    int maxcolfuncs = 0;
    for (int i = 0; i < pfock->npcol; i++) {
        compute_FD_ptr (pfock,
                        pfock->colptr_sh[i], pfock->colptr_sh[i+1] - 1,
                        ptrcol, &(pfock->colsize[i]));
        maxcolsize =
            pfock->colsize[i] > maxcolsize ? pfock->colsize[i] : maxcolsize;
        int nfuncs = pfock->colptr_f[i + 1] - pfock->colptr_f[i];
        maxcolfuncs = nfuncs > maxcolfuncs ? nfuncs : maxcolfuncs;
        if (i == mycol) {
            pfock->sizemycol = pfock->colsize[i];
            init_FD_load (pfock, ptrcol,
                          &(pfock->loadcol), &(pfock->sizeloadcol));  
        }
        for (int j = pfock->colptr_sh[i]; j < pfock->colptr_sh[i+1]; j++) {
            pfock->colpos[j] = ptrcol[j];
            for (int k = pfock->shellptr[j]; k < pfock->shellptr[j+1]; k++) {
                int sh = pfock->shellid[k];
                pfock->colptr[count] = ptrcol[sh];
                count++;
            }
        }
    }
    PFOCK_FREE(ptrrow);
    PFOCK_FREE(ptrcol);
    pfock->maxrowsize = maxrowsize;
    pfock->maxcolsize = maxcolsize;
    pfock->maxrowfuncs = maxrowfuncs;
    pfock->maxcolfuncs = maxcolfuncs;
    int sizeX1 = maxrowfuncs * maxrowsize;
    int sizeX2 = maxcolfuncs * maxcolsize;
    int sizeX3 = maxrowsize * maxcolsize;
    pfock->sizeX1 = sizeX1;
    pfock->sizeX2 = sizeX2;
    pfock->sizeX3 = sizeX3;
    //if (myrank == 0) {
    //    printf("  FD size (%d %d %d %d)\n",
    //        maxrowfuncs, maxcolfuncs, maxrowsize, maxcolsize);
    //}
    
    // D buf
    pfock->D1 = (double **)PFOCK_MALLOC(sizeof(double *) * pfock->max_numdmat2);
    pfock->D2 = (double **)PFOCK_MALLOC(sizeof(double *) * pfock->max_numdmat2); 
    pfock->D3 = (double **)PFOCK_MALLOC(sizeof(double *) * pfock->max_numdmat2);
    if (NULL == pfock->D1 ||
        NULL == pfock->D2 ||
        NULL == pfock->D3) {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    for (int i = 0; i < pfock->max_numdmat2; i++) {
        pfock->D1[i] = (double *)PFOCK_MALLOC(sizeof(double) * sizeX1);
        pfock->D2[i] = (double *)PFOCK_MALLOC(sizeof(double) * sizeX2); 
        pfock->D3[i] = (double *)PFOCK_MALLOC(sizeof(double) * sizeX3);
        pfock->mem_cpu += 1.0 * sizeof(double) * (sizeX1 + sizeX2 + sizeX3);
        if (NULL == pfock->D1[i] ||
            NULL == pfock->D2[i] ||
            NULL == pfock->D3[i]) {
            PFOCK_PRINTF (1, "memory allocation failed\n");
            return PFOCK_STATUS_ALLOC_FAILED;
        }
    }

    
    // F buf
    int nthreads = pfock->nthreads;
    char *ncpu_str = getenv("nCPU_F");
    int ncpu_f;
    if (ncpu_str == NULL) {
        ncpu_f = 1;
    } else {
        ncpu_f = atoi(ncpu_str);
        if (ncpu_f <= 0 || ncpu_f > nthreads) {
            ncpu_f = 1;
        }
    }
    
    int sizeX4 = maxrowfuncs * maxcolfuncs;
    int sizeX6 = maxrowsize * maxcolfuncs;
    int sizeX5 = maxrowfuncs * maxcolsize;
    pfock->sizeX4 = sizeX4;
    pfock->sizeX5 = sizeX5;
    pfock->sizeX6 = sizeX6;
    pfock->ncpu_f = ncpu_f;
    int numF = pfock->numF = (nthreads + ncpu_f - 1)/ncpu_f;
    // allocation

    pfock->F1 = (double *)PFOCK_MALLOC(sizeof(double) * sizeX1 *
        numF * pfock->max_numdmat2);
    pfock->F2 = (double *)PFOCK_MALLOC(sizeof(double) * sizeX2 *
        numF * pfock->max_numdmat2); 
    pfock->F3 = (double *)PFOCK_MALLOC(sizeof(double) * sizeX3 *
        pfock->max_numdmat2);
    pfock->F4 = (double *)PFOCK_MALLOC(sizeof(double) * sizeX4 *
        numF * pfock->max_numdmat2);
    pfock->F5 = (double *)PFOCK_MALLOC(sizeof(double) * sizeX5 *
        numF * pfock->max_numdmat2); 
    pfock->F6 = (double *)PFOCK_MALLOC(sizeof(double) * sizeX6 *
        numF * pfock->max_numdmat2);
    pfock->mem_cpu += 1.0 * sizeof(double) *
        (((double)sizeX1 + sizeX2 + sizeX4 + sizeX5 + sizeX6) * 
        numF + sizeX3) * pfock->max_numdmat2;
    if (NULL == pfock->F1 ||
        NULL == pfock->F2 ||
        NULL == pfock->F3 ||
        NULL == pfock->F4 ||
        NULL == pfock->F5 ||
        NULL == pfock->F6) {
        PFOCK_PRINTF (1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    } 

    pfock->ldX1 = maxrowsize;
    pfock->ldX2 = maxcolsize;
    pfock->ldX3 = maxcolsize;
    pfock->ldX4 = maxcolfuncs;
    pfock->ldX5 = maxcolsize;
    pfock->ldX6 = maxcolfuncs;        
    return PFOCK_STATUS_SUCCESS;
}


static void destroy_buffers (PFock_t pfock)
{
    for (int i = 0; i < pfock->max_numdmat2; i++) {
        GA_Destroy(pfock->ga_D1[i]);
        GA_Destroy(pfock->ga_D2[i]);
        GA_Destroy(pfock->ga_D3[i]);
        GA_Destroy(pfock->ga_F1[i]);
        GA_Destroy(pfock->ga_F2[i]);
        GA_Destroy(pfock->ga_F3[i]);
    }
    PFOCK_FREE(pfock->ga_D1);
    PFOCK_FREE(pfock->ga_D2);
    PFOCK_FREE(pfock->ga_D3);
    PFOCK_FREE(pfock->ga_F1);
    PFOCK_FREE(pfock->ga_F2);
    PFOCK_FREE(pfock->ga_F3);
        
    PFOCK_FREE(pfock->rowpos);
    PFOCK_FREE(pfock->colpos);
    PFOCK_FREE(pfock->rowptr);
    PFOCK_FREE(pfock->colptr);
    PFOCK_FREE(pfock->loadrow);
    PFOCK_FREE(pfock->loadcol);
    PFOCK_FREE(pfock->rowsize);
    PFOCK_FREE(pfock->colsize);

    for (int i = 0; i < pfock->max_numdmat2; i++) {
        PFOCK_FREE(pfock->D1[i]);
        PFOCK_FREE(pfock->D2[i]);
        PFOCK_FREE(pfock->D3[i]);
    }
    PFOCK_FREE(pfock->D1);
    PFOCK_FREE(pfock->D2);
    PFOCK_FREE(pfock->D3);
    PFOCK_FREE(pfock->F1);
    PFOCK_FREE(pfock->F2);
    PFOCK_FREE(pfock->F3);
    PFOCK_FREE(pfock->F4);
    PFOCK_FREE(pfock->F5);
    PFOCK_FREE(pfock->F6);    
}


PFockStatus_t init_GA(int nbf, int nprow, int npcol,
                      int num_dmat, int sizeheap, int sizestack)
{
    int maxrowsize = (nbf + nprow - 1)/nprow;
    int maxcolsize = (nbf + npcol - 1)/npcol;    
    int heap = num_dmat * 5 * maxrowsize * maxcolsize;
    int stack = heap;
    heap += sizeheap;
    stack += sizestack;

    GA_Initialize();
    if (!MA_initialized()) {
        if (!MA_init(C_DBL, heap, stack)) {
            return PFOCK_STATUS_INIT_FAILED;
        }
    }
    
    return PFOCK_STATUS_SUCCESS;
}


void finalize_GA()
{    
    GA_Terminate();
}


PFockStatus_t PFock_create(BasisSet_t basis, int nprow, int npcol, int ntasks,
                           double tolscr, int max_numdmat, int symm,
                           PFock_t *_pfock)
{
    // allocate pfock
    PFock_t pfock = (PFock_t)PFOCK_MALLOC(sizeof(struct PFock));    
    if (NULL == pfock) {
        PFOCK_PRINTF(1, "Failed to allocate memory: %lld\n",
            sizeof(struct PFock));
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    memset(pfock, 0, sizeof(PFock_t));
    
    // check if MPI is initialized
    int flag;    
    MPI_Initialized(&flag);
    if (!flag) {
        PFOCK_PRINTF(1, "MPI_Init() or MPI_Init_thread()"
                     " has not been called\n");
        return PFOCK_STATUS_INIT_FAILED;        
    }
    int nprocs;
    int myrank;   
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);         
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // initialization
    pfock->nosymm = (symm == 0 ? 1 : 0);
    pfock->maxnfuncs = CInt_getMaxShellDim (basis);
    pfock->nbf = CInt_getNumFuncs (basis);
    pfock->nshells = CInt_getNumShells (basis);
    pfock->natoms = CInt_getNumAtoms (basis);
    pfock->nthreads = omp_get_max_threads ();
    pfock->mem_cpu = 0.0;
    omp_set_num_threads (pfock->nthreads);
    
    // check inputs
    if (nprow <= 0 || nprow > pfock->nshells ||
        npcol <= 0 || npcol > pfock->nshells ||
        (nprow * npcol) > nprocs) {
        PFOCK_PRINTF(1, "Invalid nprow or npcol\n");
        return PFOCK_STATUS_INVALID_VALUE;
    } else {
        pfock->nprow= nprow;
        pfock->npcol = npcol;
        pfock->nprocs = nprow * npcol;
    }
    if (tolscr < 0.0) {
        PFOCK_PRINTF(1, "Invalid screening threshold\n");
        return PFOCK_STATUS_INVALID_VALUE;
    } else {
        pfock->tolscr = tolscr;
        pfock->tolscr2 = tolscr * tolscr;
    }
    if (max_numdmat <= 0) {
        PFOCK_PRINTF(1, "Invalid number of density matrices\n");
        return PFOCK_STATUS_INVALID_VALUE;
    } else {
        pfock->max_numdmat = max_numdmat;
        // Ewww. This is max_numdmat if nosymm is 0 and 2 * max_numdmat otherwise.
        pfock->max_numdmat2 = (pfock->nosymm + 1) * max_numdmat;
    }
          
    // init global arrays
    init_GA(pfock->nbf, nprow, npcol, pfock->max_numdmat2, 0, 0);
        
    // set tasks
    int minnshells = (nprow > npcol ? nprow : npcol);
    minnshells = pfock->nshells/minnshells;
    if (ntasks >= minnshells) {
        pfock->nbp_p = minnshells;
    } else if (ntasks <= 0) {
        pfock->nbp_p = 4;
        pfock->nbp_p = MIN (pfock->nbp_p, minnshells);
    } else {
        pfock->nbp_p = ntasks;
    }
    pfock->nbp_row = pfock->nbp_col = pfock->nbp_p;
       
    // functions starting positions of shells
    pfock->f_startind =
        (int *)PFOCK_MALLOC(sizeof(int) * (pfock->nshells + 1));
    pfock->mem_cpu += sizeof(int) * (pfock->nshells + 1);   
    if (NULL == pfock->f_startind) {
        PFOCK_PRINTF(1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }  
    for (int i = 0; i < pfock->nshells; i++) {
        pfock->f_startind[i] = CInt_getFuncStartInd(basis, i);
    }
    pfock->f_startind[pfock->nshells] = pfock->nbf;

    // shells starting positions of atoms
    pfock->s_startind =
        (int *)PFOCK_MALLOC(sizeof(int) * (pfock->natoms + 1));
    pfock->mem_cpu += sizeof(int) * (pfock->natoms + 1); 
    if (NULL == pfock->s_startind) {
        PFOCK_PRINTF(1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    for (int i = 0; i < pfock->natoms; i++) {
        pfock->s_startind[i] = CInt_getAtomStartInd(basis, i);
    }
    pfock->s_startind[pfock->natoms] = pfock->nshells;


    PFockStatus_t ret;
    // init comm
    if ((ret = init_fock(pfock)) != PFOCK_STATUS_SUCCESS) {
        return ret;
    }

    // init scheduler
    if (init_taskq(pfock) != 0) {
        PFOCK_PRINTF(1, "task queue initialization failed\n");
        return PFOCK_STATUS_INIT_FAILED;
    }

    // schwartz screening    
    if (myrank == 0) {
        //PFOCK_INFO("screening ...\n");
    }
    double t1 = MPI_Wtime();
    if (schwartz_screening(pfock, basis) != 0) {
        PFOCK_PRINTF (1, "schwartz screening failed\n");
        return PFOCK_STATUS_INIT_FAILED;
    }
    double t2 = MPI_Wtime();
    if (myrank == 0) {
        //PFOCK_INFO("takes %.3lf secs\n", t2 - t1);
    }

    // repartition
    if ((ret = repartition_fock(pfock)) != PFOCK_STATUS_SUCCESS) {
        return ret;
    }

    // init global arrays
    if ((ret = create_GA(pfock)) != PFOCK_STATUS_SUCCESS) {
        return ret;
    }

    // create local buffers
    if ((ret = create_buffers(pfock)) != PFOCK_STATUS_SUCCESS) {
        return ret;
    }

    if ((ret = create_FD_GArrays(pfock)) != PFOCK_STATUS_SUCCESS) {
        return ret;
    }

    //CInt_createERD(basis, &(pfock->erd), pfock->nthreads);

    // statistics
    pfock->mpi_timepass
        = (double *)PFOCK_MALLOC(sizeof(double) * pfock->nprocs);
    pfock->mpi_timereduce
        = (double *)PFOCK_MALLOC(sizeof(double) * pfock->nprocs);
    pfock->mpi_timeinit
        = (double *)PFOCK_MALLOC(sizeof(double) * pfock->nprocs);
    pfock->mpi_timecomp
        = (double *)PFOCK_MALLOC(sizeof(double) * pfock->nprocs);
    pfock->mpi_timegather
        = (double *)PFOCK_MALLOC(sizeof(double) * pfock->nprocs);
    pfock->mpi_timescatter
        = (double *)PFOCK_MALLOC(sizeof(double) * pfock->nprocs);
    pfock->mpi_usq
        = (double *)PFOCK_MALLOC(sizeof(double) * pfock->nprocs);
    pfock->mpi_uitl
        = (double *)PFOCK_MALLOC(sizeof(double) * pfock->nprocs);
    pfock->mpi_steals
        = (double *)PFOCK_MALLOC(sizeof(double) * pfock->nprocs);
    pfock->mpi_stealfrom
        = (double *)PFOCK_MALLOC(sizeof(double) * pfock->nprocs);
    pfock->mpi_ngacalls
        = (double *)PFOCK_MALLOC(sizeof(double) * pfock->nprocs);
    pfock->mpi_volumega
        = (double *)PFOCK_MALLOC(sizeof(double) * pfock->nprocs);
    if (pfock->mpi_timepass == NULL ||
        pfock->mpi_timereduce == NULL ||
        pfock->mpi_timeinit == NULL ||
        pfock->mpi_timecomp == NULL ||        
        pfock->mpi_usq == NULL ||
        pfock->mpi_uitl == NULL ||
        pfock->mpi_steals == NULL ||
        pfock->mpi_stealfrom == NULL ||
        pfock->mpi_ngacalls == NULL ||
        pfock->mpi_volumega == NULL ||
        pfock->mpi_timegather == NULL ||
        pfock->mpi_timescatter == NULL) {
        PFOCK_PRINTF(1, "memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    
    pfock->committed = 0;
    *_pfock = pfock;
    
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_destroy(PFock_t pfock)
{
    PFOCK_FREE(pfock->blkrowptr_sh);
    PFOCK_FREE(pfock->blkcolptr_sh);
    PFOCK_FREE(pfock->rowptr_sh);
    PFOCK_FREE(pfock->colptr_sh);
    PFOCK_FREE(pfock->rowptr_f);
    PFOCK_FREE(pfock->colptr_f);
    PFOCK_FREE(pfock->rowptr_blk);
    PFOCK_FREE(pfock->colptr_blk);
    PFOCK_FREE(pfock->FT_block);
    PFOCK_FREE(pfock->f_startind);
    PFOCK_FREE(pfock->s_startind);

    //CInt_destroyERD(pfock->erd);
    clean_taskq(pfock);
    clean_screening(pfock);
    destroy_GA(pfock);
    destroy_buffers(pfock);
    // I thought we should call that but apparently this
    // breaks everything. -JFG 
    //finalize_GA();

    PFOCK_FREE(pfock->mpi_timepass);
    PFOCK_FREE(pfock->mpi_timereduce);
    PFOCK_FREE(pfock->mpi_timeinit);
    PFOCK_FREE(pfock->mpi_timecomp);
    PFOCK_FREE(pfock->mpi_timegather);
    PFOCK_FREE(pfock->mpi_timescatter);
    PFOCK_FREE(pfock->mpi_usq);
    PFOCK_FREE(pfock->mpi_uitl);
    PFOCK_FREE(pfock->mpi_steals);
    PFOCK_FREE(pfock->mpi_stealfrom);
    PFOCK_FREE(pfock->mpi_ngacalls);
    PFOCK_FREE(pfock->mpi_volumega);
    
    PFOCK_FREE(pfock);
  
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_setNumDenMat(int numdmat, PFock_t pfock)
{
    if (pfock->committed == 1) {
        PFOCK_PRINTF(1, "Can't change number of matrices"
                      " after PFock_commitDenMats() is called.\n");
        return PFOCK_STATUS_EXECUTION_FAILED;
    }
    if (numdmat <= 0 || numdmat > pfock->max_numdmat) {
        PFOCK_PRINTF(1, "Invalid number of density matrices\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
    int numdmat2 = numdmat * (pfock->nosymm + 1);
    pfock->num_dmat = numdmat;
    pfock->num_dmat2 = numdmat2;
    
    return PFOCK_STATUS_EXECUTION_FAILED;
}


PFockStatus_t PFock_putDenMat(int rowstart, int rowend,
                              int colstart, int colend,
                              int stride, double *dmat,
                              int index, PFock_t pfock)
{
    int lo[2];
    int hi[2];
    int ld[1];

    if (pfock->committed == 1) {
        PFOCK_PRINTF (1, "Can't change density matrix"
                      " after PFock_commitDenMats() is called.\n");
        return PFOCK_STATUS_EXECUTION_FAILED;
    }
    if (index < 0 || index >= pfock->num_dmat) {
        PFOCK_PRINTF (1, "Invalid index\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
    
    lo[0] = rowstart;
    hi[0] = rowend;    
    lo[1] = colstart;
    hi[1] = colend;
    ld[0] = stride;
    
    NGA_Put(pfock->ga_D[index], lo, hi, (void *)dmat, ld);
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_putDenMatGA(int ga, int index, PFock_t pfock)
{
    GA_Copy(ga, pfock->ga_D[index]);
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_fillDenMat(double value, int index,
                               PFock_t pfock)
{
    if (pfock->committed == 1) {
        PFOCK_PRINTF (1, "Can't change density matrix"
                      " after PFock_commitDenMats() is called.\n");
        return PFOCK_STATUS_EXECUTION_FAILED;
    }
    if (index < 0 || index >= pfock->num_dmat) {
        PFOCK_PRINTF (1, "Invalid index\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
    
    GA_Fill(pfock->ga_D[index], &value);
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_commitDenMats(PFock_t pfock)
{
    GA_Sync();
    pfock->committed = 1;
    if (pfock->nosymm == 1) {
        for (int i = 0; i < pfock->num_dmat; i++) {
            GA_Transpose(pfock->ga_D[i], pfock->ga_D[i + pfock->num_dmat]);
        }
    }

    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_sync(PFock_t pfock)
{
    GA_Sync();
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_getMat(PFock_t pfock, PFockMatType_t type,
                           int index,
                           int rowstart, int rowend,
                           int colstart, int colend,
                           int stride, double *mat)
{
    int lo[2];
    int hi[2];
    int ld[1];
    int *ga;
    if (index < 0 || index >= pfock->max_numdmat) {
        PFOCK_PRINTF (1, "Invalid index\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
#ifdef _SCF_
    if (PFOCK_MAT_TYPE_J == type || PFOCK_MAT_TYPE_K == type) {
        PFOCK_PRINTF (1, "Invalid matrix type\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
#endif    

    lo[0] = rowstart;
    hi[0] = rowend;    
    lo[1] = colstart;
    hi[1] = colend;
    ld[0] = stride;
    
    ga = pfock->gatable[type];
    NGA_Get(ga[index], lo, hi, mat, ld);

#ifndef __SCF__
    if (PFOCK_MAT_TYPE_F == type) {
        int sizerow = rowend - rowstart + 1;
        int sizecol = colend - colstart + 1;
        double *K = (double *)PFOCK_MALLOC(sizerow * sizecol * sizeof(double));
        if (NULL == K) {
            PFOCK_PRINTF(1, "Failed to allocate memory: %lld\n",
                sizerow * sizecol * sizeof(double));
            return PFOCK_STATUS_ALLOC_FAILED;
        }
        int ga_K = pfock->ga_K[index];
        NGA_Get(ga_K, lo, hi, K, &stride);
        #pragma omp parallel for
        for (int i = 0; i < sizerow; i++) {
            #pragma simd
            for (int j = 0; j < sizecol; j++) {
                mat[i * stride + j] += K[i * sizecol + j];
            }
        }
        PFOCK_FREE(K);
    }    
#endif

    return PFOCK_STATUS_SUCCESS;    
}


PFockStatus_t PFock_getMatGA(PFock_t pfock, PFockMatType_t type,
                             int index, int ga)
{
    int *my_ga = pfock->gatable[type];
    GA_Copy(my_ga[index], ga);
    
 #ifndef __SCF__
    if (PFOCK_MAT_TYPE_F == type) {
        int ga_K = pfock->ga_K[index];
        double fone = 1.0;
        double fzero = 0.0;
        GA_Add(&fone, ga_K, &fzero, ga, ga);
    }    
#endif

    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_getLocalMatInds(PFock_t pfock,
                                    int *rowstart, int *rowend,
                                    int *colstart, int *colend)
{
    int lo[2];
    int hi[2];
    int myrank;  
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    NGA_Distribution(pfock->ga_D[0], myrank, lo, hi);
    *rowstart = lo[0];
    *rowend = hi[0];
    *colstart = lo[1];
    *colend = hi[1];

    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_getLocalMatPtr(PFock_t pfock,
                                   PFockMatType_t type, int index,
                                   int *rowstart, int *rowend,
                                   int *colstart, int *colend,
                                   int *stride, double **mat)
{
    int lo[2];
    int hi[2];
    int myrank;
    int *ga;

    if (index < 0 || index >= pfock->max_numdmat)
    {
        PFOCK_PRINTF (1, "Invalid index\n");
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


PFockStatus_t PFock_getMatGAHandle(PFock_t pfock,
                                   PFockMatType_t type, int index,
                                   int *ga)
{
    if (index < 0 || index >= pfock->max_numdmat)
    {
        PFOCK_PRINTF (1, "Invalid index\n");
        return PFOCK_STATUS_INVALID_VALUE;
    }
    
    int *g = pfock->gatable[type];
    *ga = g[index];
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_computeFock(BasisSet_t basis,
                                PFock_t pfock)
{
#ifdef GA_NB
    ga_nbhdl_t nbhdlF1;
    ga_nbhdl_t nbhdlF2;
    ga_nbhdl_t nbhdlF3;
#endif
    struct timeval tv1;
    struct timeval tv2;
    struct timeval tv3;
    struct timeval tv4; 
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    pfock->committed = 0;      
    pfock->timepass = 0.0;
    pfock->timereduce = 0.0;
    pfock->timeinit = 0.0;
    pfock->timecomp = 0.0;
    pfock->timegather = 0.0;
    pfock->timescatter = 0.0;
    pfock->usq = 0.0;
    pfock->uitl = 0.0;
    pfock->steals = 0.0;
    pfock->stealfrom = 0.0;
    pfock->ngacalls = 0.0;
    pfock->volumega = 0.0;
    int my_sshellrow = pfock->sshell_row;
    int my_sshellcol = pfock->sshell_col;
    int myrow = myrank/pfock->npcol;
    int mycol = myrank%pfock->npcol;
    int sizeX1 = pfock->sizeX1;
    int sizeX2 = pfock->sizeX2;
    int sizeX3 = pfock->sizeX3;
    int sizeX4 = pfock->sizeX4;    
    int sizeX5 = pfock->sizeX5;
    int sizeX6 = pfock->sizeX6;
    double *D1[pfock->num_dmat2];
    double *D2[pfock->num_dmat2];
    double *D3[pfock->num_dmat2];
    double *F1 = pfock->F1;
    double *F2 = pfock->F2;
    double *F3 = pfock->F3;
    double *F4 = pfock->F4;
    double *F5 = pfock->F5;
    double *F6 = pfock->F6;
    int maxrowsize = pfock->maxrowsize;
    int maxcolfuncs = pfock->maxcolfuncs;
    int maxcolsize = pfock->maxcolsize;
    int ldX1 = maxrowsize;
    int ldX2 = maxcolsize;
    int ldX3 = maxcolsize;
    int ldX4 = maxcolfuncs;
    int ldX5 = maxcolsize;
    int ldX6 = maxcolfuncs;
    double dzero = 0.0;
    double done = 1.0;
    int lo[2];
    int hi[2];

    gettimeofday (&tv1, NULL);    
    gettimeofday (&tv3, NULL);
    for (int i = 0; i < pfock->num_dmat2; i++) {
        GA_Fill(pfock->ga_F[i], &dzero);
    #ifndef __SCF__
        GA_Fill(pfock->ga_K[i], &dzero);
    #endif
        GA_Fill(pfock->ga_F1[i], &dzero);
        GA_Fill(pfock->ga_F2[i], &dzero);
        GA_Fill(pfock->ga_F3[i], &dzero);
    }
    // local my D
    load_local_bufD(pfock);
    lo[0] = myrank;
    hi[0] = myrank;
    lo[1] = 0;
    for (int i = 0; i < pfock->num_dmat2; i++) {
        int ldD;
        hi[1] = sizeX1 - 1;
        NGA_Access(pfock->ga_D1[i], lo, hi, &D1[i], &ldD);
        hi[1] = sizeX2 - 1;
        NGA_Access(pfock->ga_D2[i], lo, hi, &D2[i], &ldD);
        hi[1] = sizeX3 - 1;
        NGA_Access(pfock->ga_D3[i], lo, hi, &D3[i], &ldD);
    }
    gettimeofday(&tv4, NULL);
    pfock->timegather += (tv4.tv_sec - tv3.tv_sec) +
        (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
    pfock->ngacalls += 3;
    pfock->volumega += (sizeX1 + sizeX2 + sizeX3) * sizeof(double);
    
    gettimeofday (&tv3, NULL);   
    reset_F(pfock->numF, pfock->num_dmat2, F1, F2, F3, F4, F5, F6,
            sizeX1, sizeX2, sizeX3, sizeX4, sizeX5, sizeX6);
    gettimeofday (&tv4, NULL);
    pfock->timeinit += (tv4.tv_sec - tv3.tv_sec) +
        (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
    
    /* own part */
    reset_taskq(pfock);
    int task;
    while ((task = taskq_next (pfock, myrow, mycol, 1)) < pfock->ntasks) {
        int rowid = task/pfock->nblks_col;
        int colid = task%pfock->nblks_col;
        int startM = pfock->blkrowptr_sh[pfock->sblk_row + rowid];
        int endM = pfock->blkrowptr_sh[pfock->sblk_row + rowid + 1] - 1;
        int startP = pfock->blkcolptr_sh[pfock->sblk_col + colid];
        int endP = pfock->blkcolptr_sh[pfock->sblk_col + colid + 1] - 1;
        gettimeofday (&tv3, NULL);       
        //fock_task(basis, pfock->erd, pfock->ncpu_f, pfock->num_dmat2,
        fock_task(basis, pfock->ncpu_f, pfock->num_dmat2,
                  pfock->shellptr, pfock->shellvalue,
                  pfock->shellid, pfock->shellrid,
                  pfock->f_startind,
                  pfock->rowpos, pfock->colpos,
                  pfock->rowptr, pfock->colptr,
                  pfock->tolscr2,
                  my_sshellrow, my_sshellcol,
                  startM, endM, startP, endP,
                  D1, D2, D3, F1, F2, F3, F4, F5, F6,
                  ldX1, ldX2, ldX3, ldX4, ldX5, ldX6,
                  sizeX1, sizeX2, sizeX3,
                  sizeX4, sizeX5, sizeX6,
                  &(pfock->uitl), &(pfock->usq),pfock->nosymm);
        gettimeofday (&tv4, NULL);
        pfock->timecomp += (tv4.tv_sec - tv3.tv_sec) +
                    (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
    } /* own part */

    gettimeofday (&tv3, NULL);        
    // reduction on CPU
    reduce_F(pfock->numF, pfock->num_dmat2, F1, F2, F3, F4, F5, F6,
             sizeX1, sizeX2, sizeX3,
             sizeX4, sizeX5, sizeX6,
             maxrowsize, maxcolsize,
             pfock->nfuncs_row, pfock->nfuncs_col,
             pfock->rowpos[my_sshellrow],
             pfock->colpos[my_sshellcol],
             ldX3, ldX4, ldX5, ldX6);   
    lo[0] = myrank;
    hi[0] = myrank;
    lo[1] = 0;    
    for (int i = 0; i < pfock->num_dmat2; i++) {
        hi[1] = sizeX1 - 1;
#ifdef GA_NB
        // save results for local intergrals
        NGA_NbAcc(pfock->ga_F1[i], lo, hi, &F1[i * sizeX1],
                  &sizeX1, &done, &nbhdlF1);
        hi[1] = sizeX2 - 1;
        NGA_NbAcc(pfock->ga_F2[i], lo, hi, &F2[i * sizeX2],
                  &sizeX2, &done, &nbhdlF2); 
        hi[1] = sizeX3 - 1;
        NGA_NbAcc(pfock->ga_F3[i], lo, hi, &F3[i * sizeX3],
                  &sizeX3, &done, &nbhdlF3); 
#else
        // save results for local intergrals
        NGA_Acc(pfock->ga_F1[i], lo, hi, &F1[i * sizeX1], &sizeX1, &done);
        hi[1] = sizeX2 - 1;
        NGA_Acc(pfock->ga_F2[i], lo, hi, &F2[i * sizeX2], &sizeX2, &done); 
        hi[1] = sizeX3 - 1;
        NGA_Acc(pfock->ga_F3[i], lo, hi, &F3[i * sizeX3], &sizeX3, &done);
#endif
    }
    pfock->ngacalls += 3;
    pfock->volumega += (sizeX1 + sizeX2 + sizeX3) * sizeof(double);
    gettimeofday (&tv4, NULL);
    pfock->timereduce += (tv4.tv_sec - tv3.tv_sec) +
        (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;

    // to steal
    pfock->steals = 0;
    pfock->stealfrom = 0;
#ifdef __DYNAMIC__
#ifdef GA_NB
    ga_nbhdl_t nbhdlD1;
    ga_nbhdl_t nbhdlD2;
    ga_nbhdl_t nbhdlD3;
#endif
    double **D1_task;
    double **D2_task;
    double **VD1 = pfock->D1;
    double **VD2 = pfock->D2;
    double **VD3 = pfock->D3;
    int prevrow = myrow;
    int prevcol = mycol;   
    /* steal tasks */
    for (int idx = 0; idx < pfock->nprocs - 1; idx++) {
        int vpid = (myrank + idx + 1)%pfock->nprocs;
        int vrow = vpid/pfock->npcol;
        int vcol = vpid%pfock->npcol;
        int vsblk_row = pfock->rowptr_blk[vrow];
        int vsblk_col = pfock->colptr_blk[vcol];
        int vnblks_col = pfock->colptr_blk[vcol + 1] - vsblk_col;       
        int vnfuncs_row = pfock->rowptr_f[vrow + 1] - pfock->rowptr_f[vrow];
        int vnfuncs_col = pfock->colptr_f[vcol + 1] - pfock->colptr_f[vcol];
        int vsshellrow = pfock->rowptr_sh[vrow];   
        int vsshellcol = pfock->colptr_sh[vcol];
        int stealed = 0;
        int task;
        while ((task = taskq_next(pfock, vrow, vcol, 1)) < pfock->ntasks) {
            gettimeofday (&tv3, NULL);
            if (0 == stealed) {
                if (vrow != prevrow && vrow != myrow) {
                    D1_task = VD1;
                    lo[0] = vpid;
                    hi[0] = vpid;
                    lo[1] = 0;
                    hi[1] = sizeX1 - 1;
                    for (int i = 0; i < pfock->num_dmat2; i++) {
                    #ifdef GA_NB    
                        NGA_NbGet(pfock->ga_D1[i], lo, hi,
                                  VD1[i], &sizeX1, &nbhdlD1);
                    #else
                        NGA_Get(pfock->ga_D1[i], lo, hi,  VD1[i], &sizeX1);
                    #endif
                        pfock->ngacalls += 1;
                        pfock->volumega += sizeX1 * sizeof(double);
                    }                 
                } else if (vrow == myrow) {
                    D1_task = D1;
                }  
                if (vcol != prevcol && vcol != mycol) {
                    D2_task =  VD2;
                    lo[0] = vpid;
                    hi[0] = vpid;
                    lo[1] = 0;
                    hi[1] = sizeX2 - 1;
                    for (int i = 0; i < pfock->num_dmat2; i++) {
                    #ifdef GA_NB
                        NGA_NbGet(pfock->ga_D2[i], lo, hi,
                                  VD2[i], &sizeX2, &nbhdlD2);
                    #else
                        NGA_Get(pfock->ga_D2[i], lo, hi, VD2[i], &sizeX2);               
                    #endif
                        pfock->ngacalls += 1;
                        pfock->volumega += sizeX2 * sizeof(double);
                    }                  
                } else if (vcol == mycol) {
                    D2_task = D2;
                }
                lo[0] = vpid;
                hi[0] = vpid;
                lo[1] = 0;
                hi[1] = sizeX3 - 1;
                for (int i = 0; i < pfock->num_dmat2; i++) {
                #ifdef GA_NB
                    NGA_NbGet(pfock->ga_D3[i], lo, hi, 
                              VD3[i], &sizeX3, &nbhdlD3);
                #else
                    NGA_Get(pfock->ga_D3[i], lo, hi,  VD3[i], &sizeX3);
                #endif
                    pfock->ngacalls += 1;
                    pfock->volumega += sizeX3 * sizeof(double);
                }
            #ifdef GA_NB    
                // wait for last NbAcc F
                NGA_NbWait(&nbhdlF1);
                NGA_NbWait(&nbhdlF2);
                NGA_NbWait(&nbhdlF3);
            #endif    
                // init F bufs
                reset_F(pfock->numF, pfock->num_dmat2, F1, F2, F3, F4, F5, F6,
                        sizeX1, sizeX2, sizeX3, sizeX4, sizeX5, sizeX6);
            #ifdef GA_NB    
                // wait for NbGet
                if (vrow != prevrow && vrow != myrow) {
                    NGA_NbWait(&nbhdlD1);
                }
                if (vcol != prevcol && vcol != mycol) {
                    NGA_NbWait(&nbhdlD2);
                }
                NGA_NbWait(&nbhdlD3);
            #endif    
                pfock->stealfrom++;
            }
            gettimeofday (&tv4, NULL);
            pfock->timeinit += (tv4.tv_sec - tv3.tv_sec) +
                   (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
            int rowid = task/vnblks_col;
            int colid = task%vnblks_col;
            // compute task
            int startM = pfock->blkrowptr_sh[vsblk_row + rowid];
            int endM = pfock->blkrowptr_sh[vsblk_row + rowid + 1] - 1;
            int startP = pfock->blkcolptr_sh[vsblk_col + colid];
            int endP = pfock->blkcolptr_sh[vsblk_col + colid + 1] - 1;

            gettimeofday (&tv3, NULL);
            // Should we modify this as well?
            //fock_task(basis, pfock->erd, pfock->ncpu_f, pfock->num_dmat2,
            fock_task(basis, pfock->ncpu_f, pfock->num_dmat2,
                      pfock->shellptr, pfock->shellvalue,
                      pfock->shellid, pfock->shellrid,
                      pfock->f_startind,
                      pfock->rowpos, pfock->colpos,
                      pfock->rowptr, pfock->colptr,
                      pfock->tolscr2,
                      vsshellrow, vsshellcol, startM, endM, startP, endP,
                      D1_task, D2_task, VD3, F1, F2, F3, F4, F5, F6,
                      ldX1, ldX2, ldX3, ldX4, ldX5, ldX6,
                      sizeX1, sizeX2, sizeX3, sizeX4, sizeX5, sizeX6,
                      &(pfock->uitl), &(pfock->usq), pfock->nosymm);
            gettimeofday (&tv4, NULL);
            pfock->timecomp += (tv4.tv_sec - tv3.tv_sec) +
                        (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
            pfock->steals++;
            stealed = 1;
        }
        gettimeofday (&tv3, NULL);
        if (1 == stealed) {
            // reduction
            reduce_F(pfock->numF, pfock->num_dmat2, F1, F2, F3, F4, F5, F6,
                     sizeX1, sizeX2, sizeX3,
                     sizeX4, sizeX5, sizeX6,
                     maxrowsize, maxcolsize,
                     vnfuncs_row, vnfuncs_col,
                     pfock->rowpos[vsshellrow],
                     pfock->colpos[vsshellcol],
                     ldX3, ldX4, ldX5, ldX6);
            lo[1] = 0;
            hi[1] = sizeX1 - 1;
            if (vrow != myrow) {
                lo[0] = vpid;
                hi[0] = vpid;
                for (int i = 0; i < pfock->num_dmat2; i++) {
                #ifdef GA_NB
                    NGA_NbAcc(pfock->ga_F1[i], lo, hi,
                              &F1[i * sizeX1], &sizeX1, &done, &nbhdlF1);
                #else
                    NGA_Acc(pfock->ga_F1[i], lo, hi,
                            &F1[i * sizeX1], &sizeX1, &done);
                #endif
                    pfock->ngacalls += 1;
                    pfock->volumega += sizeX1 * sizeof(double);
                }
            } else {
                lo[0] = myrank;
                hi[0] = myrank;
                for (int i = 0; i < pfock->num_dmat2; i++) {
                #ifdef GA_NB    
                    NGA_NbAcc(pfock->ga_F1[i], lo, hi,
                              &F1[i * sizeX1], &sizeX1, &done, &nbhdlF1);
                #else
                    NGA_Acc(pfock->ga_F1[i], lo, hi,
                            &F1[i * sizeX1], &sizeX1, &done);
                #endif             
                }
            }
            lo[1] = 0;
            hi[1] = sizeX2 - 1;
            if (vcol != mycol) {
                lo[0] = vpid;
                hi[0] = vpid;
                for (int i = 0; i < pfock->num_dmat2; i++) {
                #ifdef GA_NB    
                    NGA_NbAcc(pfock->ga_F2[i], lo, hi,
                              &F2[i * sizeX2], &sizeX2, &done, &nbhdlF2);
                #else
                    NGA_Acc(pfock->ga_F2[i], lo, hi,
                            &F2[i * sizeX2], &sizeX2, &done);
                #endif
                    pfock->ngacalls += 1;
                    pfock->volumega += sizeX2 * sizeof(double);
                }
            } else {
                lo[0] = myrank;
                hi[0] = myrank;
                for (int i = 0; i < pfock->num_dmat2; i++) {
                #ifdef GA_NB
                    NGA_NbAcc(pfock->ga_F2[i], lo, hi,
                              &F2[i * sizeX2], &sizeX2, &done, &nbhdlF2);
                #else
                    NGA_Acc(pfock->ga_F2[i], lo, hi,
                            &F2[i * sizeX2], &sizeX2, &done);
                #endif
                }
            }
            lo[0] = vpid;
            hi[0] = vpid;
            lo[1] = 0;
            hi[1] = sizeX3 - 1;
            for (int i = 0; i < pfock->num_dmat2; i++) {
            #ifdef GA_NB
                NGA_NbAcc(pfock->ga_F3[i], lo, hi,
                          &F3[i * sizeX3], &sizeX3, &done, &nbhdlF3);
            #else
                NGA_Acc(pfock->ga_F3[i], lo, hi,
                        &F3[i * sizeX3], &sizeX3, &done);
            #endif
                pfock->ngacalls += 1;
                pfock->volumega += sizeX3 * sizeof(double);
            }
            prevrow = vrow;
            prevcol = vcol;
        }
        gettimeofday (&tv4, NULL);
        pfock->timereduce += (tv4.tv_sec - tv3.tv_sec) +
                        (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
    } /* steal tasks */    
#endif /* #ifdef __DYNAMIC__ */

#ifdef GA_NB
    // wait for last NbAcc F
    NGA_NbWait (&nbhdlF1);
    NGA_NbWait (&nbhdlF2);
    NGA_NbWait (&nbhdlF3);
#endif
    lo[0] = myrank;
    hi[0] = myrank;
    lo[1] = 0;
    for (int i = 0; i < pfock->num_dmat2; i++) {
        hi[1] = sizeX1 - 1;
        NGA_Release(pfock->ga_D1[i], lo, hi);
        hi[1] = sizeX2 - 1;
        NGA_Release(pfock->ga_D2[i], lo, hi);
        hi[1] = sizeX3 - 1;
        NGA_Release(pfock->ga_D3[i], lo, hi);
    }
    
    gettimeofday (&tv2, NULL);
    pfock->timepass = (tv2.tv_sec - tv1.tv_sec) +
               (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;    
    GA_Sync();
    
    gettimeofday (&tv3, NULL);
    store_local_bufF (pfock);
    gettimeofday (&tv4, NULL);
    pfock->timescatter = (tv4.tv_sec - tv3.tv_sec) +
               (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;

    if (myrank == 0) {
        //PFOCK_INFO ("correct F ...\n");
    }
  
    if (pfock->nosymm) {
        double dhalf = 0.5;
        for (int i = 0; i < pfock->num_dmat; i++) {
            GA_Transpose(pfock->ga_F[i + pfock->num_dmat],
                         pfock->ga_D[0]);
            GA_Add(&dhalf, pfock->ga_F[i],
                   &dhalf, pfock->ga_D[0], pfock->ga_F[i]);
        #ifndef __SCF__
            GA_Transpose(pfock->ga_K[i + pfock->num_dmat],
                         pfock->ga_D[0]);
            GA_Add(&dhalf, pfock->ga_K[i],
                   &dhalf, pfock->ga_D[0], pfock->ga_K[i]);
        #endif
        }
    } else {
        // correct F
        for (int i = 0; i < pfock->num_dmat; i++) {
            GA_Symmetrize(pfock->ga_F[i]);
        #ifndef __SCF__
            GA_Symmetrize(pfock->ga_K[i]);
        #endif
        }
    }
    
    return PFOCK_STATUS_SUCCESS;
}


/*PFockStatus_t PFock_createCoreHMat(PFock_t pfock, BasisSet_t basis)
{
    int lo[2];
    int hi[2];    
    double *mat;
    int stride;
    double dzero = 0.0;

    int myrank;
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);    
    pfock->ga_H = GA_Duplicate(pfock->ga_D[0], "core hamilton mat");
    if (0 == pfock->ga_H) {
        PFOCK_PRINTF (1, "GA allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    GA_Fill(pfock->ga_H, &dzero);
    NGA_Distribution(pfock->ga_H, myrank, lo, hi);
    NGA_Access(pfock->ga_H, lo, hi, &mat, &stride);    
    compute_H(pfock, basis, pfock->sshell_row, pfock->eshell_row,
              pfock->sshell_col, pfock->eshell_col, stride, mat);
    NGA_Release_update(pfock->ga_H, lo, hi);
    GA_Sync();
    
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_destroyCoreHMat(PFock_t pfock)
{
    GA_Destroy(pfock->ga_H);
    return PFOCK_STATUS_SUCCESS;    
}


PFockStatus_t PFock_getCoreHMat(PFock_t pfock, int rowstart, int rowend,
                                int colstart, int colend,
                                int stride, double *mat)
{
    int lo[2];
    int hi[2];    
    lo[0] = rowstart;
    hi[0] = rowend;    
    lo[1] = colstart;
    hi[1] = colend;    
    NGA_Get(pfock->ga_H, lo, hi, mat, &stride);

    return PFOCK_STATUS_SUCCESS;    
}


PFockStatus_t PFock_createOvlMat(PFock_t pfock, BasisSet_t basis)
{    
    int lo[2];
    int hi[2];
    double *mat;
    int stride;
    double dzero = 0.0;

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    
    pfock->ga_S = GA_Duplicate(pfock->ga_D[0], "overlap mat");
    if (0 == pfock->ga_S) {
        PFOCK_PRINTF (1, "GA allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    // compute S
    GA_Fill(pfock->ga_S, &dzero);
    NGA_Distribution(pfock->ga_S, myrank, lo, hi);
    NGA_Access(pfock->ga_S, lo, hi, &mat, &stride);    
    //compute_S(pfock, basis, pfock->sshell_row, pfock->eshell_row,
    //          pfock->sshell_col, pfock->eshell_col, stride, mat);
    NGA_Release_update(pfock->ga_S, lo, hi);
    GA_Sync();
    
    // compute X         
    int nbf = CInt_getNumFuncs(basis);
    double *eval = (double *)malloc(nbf * sizeof (double));
    if (NULL == eval) {
        PFOCK_PRINTF (1, "Memory allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;        
    }
    int ga_tmp = GA_Duplicate(pfock->ga_D[0], "tmp mat");
    if (0 == ga_tmp) {
        PFOCK_PRINTF (1, "GA allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    int ga_tmp2 = GA_Duplicate(pfock->ga_D[0], "tmp mat");
    if (0 == ga_tmp2) {
        PFOCK_PRINTF (1, "GA allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    //my_peig(pfock->ga_S, ga_tmp, nbf, pfock->nprow, pfock->npcol, eval);
    NGA_Distribution(ga_tmp, myrank, lo, hi);
    int nfuncs_row = hi[0] - lo[0] + 1;
    int nfuncs_col = hi[1] - lo[1] + 1;
    double *blocktmp;
    double *blockS;
    int ld;
    NGA_Access(ga_tmp, lo, hi, &blocktmp, &ld);
    NGA_Access(ga_tmp2, lo, hi, &blockS, &ld);

    int lo1tmp = lo[1];
    double *lambda_vector = (double *)malloc(nfuncs_col * sizeof (double));
    assert (lambda_vector != NULL);   
    #pragma omp parallel for
    #pragma simd
    #pragma ivdep
    for (int j = 0; j < nfuncs_col; j++) {
        lambda_vector[j] = 1.0 / sqrt(eval[j + lo1tmp]);
    }
    free(eval);
    #pragma omp parallel for
    for (int i = 0; i < nfuncs_row; i++)  {
        #pragma simd
        #pragma ivdep
        for (int j = 0; j < nfuncs_col; j++) {
            blockS[i * nfuncs_col + j] =
                blocktmp[i * nfuncs_col + j] * lambda_vector[j];
        }
    }
    free(lambda_vector);
    NGA_Release(ga_tmp, lo, hi);
    NGA_Release_update(ga_tmp2, lo, hi);
    pfock->ga_X = GA_Duplicate(pfock->ga_D[0], "X mat");
    if (0 == pfock->ga_X) {
        PFOCK_PRINTF (1, "GA allocation failed\n");
        return PFOCK_STATUS_ALLOC_FAILED;
    }
    GA_Dgemm('N', 'T', nbf, nbf, nbf,
             1.0, ga_tmp2, ga_tmp, 0.0, pfock->ga_X);

    GA_Destroy(ga_tmp);
    GA_Destroy(ga_tmp2);

    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_destroyOvlMat(PFock_t pfock)
{
    GA_Destroy(pfock->ga_X);
    GA_Destroy(pfock->ga_S);

    return PFOCK_STATUS_SUCCESS;    
}


PFockStatus_t PFock_getOvlMat(PFock_t pfock, int rowstart, int rowend,
                              int colstart, int colend,
                              int stride, double *mat)
{
    int lo[2];
    int hi[2];    
    lo[0] = rowstart;
    hi[0] = rowend;    
    lo[1] = colstart;
    hi[1] = colend;    
    NGA_Get(pfock->ga_S, lo, hi, mat, &stride);

    return PFOCK_STATUS_SUCCESS;    
}


PFockStatus_t PFock_getOvlMat2(PFock_t pfock, int rowstart, int rowend,
                               int colstart, int colend,
                               int stride, double *mat)
{
    int lo[2];
    int hi[2];    
    lo[0] = rowstart;
    hi[0] = rowend;    
    lo[1] = colstart;
    hi[1] = colend;    
    NGA_Get(pfock->ga_X, lo, hi, mat, &stride);

    return PFOCK_STATUS_SUCCESS;    
}*/


PFockStatus_t PFock_getMemorySize(PFock_t pfock, double *mem_cpu)
{
    *mem_cpu = pfock->mem_cpu;  
    return PFOCK_STATUS_SUCCESS;
}


PFockStatus_t PFock_getStatistics(PFock_t pfock)
{
    int myrank;
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    
    // statistics
    MPI_Gather (&pfock->steals, 1, MPI_DOUBLE,
        pfock->mpi_steals, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&pfock->stealfrom, 1, MPI_DOUBLE,
        pfock->mpi_stealfrom, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&pfock->usq, 1, MPI_DOUBLE, 
        pfock->mpi_usq, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&pfock->uitl, 1, MPI_DOUBLE, 
        pfock->mpi_uitl, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&pfock->timepass, 1, MPI_DOUBLE, 
        pfock->mpi_timepass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&pfock->timecomp, 1, MPI_DOUBLE, 
        pfock->mpi_timecomp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&pfock->timeinit, 1, MPI_DOUBLE, 
        pfock->mpi_timeinit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&pfock->timereduce, 1, MPI_DOUBLE, 
        pfock->mpi_timereduce, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&pfock->timegather, 1, MPI_DOUBLE, 
        pfock->mpi_timegather, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&pfock->timescatter, 1, MPI_DOUBLE, 
        pfock->mpi_timescatter, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&pfock->volumega, 1, MPI_DOUBLE, 
        pfock->mpi_volumega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&pfock->ngacalls, 1, MPI_DOUBLE, 
        pfock->mpi_ngacalls, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
        double total_timepass;
        double max_timepass;
        double total_timereduce;
        double total_timeinit;
        double total_timecomp;
        double total_timegather;
        double total_timescatter;
        double total_usq;
        double max_usq;
        double total_uitl;
        double max_uitl;
        double total_steals;
        double total_stealfrom;
        double total_ngacalls;
        double total_volumega;
        for (int i = 0; i < pfock->nprocs; i++) {
            total_timepass += pfock->mpi_timepass[i];
            max_timepass =
                max_timepass < pfock->mpi_timepass[i] ?
                    pfock->mpi_timepass[i] : max_timepass;            
            total_usq += pfock->mpi_usq[i];
            max_usq = 
                max_usq < pfock->mpi_usq[i] ? pfock->mpi_usq[i] : max_usq;          
            total_uitl += pfock->mpi_uitl[i];
            max_uitl = 
                max_uitl < pfock->mpi_uitl[i] ? pfock->mpi_uitl[i] : max_uitl;           
            total_steals += pfock->mpi_steals[i];
            total_stealfrom += pfock->mpi_stealfrom[i];
            total_timecomp += pfock->mpi_timecomp[i];
            total_timeinit += pfock->mpi_timeinit[i];
            total_timereduce += pfock->mpi_timereduce[i];
            total_timegather += pfock->mpi_timegather[i];
            total_timescatter += pfock->mpi_timescatter[i];
            total_ngacalls += pfock->mpi_ngacalls[i];
            total_volumega += pfock->mpi_volumega[i];
        }
        double tsq = pfock->nshells;
        tsq = ((tsq + 1) * tsq/2.0 + 1) * tsq * (tsq + 1)/4.0;
        printf("    PFock Statistics:\n");
        printf("      average totaltime   = %.3g\n"
               "      average timegather  = %.3g\n"
               "      average timeinit    = %.3g\n"
               "      average timecomp    = %.3g\n"
               "      average timereduce  = %.3g\n"
               "      average timescatter = %.3g\n"
               "      comp/total = %.3g\n",
               total_timepass/pfock->nprocs,
               total_timegather/pfock->nprocs,
               total_timeinit/pfock->nprocs,
               total_timecomp/pfock->nprocs,
               total_timereduce/pfock->nprocs,
               total_timescatter/pfock->nprocs,
               total_timecomp/total_timepass);
        printf("      usq = %.4g (lb = %.3g)\n"
               "      uitl = %.4g (lb = %.3g)\n"
               "      nsq = %.4g (screening = %.3g)\n",
               total_usq, max_usq/(total_usq/pfock->nprocs),
               total_uitl, max_uitl/(total_uitl/pfock->nprocs),
               tsq, total_usq/tsq);
        printf("      load blance = %.3lf\n",
               max_timepass/(total_timepass/pfock->nprocs));
        printf("      steals = %.3g (average = %.3g)\n"
               "      stealfrom = %.3g (average = %.3g)\n"
               "      GAcalls = %.3g\n"
               "      GAvolume %.3g MB\n",
               total_steals, total_steals/pfock->nprocs,
               total_stealfrom, total_stealfrom/pfock->nprocs,
               total_ngacalls/pfock->nprocs,
               total_volumega/pfock->nprocs/1024.0/1024.0);
    }
    
    return PFOCK_STATUS_SUCCESS;
}
