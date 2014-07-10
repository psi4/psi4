#ifndef __PFOCK_TYPE_H__
#define __PFOCK_TYPE_H__


#include <omp.h>
#include <ga.h>

#include "CInt.h"


struct PFock
{
    int nthreads;
    
    // screening
    int nnz;
    double *shellvalue;
    int *shellptr;
    int *shellid;
    int *shellrid;
    double maxvalue;
    double tolscr;
    double tolscr2;

    // problem parameters
    int nbf;
    int nshells;
    int natoms;
    int maxnfuncs;

    int nprocs;
    int nprow; // np per row    
    int npcol; // np per col
    // task size
    int nbp_p;
    int nbp_row;
    int nbp_col;
    // pointers to shells
    int *blkrowptr_sh;
    int *blkcolptr_sh;
    // local number of taskes
    int ntasks;
    // pointers to blocks
    int *rowptr_blk;
    int *colptr_blk;
    int *rowptr_sh;
    int *colptr_sh;
    int *rowptr_f;
    int *colptr_f;
    // blks
    int sblk_row;
    int eblk_row;
    int sblk_col;
    int eblk_col;
    int nblks_row;
    int nblks_col;
    // shells
    int sshell_row;
    int sshell_col;
    int eshell_row;
    int eshell_col;
    int nshells_row;
    int nshells_col;
    // functions
    int sfunc_row;
    int sfunc_col;
    int efunc_row;
    int efunc_col;
    int nfuncs_row;
    int nfuncs_col;
    int sizemyrow;
    int sizemycol;
    
    //task queue
    int ga_taskid;
    int icount;

    // integrals
    ERD_t *erd;
    omp_lock_t lock;
    int *f_startind;
    int *s_startind;

    // arrays and buffers
    int ga_screening;
    int maxnumdmat;
    int maxnumdmat2;
    int numdmat;
    int numdmat2;
    int nosymm;
    int committed;
    int *rowptr;
    int *colptr;
    int *rowpos;
    int *colpos;
    int *rowsize;
    int *colsize;
    int *loadrow;
    int sizeloadrow;
    int *loadcol;
    int sizeloadcol;

    // global arrays
    int *ga_F;
    int *ga_D;
    int *ga_J;
    int *ga_K;
    int *gatable[4];

    double *FT_block;
    double *FT_block2;
    // local buffers
    double **D1;
    double **D2;
    double **D3;
    double **VD1;
    double **VD2;
    double **VD3;
    double **F1;
    double **F2;
    double **F3;
    double ***VF1;
    double ***VF2;
    double ***VF3;
    // memory    
    double **VX1;
    double **VX2;
    double **VX3;

    // GA for buf D and F
    int maxrowfuncs;
    int maxcolfuncs;
    int maxrowsize;
    int maxcolsize;
    int *ga_bufD1;
    int *ga_bufD2;
    int *ga_bufD3;
    int *ga_bufF1;
    int *ga_bufF2;
    int *ga_bufF3;
    int *ga_bufX1;
    int *ga_bufX2;
    int *ga_bufX3;
    int lo_D1[2];
    int hi_D1[2];
    int hi_F1[2];
    int lo_D2[2];
    int hi_D2[2];
    int hi_F2[2];
    int lo_D3[2];
    int hi_D3[2];
    int hi_F3[2];
    ga_nbhdl_t *nbD1;
    ga_nbhdl_t *nbD2;
    ga_nbhdl_t *nbD3;
    ga_nbhdl_t *nbF1;
    ga_nbhdl_t *nbF2;
    ga_nbhdl_t *nbF3;
    
    // statistics
    double *mpitime;    
    double *computetime;
    double *usq;
    double *uitl;
    int *steals;
    int *stealfrom;
    char str_buf[512];
};


struct Ovl
{
    int ga;
};


struct CoreH
{
    int ga;
};


#endif /* #ifndef __PFOCK_TYPE_H__ */
