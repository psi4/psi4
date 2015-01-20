#pragma once

#include <stdint.h>
#include <stdbool.h>

/*struct OED
{
    int nalpha;
    int ncoeff;
    int ncgto1;
    int ncgto2;
    int npgto1;
    int npgto2;
    int shell1;
    int shell2;
    int natoms;
    int ncsum;
    int spheric;
    int screen;  
    double x1;
    double y1;
    double z1;
    double x2;
    double y2;
    double z2;
    double *xn;
    double *yn;
    double *zn;
    double *charge;
    double *cc;
    double *alpha;
    int cc_beg[2];
    int cc_end[2];
    int imax;
    int zmax;
    double *zcore;
    double *zcore2;
    int *icore;

    int fp_memory_opt;
    int int_memory_opt;
    int *coef_offset;
    int *exp_offset;
};


struct ERD
{
    /* The number of threads used for computation /
    uint32_t nthreads;
    size_t capacity;
    double **buffer;    
    /* Used for vrrtable /
    int max_shella;
    /* 2D array /
    int **vrrtable;
#ifdef __INTEL_OFFLOAD
    int mic_numdevs;
#endif    
};*/

struct BasisSet
{
    // atom
    int natoms;
    int *eid;
    double *xn;
    double *yn;
    double *zn;   
    double *charge;
    int nelectrons;
    double **guess;

    double ene_nuc;

    // basis
    int bs_nelements;
    int bs_natoms;
    int basistype;
    int *bs_eid;
    int *bs_eptr;
    int *bs_atom_start;
    int bs_nshells;
    int bs_totnexp;
    int *bs_nexp;
    double **bs_exp;
    double **bs_cc;
    double **bs_norm;
    int *bs_momentum;
    
    // shell
    uint32_t nshells;    
    uint32_t nfunctions;    
    uint32_t *f_start_id;
    uint32_t *f_end_id;
    uint32_t *s_start_id;
    uint32_t *nexp;
    double **exp;
    double *minexp;
    double **cc;
    double **norm;
    uint32_t *momentum;
    double *xyz0;

    uint32_t maxdim; // max number of functions of a shell 
    uint32_t max_momentum;
    uint32_t max_nexp;
    uint32_t max_nexp_id;
    
    char str_buf[512];

#ifdef __INTEL_OFFLOAD
    int mic_numdevs;
#endif
};
