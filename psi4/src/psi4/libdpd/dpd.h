/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here
*/
#ifndef _psi_src_lib_libdpd_dpd_h
#define _psi_src_lib_libdpd_dpd_h

#include <cstdio>
#include <string>
#include "psi4/psifiles.h"
#include "psi4/libpsio/config.h"
 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include <vector>
#include "psi4/psi4-dec.h"

// Testing -TDC
#include "dpdmospace.h"

namespace psi {

#define T3_TIMER_ON (0)

#define DPD_BIGNUM 2147483647 /* the four-byte signed int limit */
/* #define ALL_BUF4_SORT_OOC */

struct dpdparams4{
    int nirreps;      /* No. of irreps */
    int pqnum;        /* Pair number for the row indices */
    int rsnum;        /* Pair number for the column indices */
    int *rowtot;      /* Row dimension for each irrep */
    int *coltot;      /* Column dimension for each irrep */
    int **rowidx;     /* Row index lookup array */
    int **colidx;     /* Column index lookup array */
    int ***roworb;    /* Row index -> orbital index lookup array */
    int ***colorb;    /* Column index -> orbital index lookup array */
    int *ppi;         /* Number of p indices per irrep */
    int *qpi;         /* Number of q indices per irrep */
    int *rpi;         /* Number of r indices per irrep */
    int *spi;         /* Number of s indices per irrep */
    int *poff;        /* Orbital offset for p */
    int *qoff;        /* Orbital offset for q */
    int *roff;        /* Orbital offset for r */
    int *soff;        /* Orbital offset for s */
    int *psym;        /* Orbital symmetry for index p */
    int *qsym;        /* Orbital symmetry for index q */
    int *rsym;        /* Orbital symmetry for index r */
    int *ssym;        /* Orbital symmetry for index s */
    int perm_pq;      /* Can p and q be permuted? */
    int perm_rs;      /* Can r and s be permuted? */
    int peq;          /* Can p and q be equal? */
    int res;          /* Can r and s be equal? */
    int **start13;    /* returns the starting row of irrep matrix h for orbital index p */
};

struct dpdfile4 {
    int dpdnum;                         /* dpd structure reference */
    char label[PSIO_KEYLEN];
    int filenum;
    int my_irrep;     /* Total irrep of this quantity */
    psio_address *lfiles;   /* File address for each submatrix by ROW irrep */
    dpdparams4 *params;
    int incore;
    double ***matrix;
};

struct dpdshift4 {
    int shift_type;
    int **rowtot;
    int **coltot;
    double ****matrix;
};

struct dpdbuf4 {
    int dpdnum;                         /* dpd structure reference */
    int anti;         /* Is this buffer antisymmetric? */
    dpdparams4 *params;
    dpdfile4 file;
    dpdshift4 shift;
    int **row_offset;
    int **col_offset;
    double ***matrix;
};

struct dpdtrans4 {
    double ***matrix;
    dpdshift4 shift;
    dpdbuf4 buf;
};

struct dpdparams2 {
    int nirreps;      /* No. of irreps */
    int pnum;
    int qnum;
    int *rowtot;      /* Row dimension for each submatrix */
    int *coltot;      /* Column dimension for each submatrix */
    int *rowidx;      /* Row index lookup array */
    int *colidx;      /* Column index lookup array */
    int **roworb;     /* Row index -> orbital index lookup array */
    int **colorb;     /* Column index -> orbital index lookup array */
    int *ppi;         /* Number of p indices per irrep */
    int *qpi;         /* Number of q indices per irrep */
    int *poff;        /* Orbital offset for p */
    int *qoff;        /* Orbital offset for q */
    int *psym;        /* Orbital symmetry for index p */
    int *qsym;        /* Orbital symmetry for index q */
};

struct dpdfile2 {
    int dpdnum;                         /* dpd structure reference */
    char label[PSIO_KEYLEN];
    int filenum;
    int my_irrep;
    psio_address *lfiles;
    dpdparams2 *params;
    int incore;
    double ***matrix;
};

/* DPD File4 Cache entries */
struct dpd_file4_cache_entry {
    int dpdnum;                         /* dpd structure reference */
    int filenum;                        /* libpsio unit number */
    int irrep;                          /* overall symmetry */
    int pqnum;                          /* dpd pq value */
    int rsnum;                          /* dpd rs value */
    char label[PSIO_KEYLEN];            /* libpsio TOC keyword */
    double ***matrix;                   /* pointer to irrep blocks */
    int size;                           /* size of entry in double words */
    size_t access;                /* access time */
    size_t usage;                 /* number of accesses */
    size_t priority;              /* priority level */
    int lock;                           /* auto-deletion allowed? */
    int clean;                          /* has this file4 changed? */
    dpd_file4_cache_entry *next; /* pointer to next cache entry */
    dpd_file4_cache_entry *last; /* pointer to previous cache entry */
};

/* DPD File2 Cache entries */
struct dpd_file2_cache_entry {
    dpd_file2_cache_entry():
        next(nullptr),
        last(nullptr)
    {
    }
    int dpdnum;                         /* dpd structure reference */
    int filenum;                        /* libpsio unit number */
    int irrep;                          /* overall symmetry */
    int pnum;                           /* dpd p value */
    int qnum;                           /* dpd q value */
    char label[PSIO_KEYLEN];            /* libpsio TOC keyword */
    double ***matrix;                   /* pointer to irrep blocks */
    int size;                           /* size of entry in double words */
    int clean;                          /* has this file2 changed? */
    dpd_file2_cache_entry *next; /* pointer to next cache entry */
    dpd_file2_cache_entry *last; /* pointer to previous cache entry */
};

/* DPD global parameter set */
struct dpd_data {
    int nirreps;
    int num_subspaces;
    int num_pairs;
    int *numorbs;
    int **orboff;
    int **pairtot;
    int **orbspi;
    int **orbsym;
    int **orbidx2;
    int ***pairidx;
    int ***orbs2;
    int ****pairorb;
    dpdparams2 **params2;
    dpdparams4 **params4;
};

struct thread_data {
    dpdbuf4 *CIjAb; dpdbuf4 *WAbEi; dpdbuf4 *WMbIj; int do_singles; dpdbuf4 *Dints;
    dpdfile2 *SIA; int do_doubles; dpdfile2 *FME; dpdbuf4 *WmAEf; dpdbuf4 *WMnIe;
    dpdbuf4 *SIjAb; int *occpi; int *occ_off; int *virtpi; int *vir_off;
    double omega; dpdfile2 *fIJ; dpdfile2 *fAB; int Gi; int Gj; int Gk;
    int first_ijk; int last_ijk;
    std::string outfile; int thr_id; dpdfile2 SIA_local; dpdbuf4 SIjAb_local;
    int newtrips;
};

struct dpd_gbl {
    long int memory;        /* Total memory requested by the user */
    long int memused;       /* Total memory used (cache + other) */
    long int memcache;      /* Total memory in cache (locked and unlocked) */
    long int memlocked;     /* Total memory locked in the cache */

    // The default C'tor will zero everything out properly
    dpd_gbl():
        file2_cache(nullptr),
        file4_cache(nullptr),
        file4_cache_most_recent(0),
        file4_cache_least_recent(1),
        file4_cache_lru_del(0),
        file4_cache_low_del(0)
    {}
    dpd_file2_cache_entry *file2_cache;
    dpd_file4_cache_entry *file4_cache;
    size_t file4_cache_most_recent;
    size_t file4_cache_least_recent;
    size_t file4_cache_lru_del;
    size_t file4_cache_low_del;
    int cachetype;
    int *cachefiles;
    int **cachelist;
    dpd_file4_cache_entry *file4_cache_priority;
};

/* Useful for the generalized 4-index sorting function */
enum indices {pqrs, pqsr, prqs, prsq, psqr, psrq,
              qprs, qpsr, qrps, qrsp, qspr, qsrp,
              rqps, rqsp, rpqs, rpsq, rsqp, rspq,
              sqrp, sqpr, srqp, srpq, spqr, sprq};

/* Useful for the 3-index sorting function dpd_3d_sort() */
enum pattern {abc, acb, cab, cba, bca, bac};

class PSI_API DPD{
public:

    // These used to live in the dpd_data struct
    int nirreps;
    int num_subspaces;
    int num_pairs;
    int *numorbs;
    int **orboff;
    int **pairtot;
    int **orbspi;
    int **orbsym;
    int **orbidx2;
    int ***pairidx;
    int ***orbs2;
    int ****pairorb;
    dpdparams2 **params2;
    dpdparams4 **params4;

   std::vector<DPDMOSpace> moSpaces;

    DPD(int dpd_num, int nirreps, long int memory, int cachetype,
        int *cachefiles, int **cachelist, dpd_file4_cache_entry *priority,
        int num_subspaces, std::vector<int*> &spaceArrays);
    DPD(int dpd_num, int nirreps, long int memory, int cachetype,
        int *cachefiles, int **cachelist, dpd_file4_cache_entry *priority,
        int num_subspaces, std::vector<DPDMOSpace> &moSpaces);
    DPD();

    ~DPD();

    int init(int dpd_num, int nirreps, long int memory, int cachetype,
                 int *cachefiles, int **cachelist,
                 dpd_file4_cache_entry *priority, int num_subspaces, ...);
    int init(int dpd_num, int nirreps, long int memory, int cachetype,
                 int *cachefiles, int **cachelist, dpd_file4_cache_entry *priority,
                 int num_subspaces, std::vector<int*> &spaceArrays);

    void dpd_error(const char *caller, std::string out_fname);

    double **dpd_block_matrix(size_t n, size_t m);
    void free_dpd_block(double **array, size_t n, size_t m);

    int contract222(dpdfile2 *X, dpdfile2 *Y, dpdfile2 *Z, int target_X,
                    int target_Y, double alpha, double beta);
    int contract442(dpdbuf4 *X, dpdbuf4 *Y, dpdfile2 *Z, int target_X,
                    int target_Y, double alpha, double beta);
    int contract422(dpdbuf4 *X, dpdfile2 *Y, dpdfile2 *Z, int trans_Y,
                    int trans_Z, double alpha, double beta);
    int contract244(dpdfile2 *X, dpdbuf4 *Y, dpdbuf4 *Z, int sum_X, int sum_Y,
                    int trans_Z, double alpha, double beta);
    int contract424(dpdbuf4 *X, dpdfile2 *Y, dpdbuf4 *Z, int sum_X,
                    int sum_Y, int trans_Z, double alpha, double beta);
    int contract444(dpdbuf4 *X, dpdbuf4 *Y, dpdbuf4 *Z,
                    int target_X, int target_Y, double alpha, double beta);
    int contract444_df(dpdbuf4 *B, dpdbuf4 *tau_in, dpdbuf4 *tau_out, double alpha, double beta);

    /* Need to consolidate these routines into one general function */
    int dot23(dpdfile2 *T, dpdbuf4 *I, dpdfile2 *Z,
              int transt, int transz, double alpha, double beta);
    int dot24(dpdfile2 *T, dpdbuf4 *I, dpdfile2 *Z,
              int transt, int transz, double alpha, double beta);
    int dot13(dpdfile2 *T, dpdbuf4 *I, dpdfile2 *Z,
              int transt, int transz, double alpha, double beta);
    int dot14(dpdfile2 *T, dpdbuf4 *I, dpdfile2 *Z,
              int transt, int transz, double alpha, double beta);

    int trace42_13(dpdbuf4 *A, dpdfile2 *B, int transb, double alpha, double beta);

    int file2_init(dpdfile2 *File, int filenum, int irrep, int pnum,
                   int qnum, const char *label);
    int file2_close(dpdfile2 *File);
    int file2_mat_init(dpdfile2 *File);
    int file2_mat_close(dpdfile2 *File);
    int file2_mat_rd(dpdfile2 *File);
    int file2_mat_wrt(dpdfile2 *File);
    int file2_print(dpdfile2 *File, std::string out_fname);
    int file2_mat_print(dpdfile2 *File, std::string out_fname);
    int file2_copy(dpdfile2 *InFile, int outfilenum, const char *label);
    int file2_dirprd(dpdfile2 *FileA, dpdfile2 *FileB);
    double file2_dot(dpdfile2 *FileA, dpdfile2 *FileB);
    int file2_scm(dpdfile2 *InFile, double alpha);
    double file2_dot_self(dpdfile2 *BufX);
    double file2_trace(dpdfile2 *InFile);
    int file2_axpy(dpdfile2 *FileA, dpdfile2 *FileB,
                   double alpha, int transA);
    int file2_axpbycz(dpdfile2 *FileA, dpdfile2 *FileB, dpdfile2 *FileC,
                      double a, double b, double c);


    int file4_init(dpdfile4 *File, int filenum, int irrep, int pqnum,
                   int rsnum,  const char *label);
    int file4_init_nocache(dpdfile4 *File, int filenum, int irrep, int pqnum,
                           int rsnum,  const char *label);
    int file4_close(dpdfile4 *File);
    int file4_mat_irrep_init(dpdfile4 *File, int irrep);
    int file4_mat_irrep_close(dpdfile4 *File, int irrep);
    int file4_mat_irrep_rd(dpdfile4 *File, int irrep);
    int file4_mat_irrep_wrt(dpdfile4 *File, int irrep);
    int file4_mat_irrep_row_init(dpdfile4 *File, int irrep);
    int file4_mat_irrep_row_close(dpdfile4 *File, int irrep);
    int file4_mat_irrep_row_rd(dpdfile4 *File, int irrep, int row);
    int file4_mat_irrep_row_wrt(dpdfile4 *File, int irrep, int row);
    int file4_mat_irrep_row_zero(dpdfile4 *File, int irrep, int row);
    int file4_print(dpdfile4 *File, std::string out_fname);
    int file4_mat_irrep_rd_block(dpdfile4 *File, int irrep, int start_pq,
                                 int num_pq);
    int file4_mat_irrep_wrt_block(dpdfile4 *File, int irrep, int start_pq,
                                  int num_pq);

    int buf4_init(dpdbuf4 *Buf, int inputfile, int irrep, int pqnum, int rsnum,
                  int file_pqnum, int file_rsnum, int anti, const char *label);
    int buf4_init(dpdbuf4 *Buf, int inputfile, int irrep, std::string pq, std::string rs,
                  std::string file_pq, std::string file_rs, int anti, const char *label);
    int buf4_init(dpdbuf4 *Buf, int inputfile, int irrep, std::string pq, std::string rs, int anti, const char *label);
    int pairnum(std::string);
    double buf4_trace(dpdbuf4 *Buf);
    int buf4_close(dpdbuf4 *Buf);
    int buf4_mat_irrep_init(dpdbuf4 *Buf, int irrep);
    int buf4_mat_irrep_close(dpdbuf4 *Buf, int irrep);
    int buf4_mat_irrep_rd(dpdbuf4 *Buf, int irrep);
    int buf4_mat_irrep_wrt(dpdbuf4 *Buf, int irrep);
    int buf4_print(dpdbuf4 *Buf, std::string out_fname, int print_data);
    int buf4_copy(dpdbuf4 *InBuf, int outfilenum, const char *label);
    int buf4_sort(dpdbuf4 *InBuf, int outfilenum, enum indices index,
                  int pqnum, int rsnum, const char *label);
    int buf4_sort(dpdbuf4 *InBuf, int outfilenum, enum indices index,
                    std::string pq, std::string rs, const char *label);
    int buf4_sort_ooc(dpdbuf4 *InBuf, int outfilenum, enum indices index,
                      int pqnum, int rsnum, const char *label);
    int buf4_sort_axpy(dpdbuf4 *InBuf, int outfilenum, enum indices index,
                       int pqnum, int rsnum, const char *label, double alpha);
    int buf4_axpy(dpdbuf4 *BufX, dpdbuf4 *BufY, double alpha);
    int buf4_axpbycz(dpdbuf4 *FileA, dpdbuf4 *FileB, dpdbuf4 *FileC,
                     double a, double b, double c);
    int buf4_dirprd(dpdbuf4 *BufA, dpdbuf4 *BufB);
    double buf4_dot(dpdbuf4 *BufA, dpdbuf4 *BufB);
    double buf4_dot_self(dpdbuf4 *BufX);
    int buf4_scm(dpdbuf4 *InBuf, double alpha);
    int buf4_scmcopy(dpdbuf4 *InBuf, int outfilenum, const char *label,
                     double alpha);
    int buf4_symm(dpdbuf4 *Buf);
    int buf4_symm2(dpdbuf4 *Buf1, dpdbuf4 *Buf2);
    int buf4_mat_irrep_shift13(dpdbuf4 *Buf, int irrep);
    int buf4_mat_irrep_shift31(dpdbuf4 *Buf, int irrep);
    int buf4_mat_irrep_row_init(dpdbuf4 *Buf, int irrep);
    int buf4_mat_irrep_row_close(dpdbuf4 *Buf, int irrep);
    int buf4_mat_irrep_row_zero(dpdbuf4 *Buf, int irrep, int row);
    int buf4_mat_irrep_row_rd(dpdbuf4 *Buf, int irrep, int pq);
    int buf4_mat_irrep_row_wrt(dpdbuf4 *Buf, int irrep, int pq);
    int buf4_mat_irrep_init_block(dpdbuf4 *Buf, int irrep, int num_pq);
    int buf4_mat_irrep_close_block(dpdbuf4 *Buf, int irrep, int num_pq);
    int buf4_mat_irrep_rd_block(dpdbuf4 *Buf, int irrep, int start_pq,
                                int num_pq);
    int buf4_mat_irrep_wrt_block(dpdbuf4 *Buf, int irrep, int start_pq,
                                 int num_pq);
    int buf4_dump(dpdbuf4 *DPDBuf, struct iwlbuf *IWLBuf,
                  int *prel, int *qrel, int *rrel, int *srel,
                  int bk_pack, int swap23);
    int trans4_init(dpdtrans4 *Trans, dpdbuf4 *Buf);
    int trans4_close(dpdtrans4 *Trans);
    int trans4_mat_irrep_init(dpdtrans4 *Trans, int irrep);
    int trans4_mat_irrep_close(dpdtrans4 *Trans, int irrep);
    int trans4_mat_irrep_rd(dpdtrans4 *Trans, int irrep);
    int trans4_mat_irrep_wrt(dpdtrans4 *Trans, int irrep);
    int trans4_mat_irrep_shift13(dpdtrans4 *Trans, int irrep);
    int trans4_mat_irrep_shift31(dpdtrans4 *Trans, int irrep);

    int mat4_irrep_print(double **matrix, dpdparams4 *Params,
                         int irrep, int my_irrep, std::string out_fname);

    void file2_cache_init(void);
    void file2_cache_close(void);
    void file2_cache_print(std::string out_fname);
    dpd_file2_cache_entry* file2_cache_scan(int filenum, int irrep, int pnum, int qnum, const char *label, int dpdnum);
    dpd_file2_cache_entry* dpd_file2_cache_last(void);
    int file2_cache_add(dpdfile2 *File);
    int file2_cache_del(dpdfile2 *File);
    int file4_cache_del_low(void);
    void file2_cache_dirty(dpdfile2 *File);

    void file4_cache_init(void);
    void file4_cache_close(void);
    void file4_cache_print(std::string out_fname);
    void file4_cache_print_screen(void);
    int file4_cache_get_priority(dpdfile4 *File);

    dpd_file4_cache_entry* file4_cache_scan(int filenum, int irrep, int pqnum, int rsnum, const char *label, int dpdnum);
    dpd_file4_cache_entry* file4_cache_last(void);
    int file4_cache_add(dpdfile4 *File, size_t priority);
    int file4_cache_del(dpdfile4 *File);
    dpd_file4_cache_entry* file4_cache_find_lru(void);
    int file4_cache_del_lru(void);
    void file4_cache_dirty(dpdfile4 *File);
    void file4_cache_lock(dpdfile4 *File);
    void file4_cache_unlock(dpdfile4 *File);

    void sort_3d(double ***Win, double ***Wout, int nirreps, int h, int *rowtot, int **rowidx,
                 int ***roworb, int *asym, int *bsym, int *aoff, int *boff,
                 int *cpi, int *coff, int **rowidx_out, enum pattern index, int sum);


    void T3_AAA(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk,
                dpdbuf4 *T2, dpdbuf4 *F, dpdbuf4 *E, dpdfile2 *fIJ, dpdfile2 *fAB,
                int *occpi, int *occ_off, int *virtpi, int *vir_off, double omega);

    void T3_AAB(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk,
                dpdbuf4 *T2AA, dpdbuf4 *T2AB, dpdbuf4 *T2BA, dpdbuf4 *FAA, dpdbuf4 *FAB, dpdbuf4 *FBA,
                dpdbuf4 *EAA, dpdbuf4 *EAB, dpdbuf4 *EBA, dpdfile2 *fIJ, dpdfile2 *fij,
                dpdfile2 *fAB, dpdfile2 *fab, int *aoccpi, int *aocc_off, int *boccpi, int *bocc_off,
                int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off, double omega);

    void T3_RHF(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk,
                dpdbuf4 *T2, dpdbuf4 *F, dpdbuf4 *E, dpdfile2 *fIJ, dpdfile2 *fAB,
                int *occpi, int *occ_off, int *virtpi, int *vir_off, double omega);

    void T3_RHF_ic(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk,
                   dpdbuf4 *T2, dpdbuf4 *F, dpdbuf4 *E, dpdfile2 *fIJ, dpdfile2 *fAB,
                   int *occpi, int *occ_off, int *virtpi, int *vir_off, double omega);

    void cc3_sigma_RHF(dpdbuf4 *CIjAb, dpdbuf4 *WAbEi, dpdbuf4 *WMbIj,
                       int do_singles, dpdbuf4 *Dints, dpdfile2 *SIA,
                       int do_doubles, dpdfile2 *FME, dpdbuf4 *WAmEf, dpdbuf4 *WMnIe,
                       dpdbuf4 *SIjAb, int *occpi, int *occ_off, int *virtpi, int *vir_off,
                       double omega,std::string out_fname, int newtrips);

    void cc3_sigma_RHF_ic(dpdbuf4 *CIjAb, dpdbuf4 *WAbEi, dpdbuf4 *WMbIj,
                          int do_singles, dpdbuf4 *Dints, dpdfile2 *SIA,
                          int do_doubles, dpdfile2 *FME, dpdbuf4 *WAmEf, dpdbuf4 *WMnIe,
                          dpdbuf4 *SIjAb, int *occpi, int *occ_off, int *virtpi, int *vir_off,
                          double omega, std::string out_fname, int nthreads, int newtrips);

    void cc3_sigma_UHF_AAA(dpdbuf4 *CMNEF, dpdbuf4 *WABEI, dpdbuf4 *WMBIJ,
                           int do_singles, dpdbuf4 *Dints_anti, dpdfile2 *SIA, int do_doubles,
                           dpdfile2 *FME, dpdbuf4 *WMAFE, dpdbuf4 *WMNIE, dpdbuf4 *SIJAB,
                           int *aoccpi, int *aocc_off, int *avirtpi, int *avir_off, double omega,
                           std::string out_fname);

    void cc3_sigma_UHF_BBB(dpdbuf4 *Cmnef, dpdbuf4 *Wabei, dpdbuf4 *Wmbij,
                           int do_singles, dpdbuf4 *Dijab_anti, dpdfile2 *Sia, int do_doubles,
                           dpdfile2 *Fme, dpdbuf4 *Wmafe, dpdbuf4 *Wmnie, dpdbuf4 *Sijab,
                           int *boccpi, int *bocc_off, int *bvirtpi, int *bvir_off, double omega,
                           std::string out_fname);

    void cc3_sigma_UHF_AAB(dpdbuf4 *C2AA, dpdbuf4 *C2AB, dpdbuf4 *C2BA,
                           dpdbuf4 *FAA, dpdbuf4 *FAB, dpdbuf4 *FBA,
                           dpdbuf4 *EAA, dpdbuf4 *EAB, dpdbuf4 *EBA,
                           int do_singles, dpdbuf4 *DAA, dpdbuf4 *DAB, dpdfile2 *SIA, dpdfile2 *Sia,
                           int do_doubles, dpdfile2 *FME, dpdfile2 *Fme,
                           dpdbuf4 *WMAFE, dpdbuf4 *WMaFe, dpdbuf4 *WmAfE,
                           dpdbuf4 *WMNIE, dpdbuf4 *WMnIe, dpdbuf4 *WmNiE,
                           dpdbuf4 *SIJAB, dpdbuf4 *SIjAb, int *aoccpi, int *aocc_off, int *boccpi,
                           int *bocc_off, int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off,
                           double omega, std::string out_fname);

    void cc3_sigma_UHF_BBA(dpdbuf4 *C2BB, dpdbuf4 *C2AB, dpdbuf4 *C2BA,
                           dpdbuf4 *FBB, dpdbuf4 *FAB, dpdbuf4 *FBA,
                           dpdbuf4 *EBB, dpdbuf4 *EAB, dpdbuf4 *EBA,
                           int do_singles, dpdbuf4 *DBB, dpdbuf4 *DBA, dpdfile2 *SIA, dpdfile2 *Sia,
                           int do_doubles, dpdfile2 *FME, dpdfile2 *Fme,
                           dpdbuf4 *Wmafe, dpdbuf4 *WMaFe, dpdbuf4 *WmAfE,
                           dpdbuf4 *Wmnie, dpdbuf4 *WMnIe, dpdbuf4 *WmNiE,
                           dpdbuf4 *Sijab, dpdbuf4 *SIjAb, int *aoccpi, int *aocc_off, int *boccpi,
                           int *bocc_off, int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off,
                           double omega, std::string out_fname);
}; // Dpd class


/*
 * Static variables/functions to mimic the old C machinery
 */
extern dpd_gbl dpd_main;
extern PSI_API DPD* global_dpd_;
extern int dpd_default;
extern DPD* dpd_list[2];
extern PSI_API int dpd_set_default(int dpd_num);
extern int dpd_init(int dpd_num, int nirreps, long int memory, int cachetype,
            int *cachefiles, int **cachelist, dpd_file4_cache_entry *priority,
            int num_subspaces, std::vector<int*> &spaceArrays);
extern int dpd_close(int dpd_num);
extern long int dpd_memfree(void);
extern void dpd_memset(long int memory);


}// Namespace psi

#endif /* _psi_src_lib_libdpd_dpd_h */
