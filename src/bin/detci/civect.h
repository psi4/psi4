/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/
#ifndef _psi_src_bin_detci_civect_h
#define _psi_src_bin_detci_civect_h

namespace psi { namespace detci {

//typedef unsigned long long int BIGINT;
typedef unsigned long int BIGINT;

/*
** CIVECT.H
**
** Header file for CI Vector Class
**
** David Sherrill, 15 June 1995
** Center for Computational Quantum Chemistry
**
*/

class CIvect {
   protected:
      BIGINT vectlen;            /* total vector length */ 
      unsigned long buffer_size; /* size of largest in-core chunk */
      int num_blocks;            /* number of blocks in vector */
      int icore;                 /* in-core option. 1 = whole vector in-core,
                                                    2 = symm block in-core,
                                                    0 = one subblock in-core */
      int Ms0;                   /* is this a vector for M_s=0? 1=yes, 0=no */
      int *Ia_code;              /* each block's alpha string id code */
      int *Ib_code;              /* each block's beta string id code */
      int *Ia_size;              /* num alp strings in each block */
      int *Ib_size;              /* num bet strings in each block */
      BIGINT *offset;            /* offsets for absolute numbering.  This
                                    is a word offset, not a byte offset,
                                    so unsigned long should be ok          */
      int num_alpcodes;          /* number of possible (total) alpha codes */
      int num_betcodes;          /* number of possible (total) beta codes */
      int nirreps;               /* dimension of next four arrays */
      int codes_per_irrep;       /* number of codes per irrep (for alpha) */
      int buf_per_vect;          /* number of buffers per CI vector */
      int buf_total;             /* number of total buffers (all vectors) */
      int new_first_buf;         /* after collapse, buffs get renumbered */
      int maxvect;               /* max number of CI vectors */
      int nvect;                 /* number of CI vectors */
      int nunits;                /* number of units (physical disk files) */
      int cur_vect;              /* current vector number in core */
      int cur_buf;               /* current buffer in core */
      int buf_locked;            /* is a memory buffer locked in/available?  */
      int *units;                /* file numbers */
      int *file_number;          /* unit number for given vector/block */
      unsigned long *buf_size;   /* size of each buffer on disk 
                                    (0...buf_per_vect) */
      int *buf2blk;              /* buffer number -> block number for
                                    icore=0, else buf->irrep for icore=2 */
      int *buf_offdiag;          /* is the buffer "off-diagonal"? only applies
                                    to Ms=0.  If Ms<>0, always=0 */
      int *first_ablk;           /* first blocknum with a given Ia irrep */
      int *last_ablk;            /* last blocknum with a given Ia irrep */
      int **decode;              /* gives block number for a (alp,bet) code */
                                 /* dimensions num_alpcodes * num_betcodes */
      double ***blocks;          /* a matrix for each block                */
      double *buffer;            /* pointer to buffer, same as blocks[0][0] */
      int *zero_blocks;          /* array for which blocks happen to be 0   */
      int in_file;               /* increment for how many buffers in a file */
      int extras;                /* accounts for extra buffers */
      int units_used;            /* accounts for number of unit files used */  
      int cur_unit;              /* current unit file */
      int cur_size;              /* current size of buffer */
      int first_unit;            /* first file unit number (if > 1) */ 
      
   public:
      CIvect();
      CIvect(BIGINT vl, int nb, int incor, int ms0, int *iac,
         int *ibc, int *ias, int *ibs, BIGINT *offs, int nac, int nbc, 
         int nirr, int cdperirr, int maxvect, int nunits, 
         int funit, int *fablk, int *lablk, int **dc);
      ~CIvect();

      double * buf_malloc(void);
      void set(BIGINT vl, int nb, int incor, int ms0, int *iac,
         int *ibc, int *ias, int *ibs, BIGINT *offs, int nac, int nbc, 
         int nirr, int cdperirr, int maxvect, int nunits, int funit, 
         int *fablk, int *lablk, int **dc);
      void print(FILE *outfile);
      double operator*(CIvect &b);
      void set_nvect(int i);
      void setarray(const double *a, int len);
      void max_abs_vals(int nval, int *iac, int *ibc, int *iaidx, int *ibidx,
         double *coeff, int neg_only);
      double blk_max_abs_vals(int i, int offdiag, int nval, int *iac, int *ibc,
         int *iaidx, int *ibidx, double *coeff, double minval, int neg_only);
      void det2strings(BIGINT det, int *alp_code, int *bet_code,
         int *alp_idx, int *bet_idx);
      BIGINT strings2det(int alp_code, int alp_idx,
         int bet_code, int bet_idx);
      void diag_mat_els(struct stringwr **alplist, struct stringwr
         **betlist, double *oei, double *tei, double efzc, int na, int nb, 
         int nbf, int method);
      void diag_mat_els_otf(struct stringwr **alplist, struct stringwr
         **betlist, double *oei, double *tei, double efzc, int na, int nb, 
         int nbf, int buf, int method);
      void init_vals(int ivect, int nvals, int *alplist, int *alpidx, 
         int *betlist, int *betidx, int *blknums, double *value);
      void set_vals(int ivect, int nvals, int *alplist, int *alpidx, 
         int *betlist, int *betidx, int *blknums, double *value);
      void extract_vals(int ivect, int nvals, int *alplist, int *alpidx, 
         int *betlist, int *betidx, int *blknums, double *value);
      void symnorm(double a, int vecode, int gather_vec);
      double zero_det(int iac, int ia, int ibc, int ib);
      void scale(double a, int vecode, int gather_vec);
      void symmetrize(double phase, int iblock);
      void buf_lock(double *a);
      void buf_unlock(void);
      double ** blockptr(int blknum);
      void init_io_files(bool open_old);
      void close_io_files(int keep);
      int read(int ivect, int ibuf);
      int write(int ivect, int ibuf);
      int schmidt_add(CIvect &c, int L);
      int schmidt_add2(CIvect &c, int first_vec, int last_vec, int source_vec,
          int target_vec, double *dotval, double *nrm, double *ovlpmax);
      void zero(void);
      void dcalc(int nr, int L, double **alpha, double *lambda,
         double *norm_arr, CIvect &C, CIvect &S, double *buf1, double *buf2, 
         int *root_converged, int printflag, FILE *outfile, double *E_est);
      void sigma_renorm(int nr, int L, double renorm_C, CIvect &S, 
         double *buf1, int printflag, FILE *outfile);
      double dcalc2(int rootnum, double lambda, CIvect &Hd, 
           int precon, struct stringwr **alplist, struct stringwr **betlist);
      double dcalc_evangelisti(int rootnum, int num_vecs, double lambda, 
           CIvect &Hd, CIvect &C, double *buf1, double *buf2, int precon, 
           int L, struct stringwr **alplist, struct stringwr **betlist, 
           double **alpha);
      void construct_kth_order_wf(CIvect &Hd, CIvect &S, CIvect &C, 
           struct stringwr **alplist, struct stringwr **betlist, double *buf1,
           double *buf2, int k, double *mp_energy, double **bvec_overlap,
           double *bvec_norm);
      void wigner_E2k_formula(CIvect &Hd, CIvect &S, CIvect &C, 
           struct stringwr **alplist, struct stringwr **betlist, double *buf1,
           double *buf2, int k, double *mp2k_energy, double **wfn_overlap,
           double **bvec_overlap, double *bvec_norm, int kvec_offset);
      void print_buf(FILE *outfile);
      void civ_xeay(double a, CIvect &Y, int xvect, int yvect);
      void civ_xpeay(double a, CIvect &Y, int xvect, int yvect);
      void transp_block(int iblock, double **tmparr);
      unsigned long get_max_blk_size(void);
      double checknorm(void);
      void copy(CIvect &Src, int targetvec, int srcvec);
      void restart_gather(int ivec, int nvec, int nroot, double **alpha,
         double *buffer1, double *buffer2);
      void gather(int ivec, int nvec, int nroot, double **alpha, CIvect &C);
      void restart_reord_fp(int L);
      void print_fptrs(void);
      double calc_ssq(double *buffer1, double *buffer2, 
           struct stringwr **alplist, struct stringwr **betlist, int vec_num);
      void h0block_buf_init(void);
      void h0block_buf_ols(double *norm,double *ovrlap,double *c1norm,
                           double E_est);
      void h0block_buf_precon(double *norm, int root);
      void h0block_gather_vec(int vecode);
      void h0block_gather_multivec(double *vec);
      int check_zero_block(int blocknum);
      void set_zero_block(int blocknum, int value);
      void set_zero_blocks_all(void);
      void copy_zero_blocks(CIvect &src);
      void print_zero_blocks(void);
      void scale_sigma(CIvect &Hd, CIvect &C,
        struct stringwr **alplist, struct stringwr **betlist, int i, 
        double *buf1, double *buf2);
      int read_new_first_buf(void);
      void write_new_first_buf(void);
      void set_new_first_buf(int nfb);
      int read_num_vecs(void);
      void write_num_vecs(int nv);
      void write_toc(void);
      void civect_psio_debug(void);
      void pt_correction(struct stringwr **alplist, struct stringwr
        **betlist);
      double compute_follow_overlap(int troot, int ncoef, double *coef,
        int *Iac, int *Iaridx, int *Ibc, int *Ibridx);

      friend void sigma_init(CIvect& C, CIvect &S, struct stringwr **alplist, 
         struct stringwr **betlist);
      friend void sigma(struct stringwr **alplist, struct stringwr **betlist,
         CIvect& C, CIvect& S, double *oei, double *tei, int fci, int iter);
      friend void sigma_a(struct stringwr **alplist, struct stringwr **betlist,
         CIvect& C, CIvect& S, double *oei, double *tei, int fci, int iter);
      friend void sigma_b(struct stringwr **alplist, struct stringwr **betlist,
         CIvect& C, CIvect& S, double *oei, double *tei, int fci, int iter);
      friend void sigma_c(struct stringwr **alplist, struct stringwr **betlist,
         CIvect& C, CIvect& S, double *oei, double *tei, int fci, int iter);
      friend void sigma_get_contrib(struct stringwr **alplist, struct
         stringwr **betlist, CIvect &C, CIvect &S, int **s1_contrib, 
         int **s2_contrib, int **s3_contrib);
      friend void sigma_get_contrib_rotf(CIvect &C, CIvect &S, 
         int **s1_contrib, int **s2_contrib, int **s3_contrib,
         int *Jcnt[2], int **Jij[2], int **Joij[2], int **Jridx[2],
         signed char **Jsgn[2], unsigned char **Toccs);
      friend void olsen_iter_xy(CIvect &C, CIvect &S, CIvect &Hd, double *x, 
         double *y, double *buf1, double *buf2, double E, int curvect, 
         int L, double **alpha, struct stringwr **alplist, 
         struct stringwr **betlist);
      friend void olsen_update(CIvect &C, CIvect &S, CIvect &Hd, double E, 
         double E_est, double *norm, double *c1norm, double *ovrlap, 
         double *buffer1, double *buffer2, int curr, int next, FILE *outfile, 
         int iter, struct stringwr **alplist, struct stringwr **betlist);
      friend void mitrush_update(CIvect &C, CIvect &S, double norm, double
         acur, double alast, double *buffer1, double *buffer2, int curr,
         int next);
      friend void opdm(struct stringwr **alplist, struct stringwr **betlist,
          int transdens, int dipmom,
          int Inroots, int Iroot, int Inunits, int Ifirstunit,
          int Jnroots, int Jroot, int Jnunits, int Jfirstunit,
          int targetfile, int writeflag, int printflag);
      friend void tpdm(struct stringwr **alplist, struct stringwr **betlist,
          int Inroots, int Inunits, int Ifirstunit,
          int Jnroots, int Jnunits, int Jfirstunit,
          int targetfile, int writeflag, int printflag);
};

}} // namespace psi::detci

#endif  // header guard

