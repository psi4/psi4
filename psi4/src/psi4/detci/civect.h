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
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here
*/
#ifndef _psi_src_bin_detci_civect_h
#define _psi_src_bin_detci_civect_h

#include "psi4/pybind11.h"

// Forward declarations
namespace psi { namespace detci {
typedef unsigned long int BIGINT;
struct calcinfo;
struct params;
struct H_zero_block;
struct ci_blks;
class CIWavefunction;
class CIvect;
typedef std::shared_ptr<psi::detci::CIvect> SharedCIVector;
}}

namespace psi {
namespace detci {

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
  friend class CIWavefunction;

  protected:

    // Holds pointers to relevant CI structs
    struct calcinfo *CI_CalcInfo_;
    struct params *CI_Params_;
    struct H_zero_block *CI_H0block_;

    void common_init(); /* common init func */

    BIGINT vectlen_;            /* total vector length */
    unsigned long buffer_size_; /* size of largest in-core chunk */
    int num_blocks_;            /* number of blocks in vector */
    int icore_;                 /* in-core option. 1 = whole vector in-core,
                                                   2 = symm block in-core,
                                                   0 = one subblock in-core */
    int Ms0_;                   /* is this a vector for M_s=0? 1=yes, 0=no */
    std::vector<int> Ia_code_;   /* each block's alpha string id code */
    std::vector<int> Ib_code_;   /* each block's beta string id code */
    std::vector<int> Ia_size_;   /* num alp strings in each block */
    std::vector<int> Ib_size_;   /* num bet strings in each block */
    std::vector<BIGINT> offset_; /* offsets for absolute numbering.  This
                                   is a word offset, not a byte offset,
                                   so unsigned long should be ok          */
    int num_alpcodes_;          /* number of possible (total) alpha codes */
    int num_betcodes_;          /* number of possible (total) beta codes */
    int nirreps_;               /* dimension of next four arrays */
    int codes_per_irrep_;       /* number of codes per irrep (for alpha) */
    int buf_per_vect_;          /* number of buffers per CI vector */
    int buf_total_;             /* number of total buffers (all vectors) */
    int new_first_buf_;         /* after collapse, buffs get renumbered */
    int maxvect_;               /* max number of CI vectors */
    int nvect_;                 /* number of CI vectors */
    int nunits_;                /* number of units (physical disk files) */
    int cur_vect_;              /* current vector number in core */
    int cur_buf_;               /* current buffer in core */
    int buf_locked_;            /* is a memory buffer locked in/available?  */
    std::vector<int> units_;                /* file numbers */
    std::vector<int> file_number_;          /* unit number for given vector/block */
    unsigned long *buf_size_;   /* size of each buffer on disk
                                   (0...buf_per_vect) */
    int *buf2blk_;              /* buffer number -> block number for
                                   icore=0, else buf->irrep for icore=2 */
    int *buf_offdiag_;          /* is the buffer "off-diagonal"? only applies
                                   to Ms=0.  If Ms<>0, always=0 */
    int *first_ablk_;           /* first blocknum with a given Ia irrep */
    int *last_ablk_;            /* last blocknum with a given Ia irrep */
    int **decode_;              /* gives block number for a (alp,bet) code */
                                /* dimensions num_alpcodes * num_betcodes */
    double ***blocks_;          /* a matrix for each block                */
    double *buffer_;            /* pointer to buffer, same as blocks[0][0] */
    std::vector<int> zero_blocks_;          /* array for which blocks happen to be 0   */
    int in_file_;               /* increment for how many buffers in a file */
    int extras_;                /* accounts for extra buffers */
    int units_used_;            /* accounts for number of unit files used */
    int cur_unit_;              /* current unit file */
    int cur_size_;              /* current size of buffer */
    int first_unit_;            /* first file unit number (if > 1) */
    int subgr_per_irrep_;       /* possible number of Olsen subgraphs per irrep */
    int print_lvl_;             /* print level*/
    bool fopen_;                /* Are CIVec files open? */

    double ssq(struct stringwr *alplist, struct stringwr *betlist, double **CL,
               double **CR, int nas, int nbs, int Ja_list, int Jb_list);

   public:
    CIvect();
    CIvect(BIGINT vl, int nb, int incor, int ms0, int *iac, int *ibc, int *ias,
           int *ibs, BIGINT *offs, int nac, int nbc, int nirr, int cdperirr,
           int maxvect, int nunits, int funit, int *fablk, int *lablk, int **dc,
           struct calcinfo *CI_CalcInfo, struct params *CI_Params,
           struct H_zero_block *CI_H0block, bool buf_init = true);
    CIvect(int incor, int maxvect, int nunits, int funit,
           struct ci_blks *CIblks, struct calcinfo *CI_CalcInfo,
           struct params *CI_Params, struct H_zero_block *CI_H0block,
           bool buf_init = true);
    ~CIvect();

    /// Numpy interface to the current buffer
    py::buffer_info array_interface();

    /// BLAS equivalents for CIVectors
    void axpy(double a, SharedCIVector x, int tvec, int ovec);
    void scale(double a, int tvec);
    void shift(double a, int tvec);
    void copy(SharedCIVector src, int tvec, int ovec);
    void divide(SharedCIVector denom, double min_val, int tvec, int ovec);
    void zero(void);
    double vdot(SharedCIVector b, int tvec, int ovec);
    double norm(int tvec);

    // self += scale * a * b
    void vector_multiply(double scale, SharedCIVector X, SharedCIVector Y, int tvec, int xvec, int yvec);

    /// Specific CIVector operations
    double dcalc3(double lambda, SharedCIVector Hd, int rootnum);
    void symnormalize(double a, int tvec);

    /// Disk/memory manipulation
    void init_io_files(bool open_old);
    void close_io_files(int keep);
    int read(int tvec, int ibuf);
    int write(int tvec, int ibuf);
    void buf_lock(double *a);
    void buf_unlock(void);
    double *buf_malloc(void);
    void set_nvect(int i);

    // Questionable functions and/or should be private
    void set(int incor, int maxvect, int nunits, int funit,
             struct ci_blks *CIblks);
    void set(BIGINT vl, int nb, int incor, int ms0, int *iac, int *ibc,
             int *ias, int *ibs, BIGINT *offs, int nac, int nbc, int nirr,
             int cdperirr, int maxvect, int nunits, int funit, int *fablk,
             int *lablk, int **dc);
    void print();
    double operator*(CIvect &b);
    void setarray(const double *a, BIGINT len);
    void max_abs_vals(int nval, int *iac, int *ibc, int *iaidx, int *ibidx,
                      double *coeff, int neg_only);
    double blk_max_abs_vals(int i, int offdiag, int nval, int *iac, int *ibc,
                            int *iaidx, int *ibidx, double *coeff,
                            double minval, int neg_only);
    void det2strings(BIGINT det, int *alp_code, int *bet_code, int *alp_idx,
                     int *bet_idx);
    BIGINT strings2det(int alp_code, int alp_idx, int bet_code, int bet_idx);
    void diag_mat_els(struct stringwr **alplist, struct stringwr **betlist,
                      double *oei, double *tei, double edrc, int na, int nb,
                      int nbf, int method);
    void diag_mat_els_otf(struct stringwr **alplist, struct stringwr **betlist,
                          double *oei, double *tei, double edrc, int na, int nb,
                          int nbf, int buf, int method);
    void init_vals(int ivect, int nvals, int *alplist, int *alpidx,
                   int *betlist, int *betidx, int *blknums, double *value);
    void set_vals(int ivect, int nvals, int *alplist, int *alpidx, int *betlist,
                  int *betidx, int *blknums, double *value);
    void extract_vals(int ivect, int nvals, int *alplist, int *alpidx,
                      int *betlist, int *betidx, int *blknums, double *value);
    void symnorm(double a, int vecode, int gather_vec);
    double zero_det(int iac, int ia, int ibc, int ib);
    void scale(double a, int vecode, int gather_vec);
    void symmetrize(double phase, int iblock);
    double **blockptr(int blknum);
    int schmidt_add(CIvect &c, int L);
    int schmidt_add2(CIvect &c, int first_vec, int last_vec, int source_vec,
                     int target_vec, double *dotval, double *nrm,
                     double *ovlpmax);
    void dcalc(int nr, int L, double **alpha, double *lambda, double *norm_arr,
               CIvect &C, CIvect &S, double *buf1, double *buf2,
               int *root_converged, int printflag,
               double *E_est);
    void sigma_renorm(int nr, int L, double renorm_C, CIvect &S, double *buf1,
                      int printflag);
    double dcalc2(int rootnum, double lambda, CIvect &Hd, int precon,
                  struct stringwr **alplist, struct stringwr **betlist);
    double dcalc_evangelisti(int rootnum, int num_vecs, double lambda,
                             CIvect &Hd, CIvect &C, double *buf1, double *buf2,
                             int precon, int L, struct stringwr **alplist,
                             struct stringwr **betlist, double **alpha);
    void construct_kth_order_wf(CIvect &Hd, CIvect &S, CIvect &C,
                                struct stringwr **alplist,
                                struct stringwr **betlist, double *buf1,
                                double *buf2, int k, double *mp_energy,
                                double **bvec_overlap, double *bvec_norm);
    void wigner_E2k_formula(CIvect &Hd, CIvect &S, CIvect &C,
                            struct stringwr **alplist,
                            struct stringwr **betlist, double *buf1,
                            double *buf2, int k, double *mp2k_energy,
                            double **wfn_overlap, double **bvec_overlap,
                            double *bvec_norm, int kvec_offset);
    void print_buf();
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
    double calc_ssq(double *buffer1, double *buffer2, struct stringwr **alplist,
                    struct stringwr **betlist, int vec_num);
    void h0block_buf_init(void);
    void h0block_buf_ols(double *norm, double *ovrlap, double *c1norm,
                         double E_est);
    void h0block_buf_precon(double *norm, int root);
    void h0block_gather_vec(int vecode);
    void h0block_gather_multivec(double *vec);
    int check_zero_block(int blocknum);
    void set_zero_block(int blocknum, int value);
    void set_zero_blocks_all(void);
    void copy_zero_blocks(CIvect &src);
    void print_zero_blocks(void);
    void scale_sigma(CIvect &Hd, CIvect &C, struct stringwr **alplist,
                     struct stringwr **betlist, int i, double *buf1,
                     double *buf2);
    int read_new_first_buf(void);
    void write_new_first_buf(void);
    void set_new_first_buf(int nfb);
    int read_num_vecs(void);
    void write_num_vecs(int nv);
    void write_toc(void);
    void civect_psio_debug(void);
    void pt_correction(struct stringwr **alplist, struct stringwr **betlist);
    double compute_follow_overlap(int troot, int ncoef, double *coef, int *Iac,
                                  int *Iaridx, int *Ibc, int *Ibridx);

    void calc_hd_block(struct stringwr *alplist, struct stringwr *betlist,
                       double **H0, double *oei, double *tei, double edrc,
                       int nas, int nbs, int na, int nb, int nbf);
    void calc_hd_block_ave(struct stringwr *alplist, struct stringwr *betlist,
                           double **H0, double *tf_oei, double *tei,
                           double edrc, int nas, int nbs, int na, int nb,
                           int nbf);
    void calc_hd_block_z_ave(struct stringwr *alplist, struct stringwr *betlist,
                             double **H0, double pert_param, double *tei,
                             double edrc, int nas, int nbs, int na, int nb,
                             int nbf);
    void calc_hd_block_orbenergy(struct stringwr *alplist,
                                 struct stringwr *betlist, double **H0,
                                 double *oei, double *tei, double edrc, int nas,
                                 int nbs, int na, int nb, int nbf);
    void calc_hd_block_mll(struct stringwr *alplist, struct stringwr *betlist,
                           double **H0, double *oei, double *tei, double edrc,
                           int nas, int nbs, int na, int nb, int nbf);
    void calc_hd_block_evangelisti(struct stringwr **alplist,
                                   struct stringwr **betlist,
                                   struct stringwr *alplist_local,
                                   struct stringwr *betlist_local, double **H0,
                                   double *tf_oei, double *tei, double edrc,
                                   int nas, int nbs, int na, int nb, int nbf);
};
}
}  // namespace psi::detci

#endif  // header guard
