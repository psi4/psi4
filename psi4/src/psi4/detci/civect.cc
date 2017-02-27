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
**  \ingroup DETCI
**  \brief Code for the CI vector class
**
** David Sherrill, 15 June 1995
** Center for Comptuational Quantum Chemistry
**
** Update Feb 1996 to ensure everything works symmetry block at a time
** Rewrite blk_xxx members to buf_xxx to avoid confusion with RAS blocks
**
** Current working assumption for Ms=0: try to fix it so that no disk
** space is required for redundant buffers, but once in memory, assume
** that a buffer has been transposed to give redundant information.  That
** is, assume complete core storage for whole buffer, but don't write
** redundant buffers to disk.
** Modification: actually, don't store redundant _buffers_, but store
** redundant blocks if needed.
**
*/


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/pybind11.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/ci_tol.h"
#include "psi4/detci/civect.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "psi4/libmints/vector.h"

namespace psi { namespace detci {

extern void transp_sigma(double **a, int rows, int cols, int phase);
extern void xey(double *x, double *y, int size);
extern void xeay(double *x, double a, double *y, int size);
extern void xpeay(double *x, double a, double *y, int size);
extern void xpey(double *x, double *y, int size);
extern void xeax(double *x, double a, int size);
extern void xexmy(double *x, double *y, int size);
extern double calc_d2(double *target, double lambda, double *Hd,
   int size, int precon);
extern double calc_mpn_vec(double *target, double energy, double *Hd,
   int size, double sign1, double sign2, int precon);
extern void xeaxmy(double *x, double *y, double a, int size);
extern void xeaxpby(double *x, double *y, double a, double b, int size);
extern void xexy(double *x, double *y, int size);
extern int calc_orb_diff(int cnt, unsigned char *I, unsigned char *J,
   int *I_alpha_diff, int *J_alpha_diff, int *sign, int *same, int extended);


#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

CIvect::CIvect()  // Default constructor
{
    common_init();
}

CIvect::CIvect(BIGINT vl, int nb, int incor, int ms0, int *iac, int *ibc,
               int *ias, int *ibs, BIGINT *offs, int nac, int nbc, int nirr,
               int cdpirr, int mxv, int nu, int funit, int *fablk, int *lablk,
               int **dc, struct calcinfo *CI_CalcInfo, struct params *CI_Params,
               struct H_zero_block *CI_H0block, bool buf_init) {
    common_init();
    CI_CalcInfo_ = CI_CalcInfo;
    CI_Params_ = CI_Params;
    CI_H0block_ = CI_H0block;

    set(vl, nb, incor, ms0, iac, ibc, ias, ibs, offs, nac, nbc, nirr, cdpirr,
        mxv, nu, funit, fablk, lablk, dc);

    if (buf_init) {
        buffer_ = buf_malloc();
        blocks_[0][0] = buffer_;
        buf_lock(buffer_);
    }
}
CIvect::CIvect(int incor, int maxvect, int nunits, int funit,
               struct ci_blks *CIblks, struct calcinfo *CI_CalcInfo,
               struct params *CI_Params, struct H_zero_block *CI_H0block,
               bool buf_init) {
    common_init();
    CI_CalcInfo_ = CI_CalcInfo;
    CI_Params_ = CI_Params;
    CI_H0block_ = CI_H0block;

    set(CIblks->vectlen, CIblks->num_blocks, incor, CIblks->Ms0,
        CIblks->Ia_code, CIblks->Ib_code, CIblks->Ia_size, CIblks->Ib_size,
        CIblks->offset, CIblks->num_alp_codes, CIblks->num_bet_codes,
        CIblks->nirreps, CIblks->subgr_per_irrep, maxvect, nunits, funit,
        CIblks->first_iablk, CIblks->last_iablk, CIblks->decode);

    if (buf_init) {
        buffer_ = buf_malloc();
        blocks_[0][0] = buffer_;
        buf_lock(buffer_);
    }
}

void CIvect::common_init(void) {
    vectlen_ = 0;
    num_blocks_ = 0;
    icore_ = 1;
    Ms0_ = 0;
    num_alpcodes_ = 0;
    num_betcodes_ = 0;
    nirreps_ = 0;
    codes_per_irrep_ = 0;
    buf_per_vect_ = 0;
    buf_total_ = 0;
    new_first_buf_ = 0;
    maxvect_ = 0;
    nvect_ = 0;
    nunits_ = 0;
    cur_vect_ = -1;
    cur_buf_ = -1;
    buffer_size_ = 0;
    buf_size_ = NULL;
    buf2blk_ = NULL;
    buf_offdiag_ = NULL;
    first_ablk_ = NULL;
    last_ablk_ = NULL;
    decode_ = NULL;
    blocks_ = NULL;
    buf_locked_ = 0;
    buffer_ = NULL;
    in_file_ = 0;
    extras_ = 0;
    units_used_ = 0;
    cur_unit_ = 0;
    cur_size_ = 0;
    first_unit_ = 0;
    print_lvl_ = 0;
    fopen_ = false;
}

void CIvect::set(int incor, int maxvect, int nunits, int funit,
                 struct ci_blks *CIblks) {
    set(CIblks->vectlen, CIblks->num_blocks, incor, CIblks->Ms0,
        CIblks->Ia_code, CIblks->Ib_code, CIblks->Ia_size, CIblks->Ib_size,
        CIblks->offset, CIblks->num_alp_codes, CIblks->num_bet_codes,
        CIblks->nirreps, CIblks->subgr_per_irrep, maxvect, nunits, funit,
        CIblks->first_iablk, CIblks->last_iablk, CIblks->decode);
}

void CIvect::set(BIGINT vl, int nb, int incor, int ms0, int *iac, int *ibc,
                 int *ias, int *ibs, BIGINT *offs, int nac, int nbc, int nirr,
                 int cdpirr, int mxv, int nu, int fu, int *fablk, int *lablk,
                 int **dc) {
    int i, j, k;
    int maxrows = 0, maxcols = 0;
    unsigned long bufsize, maxbufsize;
    unsigned long size;
    // static int first=1;
    /* int in_file, extras, units_used, cur_unit; */

    vectlen_ = vl;
    num_blocks_ = nb;
    icore_ = incor;
    Ms0_ = ms0;
    nirreps_ = nirr;
    codes_per_irrep_ = cdpirr;
    maxvect_ = mxv;
    nvect_ = 1;
    nunits_ = nu;
    if (nunits_) units_.resize(nunits_);
    first_unit_ = fu;
    for (i = 0; i < nunits_; i++) units_[i] = fu + i;

    Ia_code_.resize(nb);
    Ib_code_.resize(nb);
    Ia_size_.resize(nb);
    Ib_size_.resize(nb);
    offset_.resize(nb);
    // offset_ = (BIGINT *)malloc(nb * sizeof(BIGINT));

    for (i = 0; i < nb; i++) {
        Ia_code_[i] = iac[i];
        Ib_code_[i] = ibc[i];
        Ia_size_[i] = ias[i];
        Ib_size_[i] = ibs[i];
        offset_[i] = offs[i];
    }

    num_alpcodes_ = nac;
    num_betcodes_ = nbc;

    first_ablk_ = init_int_array(nirr);
    last_ablk_ = init_int_array(nirr);
    for (i = 0; i < nirr; i++) {
        first_ablk_[i] = fablk[i];
        last_ablk_[i] = lablk[i];
    }

    decode_ = init_int_matrix(nac, nbc);
    for (i = 0; i < nac; i++) {
        for (j = 0; j < nbc; j++) {
            decode_[i][j] = dc[i][j];
        }
    }

    if (icore_ == 1) { /* whole vector in-core */
        buf_per_vect_ = 1;
        buf_total_ = maxvect_;
        buf_size_ = (unsigned long *)malloc(buf_per_vect_ * sizeof(unsigned long));
        buf2blk_ = init_int_array(buf_per_vect_);
        buf_offdiag_ = init_int_array(buf_per_vect_);
        for (i = 0; i < buf_per_vect_; i++) buf_size_[i] = 0;
        if (maxvect_ < nunits_) nunits_ = maxvect_;
        if (nvect_ > maxvect_) nvect_ = maxvect_;
        size = vectlen_; /* may want to change for Ms=0 later */
        for (i = 0; i < buf_per_vect_; i++) {
            buf2blk_[i] = -1;
            buf_size_[i] = size;
        }
    }

    if (icore_ == 2) { /* whole symmetry block in-core */
        /* figure out how many buffers per vector */
        for (i = 0, buf_per_vect_ = 0; i < nirreps_; i++) {
            j = first_ablk_[i];
            if (j < 0) continue;
            if (!Ms0_ || CI_CalcInfo_->ref_sym == 0)
                buf_per_vect_++;
            else if (Ia_code_[j] / codes_per_irrep_ >
                     Ib_code_[j] / codes_per_irrep_)
                buf_per_vect_++;
        }

        buf2blk_ = init_int_array(buf_per_vect_);
        buf_offdiag_ = init_int_array(buf_per_vect_);

        for (i = 0, j = 0; i < nirreps_; i++) {
            k = first_ablk_[i];
            if (k < 0) continue;
            if (!Ms0_ || CI_CalcInfo_->ref_sym == 0) {
                buf2blk_[j] = i;
                j++;
            } else if (Ia_code_[k] / codes_per_irrep_ >
                       Ib_code_[k] / codes_per_irrep_) {
                buf_offdiag_[j] = 1;
                buf2blk_[j] = i;
                j++;
            }
        }

        buf_total_ = maxvect_ * buf_per_vect_;
        buf_size_ = (unsigned long *)malloc(buf_per_vect_ * sizeof(unsigned long));
        for (i = 0; i < buf_per_vect_; i++) {
            buf_size_[i] = 0;
            j = buf2blk_[i];
            for (k = first_ablk_[j]; k <= last_ablk_[j]; k++) {
                buf_size_[i] += (unsigned long)Ia_size_[k] * (unsigned long)Ib_size_[k];
            }
        }
    } /* end icore==2 */

    if (icore_ == 0) { /* one subblock in-core */

        /* figure out how many blocks are stored */
        buf_per_vect_ = 0;
        if (Ms0_) {
            for (i = 0; i < num_blocks_; i++) {
                if (Ia_code_[i] >= Ib_code_[i] && Ia_size_[i] > 0 &&
                    Ib_size_[i] > 0)
                    buf_per_vect_++;
            }
        } else {
            for (i = 0; i < num_blocks_; i++) {
                if (Ia_size_[i] > 0 && Ib_size_[i] > 0) buf_per_vect_++;
            }
        }

        buf_total_ = buf_per_vect_ * maxvect_;
        buf2blk_ = init_int_array(buf_per_vect_);
        buf_offdiag_ = init_int_array(buf_per_vect_);
        buf_size_ = (unsigned long *)malloc(buf_per_vect_ * sizeof(unsigned long));

        if (Ms0_) {
            for (i = 0, j = 0; i < num_blocks_; i++) {
                if (Ia_code_[i] >= Ib_code_[i] && Ia_size_[i] > 0 &&
                    Ib_size_[i] > 0) {
                    buf2blk_[j] = i;
                    buf_size_[j] = (unsigned long)Ia_size_[i] * (unsigned long)Ib_size_[i];
                    if (Ia_code_[i] != Ib_code_[i]) buf_offdiag_[j] = 1;
                    j++;
                }
            }
        } else {
            for (i = 0, j = 0; i < num_blocks_; i++) {
                if (Ia_size_[i] > 0 && Ib_size_[i] > 0) {
                    buf2blk_[j] = i;
                    buf_size_[j] = (unsigned long)Ia_size_[i] * (unsigned long)Ib_size_[i];
                    j++;
                }
            }
        }

    } /* end icore==0 */

    file_number_.resize(buf_total_);

    if (nunits_) {
        in_file_ = 0;
        extras_ = buf_total_ % nunits_;
        units_used_ = 0;
        cur_unit_ = units_[0];

        for (i = 0; i < buf_total_; i++) {
            if (in_file_ + 1 <= buf_total_ / nunits_) {
                file_number_[i] = cur_unit_;
                in_file_++;
            }

            else if ((in_file_ == buf_total_ / nunits_) && extras_) {
                file_number_[i] = cur_unit_;
                extras_--;
                in_file_++;
            }

            else {
                units_used_++;
                cur_unit_ = units_[units_used_];
                file_number_[i] = cur_unit_;
                in_file_ = 1;
            }
        } /* end loop over buffers */
    }

    // do next step separately now to control OPEN_NEW vs OPEN_OLD
    // init_io_files();

    /*
       outfile->Printf("num_blocks_ = %d\n", num_blocks_);
       for (i=0; i<buf_total_; i++)
          outfile->Printf("file_offset_[%d] = %lu\n ", i, file_offset_[i]);
       for (i=0; i<maxvect_; i++)
          outfile->Printf("zero_block_offset_[%d] = %lu\n ", i,zero_block_offset_[i]);
       for (i=0; i<buf_per_vect_; i++)
          outfile->Printf("buf_size_[%d] = %lu\n ", i, buf_size_[i]);
    */

    // Set up the flags for zero blocks...at first, all blocks are all 0's
    // but we will put '0' indicating that the program should assume
    // nonzero.  It is hard to put this stuff in correctly at this point,
    // and dangerous to throw out any nonzero data, so we will only set
    // the blocks to all zero when we know we want them to be treated as
    // all 0's for certain.
    zero_blocks_.resize(num_blocks_);

    // Figure out the buffer size and allocate some pointers (but not storage)
    blocks_ = (double ***)malloc(num_blocks_ * sizeof(double **));

    if (icore_ == 1) { /* everything is in-core */
        for (i = 0; i < num_blocks_; i++) {
            if (Ia_size_[i])
                blocks_[i] = (double **)malloc(Ia_size_[i] * sizeof(double *));
            else
                blocks_[i] = (double **)malloc(sizeof(double *));
        }
        buffer_size_ = vectlen_;
    }                       /* end icore==1 */
    else if (icore_ == 2) { /* one symmetry block is held at a time */
                            /* figure out which symmetry block is largest */
        for (i = 0, maxbufsize = 0; i < nirreps_; i++) {
            for (j = first_ablk_[i], bufsize = 0; j <= last_ablk_[i]; j++) {
                bufsize +=
                    (unsigned long)Ia_size_[j] * (unsigned long)Ib_size_[j];
            }
            if (bufsize > maxbufsize) maxbufsize = bufsize;
        }
        for (i = 0; i < num_blocks_; i++) {
            if (Ia_size_[i])
                blocks_[i] = (double **)malloc(Ia_size_[i] * sizeof(double *));
            else
                blocks_[i] = (double **)malloc(sizeof(double *));
        }
        buffer_size_ = maxbufsize;
    }                       /* end icore==2 */
    else if (icore_ == 0) { /* only one subblock in core at once */
        for (i = 0, maxbufsize = 0; i < num_blocks_; i++) {
            if (Ia_size_[i] > maxrows) maxrows = Ia_size_[i];
            if (Ib_size_[i] > maxcols) maxcols = Ib_size_[i];
            if (Ia_size_[i])
                blocks_[i] = (double **)malloc(Ia_size_[i] * sizeof(double *));
            else
                blocks_[i] = (double **)malloc(sizeof(double *));
            bufsize = (unsigned long)Ia_size_[i] * (unsigned long)Ib_size_[i];
            if (bufsize > maxbufsize) maxbufsize = bufsize;
        }
        // CDS 11/5/97: Revise buffer_size, the size of the biggest buffer
        // buffer_size = maxrows * maxcols;   Made buffers too large
        buffer_size_ = maxbufsize;  // Why didn't I do it this way before?
    }                               /* end icore==0 */
    else {
       outfile->Printf("CIvect::set(): unrecognized option for icore = %d\n", icore_);
        return;
    }

    // MLL 5/7/98: Want to know the subblock length of a vector //
    //  if (print_lvl_) {
    //     outfile->Printf("\n CI vector/subblock length = %ld\n", buffer_size);
    //
    //     }
}

CIvect::~CIvect() {
    if (num_blocks_) {
        if (buf_locked_) free(buffer_);
        for (int i = 0; i < num_blocks_; i++) {
            free(blocks_[i]);
        }

        free(blocks_);
        // free(zero_blocks_);
        // free(Ia_code_);
        // free(Ib_code_);
        // free(Ia_size_);
        // free(Ib_size_);
        // free(units_);
        // free(file_number_);
        free(buf_size_);
        free(buf2blk_);
        free(buf_offdiag_);
        free(first_ablk_);
        free(last_ablk_);
        free_int_matrix(decode_);
        // free(offset_);
    }
}

/*
** CIvect::buf_malloc
**
** Function malloc's memory for a buffer of the CI vector.  A buffer
** is the appropriate size to hold either one entire CI vector,
** a symmetry block of the CI vector, or a subblock of the vector,
** depending on the value of icore.
** Each subblock is a matrix, but the matrices are allocated as
** contiguous blocks of memory.
**
** Parameters:
**    none
**
** Returns:
**    pointer to memory buffer (double *)
*/
double *CIvect::buf_malloc(void) {
    double *tmp;

    tmp = init_array(buffer_size_);

    return (tmp);
}

/*
** CIvect::set_nvect
**
** Sets the value of nvect.  Usually set in CIvect::set().  Only need
** to call if it must be changed in an unusual way (i.e. resetting
** the Davidson subspace).
*/
void CIvect::set_nvect(int i) { nvect_ = i; }
/*
** CIvect::vdot
**
** Function returns the scalar product of two CI vectors.
** Assumes that diagonal blocks are full for Ms=0 cases
** tvec and ovec are this and other vector numbers, respectively.
*/

double CIvect::vdot(SharedCIVector b, int tvec, int ovec) {
    tvec = (tvec == -1) ? cur_vect_ : tvec;
    ovec = (ovec == -1) ? b->cur_vect_ : ovec;

    double dotprod = 0.0;
    if (Ms0_) {
        for (int buf = 0; buf < buf_per_vect_; buf++) {
            read(tvec, buf);
            b->read(ovec, buf);
            double tval = C_DDOT(buf_size_[buf], buffer_, 1, b->buffer_, 1);
            // dot_arr(buffer_, b->buffer_, buf_size_[buf], &tval);
            if (buf_offdiag_[buf]) tval *= 2.0;
            dotprod += tval;
        }
    }

    else {
        for (int buf = 0; buf < buf_per_vect_; buf++) {
            read(tvec, buf);
            b->read(ovec, buf);
            double tval = C_DDOT(buf_size_[buf], buffer_, 1, b->buffer_, 1);
            // dot_arr(buffer_, b->buffer_, buf_size_[buf], &tval);
            dotprod += tval;
        }
    }

    return dotprod;
}
void CIvect::axpy(double a, SharedCIVector X, int yvec, int xvec) {
    // Y += a * X
    for (int buf = 0; buf < buf_per_vect_; buf++) {
        X->read(xvec, buf);
        read(yvec, buf);
        C_DAXPY(buf_size_[buf], a, X->buffer_, 1, buffer_, 1);
        write(yvec, buf);
    }
}

void CIvect::scale(double a, int vec) {
    // this *= a
    for (int buf = 0; buf < buf_per_vect_; buf++) {
        read(vec, buf);
        C_DSCAL(buf_size_[buf], a, buffer_, 1);
        write(vec, buf);
    }
}
void CIvect::shift(double a, int vec) {
    // this *= a
    for (int buf = 0; buf < buf_per_vect_; buf++) {
        read(vec, buf);
        for (size_t i = 0; i < buf_size_[buf]; ++i) {
            buffer_[i] += a;
        }
        write(vec, buf);
    }
}

void CIvect::copy(SharedCIVector src, int tvec, int ovec) {
    for (int buf = 0; buf < buf_per_vect_; buf++) {
        src->read(ovec, buf);
        read(tvec, buf);
        C_DCOPY(buf_size_[buf], src->buffer_, 1, buffer_, 1);
        int blk = buf2blk_[buf];
        if ((blk >= 0) && ((zero_blocks_[blk] == 0) || (src->zero_blocks_[blk] == 0)))
            zero_blocks_[blk] = 0;
        write(tvec, buf);
    }
}
void CIvect::divide(SharedCIVector denom, double min_val, int tvec, int ovec) {
   for (int buf=0; buf<buf_per_vect_; buf++) {
      denom->read(ovec, buf);
      read(tvec, buf);
      for (size_t i = 0; i < buf_size_[buf]; ++i) {
          if (std::fabs(denom->buffer_[i]) > min_val) {
              buffer_[i] /= denom->buffer_[i];
          }
          else{
            buffer_[i] = 0;
          }
      }
      write(tvec, buf);
  }
}

py::buffer_info CIvect::array_interface() {

    // Why is this so complex other places?
    if (!buf_locked_)
        throw PSIEXCEPTION("CIVector::matrix_array_interface: No buffer is locked.");

    return py::buffer_info(buffer_, sizeof(double),
                           py::format_descriptor<double>::format(), 1,
                           {static_cast<size_t>(buffer_size_)},
                           {sizeof(double)});
}

double CIvect::dcalc3(double lambda, SharedCIVector Hd, int rootnum) {
    double normsq = 0.0;
    for (int buf = 0; buf < buf_per_vect_; buf++) {
        read(rootnum, buf);
        Hd->read(0, buf);

        double tval = 0.0;
        for (size_t i = 0; i < (size_t)buf_size_[buf]; i++) {
            double divisor = (lambda - Hd->buffer_[i]);
            if (std::fabs(divisor) > HD_MIN) {
                buffer_[i] /= divisor;
                tval += buffer_[i] * buffer_[i];
            } else {
                buffer_[i] = 0.0;
            }
        }

        if (buf_offdiag_[buf]) tval *= 2.0;
        normsq += tval;
        write(rootnum, buf);
    }
    double norm = std::sqrt(normsq);
    return (norm);
}

void CIvect::symnormalize(double a, int tvec) {
    int i, j;
    int blk, ac, bc, upper;
    double **mat, *arr, phase;

    if (!Ms0_) {
        scale(a, tvec);
        return;
    }

    if (!CI_Params_->Ms0)
        phase = 1.0;
    else
        phase = ((int)CI_Params_->S % 2) ? -1.0 : 1.0;

    if (icore_ == 1) {
        read(tvec, 0);
        for (blk = 0; blk < num_blocks_; blk++) {
            ac = Ia_code_[blk];
            bc = Ib_code_[blk];
            mat = blocks_[blk];
            if (ac == bc) { /* diagonal block */
                for (i = 0; i < Ia_size_[blk]; i++) {
                    mat[i][i] *= a;
                    for (j = 0; j < i; j++) {
                        mat[i][j] *= a;
                        mat[j][i] = mat[i][j] * phase;
                    }
                }
            }

            if (ac > bc) { /* off-diagonal block */
                // xeax(blocks_[blk][0], a, Ia_size_[blk] * Ib_size_[blk]);
                C_DSCAL(Ia_size_[blk] * Ib_size_[blk], a, blocks_[blk][0], 1);
                upper = decode_[bc][ac];
                if (upper >= 0) {
                    zero_blocks_[upper] = zero_blocks_[blk];
                    for (i = 0; i < Ia_size_[blk]; i++) {
                        for (j = 0; j < Ib_size_[blk]; j++) {
                            blocks_[upper][j][i] = mat[i][j] * phase;
                        }
                    }
                }
            }
        } /* end loop over blocks */

        write(tvec, 0);

    } /* end icore_ == 1 */

    else {
       outfile->Printf("(CIvect::symnorm): Only supports incore=1 at the moment\n");
        return;
    }
}

double CIvect::norm(int tvec) {
    tvec = (tvec == -1) ? cur_vect_ : tvec;

    double dotprod = 0.0;
    if (Ms0_) {
        for (int buf = 0; buf < buf_per_vect_; buf++) {
            read(tvec, buf);
            double tval = C_DDOT(buf_size_[buf], buffer_, 1, buffer_, 1);
            if (buf_offdiag_[buf]) tval *= 2.0;
            dotprod += tval;
        }
    }

    else {
        for (int buf = 0; buf < buf_per_vect_; buf++) {
            read(tvec, buf);
            double tval = C_DDOT(buf_size_[buf], buffer_, 1, buffer_, 1);
            dotprod += tval;
        }
    }

    return std::sqrt(dotprod);
}

void CIvect::vector_multiply(double scale, SharedCIVector X, SharedCIVector Y, int tvec, int xvec, int yvec) {
    // T += a * X * Y
    for (int buf = 0; buf < buf_per_vect_; buf++) {
        X->read(xvec, buf);
        Y->read(yvec, buf);
        read(tvec, buf);
        for (size_t i=0; i<buf_size_[buf]; i++){
            buffer_[i] += scale * X->buffer_[i] * Y->buffer_[i];
        }
        write(tvec, buf);
    }
}

/*
** CIvect::operator *
**
** Function returns the scalar product of two CI vectors.
** Assumes that diagonal blocks are full for Ms=0 cases
*/

double CIvect::operator*(CIvect &b)
{
   double dotprod=0.0, tval;
   int i, buf;

   if (Ms0_) {
      for (buf=0; buf<buf_per_vect_; buf++) {
         read(cur_vect_, buf);
         b.read(b.cur_vect_, buf);
         dot_arr(buffer_, b.buffer_, buf_size_[buf], &tval);
         if (buf_offdiag_[buf]) tval *= 2.0;
         dotprod += tval;
         }
      }

   else {
      for (buf=0; buf<buf_per_vect_; buf++) {
         read(cur_vect_, buf);
         b.read(b.cur_vect_, buf);
         dot_arr(buffer_, b.buffer_, buf_size_[buf], &tval);
         dotprod += tval;
         }
      }

   return(dotprod);
}

void CIvect::setarray(const double *a, BIGINT len) {
    double *aptr;
    BIGINT i;

    if (len > vectlen_) len = vectlen_;

    if (icore_ == 1) {
        aptr = buffer_;
        for (i = 0; i < len; i++) {
            aptr[i] = a[i];
        }
    }

    else {
        outfile->Printf("(CIvect::setarray): Invalid icore option!\n");
        outfile->Printf("   use only for icore_=1\n");
    }
}

void CIvect::max_abs_vals(int nval, int *iac, int *ibc, int *iaidx, int *ibidx,
      double *coeff, int neg_only)
{
   int i,buf,irrep;
   double minval=0.0;


   if (icore_ == 1) {
      for (i=0; i<num_blocks_; i++) {
         minval = blk_max_abs_vals(i, 0, nval, iac, ibc, iaidx, ibidx, coeff,
                    minval, neg_only);
         }
      } /* end case icore==1 */

   if (icore_ == 2) { /* symmetry block at a time */
      for (buf=0; buf<buf_per_vect_; buf++) {
         if (!read(cur_vect_, buf)) continue;
         irrep = buf2blk_[buf];
         for (i=first_ablk_[irrep]; i<=last_ablk_[irrep]; i++) {
            minval = blk_max_abs_vals(i, buf_offdiag_[buf], nval, iac, ibc,
               iaidx, ibidx, coeff, minval, neg_only);
            }
         }
      } /* end case icore==2 */

   if (icore_ == 0) { /* RAS block at a time */
      for (buf=0; buf<buf_per_vect_; buf++) {
         if (!read(cur_vect_, buf)) continue;
         i = buf2blk_[buf];
         minval = blk_max_abs_vals(i, buf_offdiag_[buf], nval, iac, ibc,
            iaidx, ibidx, coeff, minval, neg_only);
         }
      } /* end case icore==0 */

}


double CIvect::blk_max_abs_vals(int i, int offdiag, int nval, int *iac,
       int *ibc, int *iaidx, int *ibidx, double *coeff, double minval,
        int neg_only)
{
   int j,k,m,n;
   double value, abs_value;
   int iacode, ibcode;

   iacode = Ia_code_[i];
   ibcode = Ib_code_[i];
   for (j=0; j<Ia_size_[i]; j++) {
      for (k=0; k<Ib_size_[i]; k++) {
         value = blocks_[i][j][k];
         if ((value > 0.0) && (neg_only)) continue;
         abs_value = fabs(value);
         if (abs_value >= fabs(minval)) {
            for (m=0; m<nval; m++) {
               if (abs_value > fabs(coeff[m])) {
                  for (n=nval-1; n>m; n--) {
                     coeff[n] = coeff[n-1];
                     iac[n] = iac[n-1];
                     ibc[n] = ibc[n-1];
                     iaidx[n] = iaidx[n-1];
                     ibidx[n] = ibidx[n-1];
                     }
                  coeff[n] = value;
                  iac[n] = iacode;
                  ibc[n] = ibcode;
                  iaidx[n] = j;
                  ibidx[n] = k;
                  break;
                  }
               }
            CI_H0block_->spin_cp_vals = minval;
            minval = coeff[nval-1];
            }
         if (offdiag) {
            if (CI_Params_->Ms0 && ((int) CI_Params_->S % 2) &&
                (!neg_only)) value = -value;
            if (abs_value >= minval) {
               for (m=0; m<nval; m++) {
                  if (abs_value > fabs(coeff[m])) {
                     for (n=nval-1; n>m; n--) {
                        coeff[n] = coeff[n-1];
                        iac[n] = iac[n-1];
                        ibc[n] = ibc[n-1];
                        iaidx[n] = iaidx[n-1];
                        ibidx[n] = ibidx[n-1];
                        }
                     coeff[n] = value;
                     iac[n] = ibcode;
                     ibc[n] = iacode;
                     iaidx[n] = k;
                     ibidx[n] = j;
                     break;
                     }
                  }
               CI_H0block_->spin_cp_vals = minval;
               minval = coeff[nval-1];
               }
            }
         }
      }

 /*
   for (i=0; i<CI_H0block_->size+CI_H0block_->coupling_size; i++)
      outfile->Printf("CI_H0block_->H00[%d] = %lf\n",i,CI_H0block_->H00[i]);
   outfile->Printf("CI_H0block_->spin_cp_vals = %lf\n",CI_H0block_->spin_cp_vals);
   outfile->Printf("minval = %lf\n",minval);
   outfile->Printf("printed in civect::blk_max_abs_vals\n");
 */
   return(minval);
}


void CIvect::det2strings(BIGINT det, int *alp_code, int *alp_idx,
         int *bet_code, int *bet_idx)
{
   int i;

   /* determine the CI block we're in */
   for (i=0; i<num_blocks_-1; i++) {
      if (offset_[i+1] > det) break;
      }
   *alp_code = Ia_code_[i];
   *bet_code = Ib_code_[i];

   *alp_idx = (int) ((det - offset_[i]) / (BIGINT) Ib_size_[i]);
   *bet_idx = ((det - offset_[i]) % (BIGINT) Ib_size_[i]);

}

BIGINT CIvect::strings2det(int alp_code, int alp_idx,
      int bet_code, int bet_idx)
{
   int blknum;
   BIGINT addr;

   blknum = decode_[alp_code][bet_code];
   addr = offset_[blknum];
   addr += alp_idx * Ib_size_[blknum] + bet_idx;

   return(addr);
}


void CIvect::diag_mat_els(struct stringwr **alplist, struct stringwr
      **betlist, double *oei, double *tei, double edrc, int na, int nb,
      int nbf, int method)
{

   int block, buf, iac, ibc, ias, ibs, irrep;
   double minval=0.0;

   if (icore_ == 1) { /* whole vector in-core */
      for (block=0; block<num_blocks_; block++) {
         iac = Ia_code_[block];
         ibc = Ib_code_[block];
         ias = Ia_size_[block];
         ibs = Ib_size_[block];
         if (method == HD_KAVE)
           calc_hd_block_ave(alplist[iac], betlist[ibc], blocks_[block],
              oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == ORB_ENER)
           calc_hd_block_orbenergy(alplist[iac], betlist[ibc], blocks_[block],
              oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == EVANGELISTI)
           calc_hd_block_evangelisti(alplist, betlist, alplist[iac], betlist[ibc], blocks_[block],
              oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == LEININGER)
           calc_hd_block_mll(alplist[iac], betlist[ibc], blocks_[block],
              oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == HD_EXACT)
           calc_hd_block(alplist[iac], betlist[ibc], blocks_[block],
               oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == Z_HD_KAVE)
           calc_hd_block_z_ave(alplist[iac], betlist[ibc], blocks_[block],
              CI_Params_->perturbation_parameter, tei, edrc, ias, ibs, na,
              nb, nbf);
         else {
           throw PsiException("hd_ave option not recognized.",__FILE__,__LINE__);
           }

         if (CI_Params_->hd_otf && CI_H0block_->size) {
            minval = blk_max_abs_vals(block, 0,
                    (CI_H0block_->size+CI_H0block_->coupling_size), CI_H0block_->alplist,
                    CI_H0block_->betlist, CI_H0block_->alpidx, CI_H0block_->betidx,
                    CI_H0block_->H00, minval, CI_Params_->neg_only);
           }

         }
      if (!CI_Params_->hd_otf) write(0,0);

      /*
      outfile->Printf("Diagonal matrix elements\n");
      print(outfile);
      */
      } /* end icore==1 */

   else if (icore_ == 2) { /* whole symmetry block at a time */
      for (buf=0; buf<buf_per_vect_; buf++) {
         irrep = buf2blk_[buf];
         for (block=first_ablk_[irrep]; block<=last_ablk_[irrep]; block++) {
            iac = Ia_code_[block];
            ibc = Ib_code_[block];
            ias = Ia_size_[block];
            ibs = Ib_size_[block];
            if (method == HD_KAVE)
              calc_hd_block_ave(alplist[iac], betlist[ibc], blocks_[block],
                 oei, tei, edrc, ias, ibs, na, nb, nbf);
            else if (method == ORB_ENER)
              calc_hd_block_orbenergy(alplist[iac], betlist[ibc],
                 blocks_[block], oei, tei, edrc, ias, ibs, na, nb, nbf);
            else if (method == EVANGELISTI)
              calc_hd_block_evangelisti(alplist, betlist, alplist[iac], betlist[ibc],
                 blocks_[block], oei, tei, edrc, ias, ibs, na, nb, nbf);
            else if (method == LEININGER)
              calc_hd_block_mll(alplist[iac], betlist[ibc],
                 blocks_[block], oei, tei, edrc, ias, ibs, na, nb, nbf);
            else if (method == HD_EXACT)
              calc_hd_block(alplist[iac], betlist[ibc], blocks_[block],
                 oei, tei, edrc, ias, ibs, na, nb, nbf);
            else if (method == Z_HD_KAVE)
              calc_hd_block_z_ave(alplist[iac], betlist[ibc], blocks_[block],
              CI_Params_->perturbation_parameter, tei, edrc, ias, ibs, na,
              nb, nbf);
            else {
              throw PsiException("hd_ave option not recognized.",__FILE__,__LINE__);
              }

            if (CI_Params_->hd_otf && CI_H0block_->size) {
              minval = blk_max_abs_vals(block, buf_offdiag_[buf],
                    (CI_H0block_->size+CI_H0block_->coupling_size),
                    CI_H0block_->alplist, CI_H0block_->betlist,
                    CI_H0block_->alpidx, CI_H0block_->betidx, CI_H0block_->H00,
                    minval, CI_Params_->neg_only);
              }

            }
         if (!CI_Params_->hd_otf) write(0,buf);
         }

      } /* end icore==2 */

   else if (icore_ == 0) { /* one subblock at a time */
      for (buf=0; buf<buf_per_vect_; buf++) {
         block = buf2blk_[buf];
         iac = Ia_code_[block];
         ibc = Ib_code_[block];
         ias = Ia_size_[block];
         ibs = Ib_size_[block];
        if (method == HD_KAVE)
          calc_hd_block_ave(alplist[iac], betlist[ibc], blocks_[block],
             oei, tei, edrc, ias, ibs, na, nb, nbf);
        else if (method == ORB_ENER)
          calc_hd_block_orbenergy(alplist[iac], betlist[ibc], blocks_[block],
             oei, tei, edrc, ias, ibs, na, nb, nbf);
        else if (method == EVANGELISTI)
          calc_hd_block_evangelisti(alplist, betlist, alplist[iac], betlist[ibc], blocks_[block],
             oei, tei, edrc, ias, ibs, na, nb, nbf);
        else if (method == LEININGER)
          calc_hd_block_mll(alplist[iac], betlist[ibc], blocks_[block],
             oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == HD_EXACT)
           calc_hd_block(alplist[iac], betlist[ibc], blocks_[block],
              oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == Z_HD_KAVE)
           calc_hd_block_z_ave(alplist[iac], betlist[ibc], blocks_[block],
              CI_Params_->perturbation_parameter, tei, edrc, ias, ibs, na,
              nb, nbf);
         else {
           throw PsiException("hd_ave option not recognized.",__FILE__,__LINE__);
           }
         if (CI_Params_->hd_otf && CI_H0block_->size) {
            minval = blk_max_abs_vals(block, buf_offdiag_[buf],
                    (CI_H0block_->size+CI_H0block_->coupling_size),
                    CI_H0block_->alplist, CI_H0block_->betlist, CI_H0block_->alpidx,
                    CI_H0block_->betidx, CI_H0block_->H00, minval, CI_Params_->neg_only);
           }
         if (!CI_Params_->hd_otf) write(0,buf);
         }
      } /* end icore==0 */

   else {
     outfile->Printf("(diag_mat_els): Unrecognized icore_ option!\n");
      }
}


void CIvect::diag_mat_els_otf(struct stringwr **alplist, struct stringwr
      **betlist, double *oei, double *tei, double edrc, int na, int nb,
      int nbf, int buf, int method)
{

   int block, iac, ibc, ias, ibs, irrep;

   if (icore_ == 1) { /* whole vector in-core */
      for (block=0; block<num_blocks_; block++) {
         iac = Ia_code_[block];
         ibc = Ib_code_[block];
         ias = Ia_size_[block];
         ibs = Ib_size_[block];
         if (method == HD_KAVE)
           calc_hd_block_ave(alplist[iac], betlist[ibc], blocks_[block],
              oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == ORB_ENER)
           calc_hd_block_orbenergy(alplist[iac], betlist[ibc], blocks_[block],
              oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == EVANGELISTI)
           calc_hd_block_evangelisti(alplist, betlist, alplist[iac], betlist[ibc], blocks_[block],
              oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == LEININGER)
           calc_hd_block_mll(alplist[iac], betlist[ibc], blocks_[block],
              oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == HD_EXACT)
           calc_hd_block(alplist[iac], betlist[ibc], blocks_[block],
               oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == Z_HD_KAVE)
           calc_hd_block_z_ave(alplist[iac], betlist[ibc], blocks_[block],
              CI_Params_->perturbation_parameter, tei, edrc, ias, ibs, na,
              nb, nbf);
         else {
           throw PsiException("hd_ave option not recognized.",__FILE__,__LINE__);
           }
         }
      } /* end icore==1 */

   else if (icore_ == 2) { /* whole symmetry block at a time */
         irrep = buf2blk_[buf];
         for (block=first_ablk_[irrep]; block<=last_ablk_[irrep]; block++) {
            iac = Ia_code_[block];
            ibc = Ib_code_[block];
            ias = Ia_size_[block];
            ibs = Ib_size_[block];
            if (method == HD_KAVE)
              calc_hd_block_ave(alplist[iac], betlist[ibc], blocks_[block],
                 oei, tei, edrc, ias, ibs, na, nb, nbf);
            else if (method == ORB_ENER)
              calc_hd_block_orbenergy(alplist[iac], betlist[ibc],
                 blocks_[block], oei, tei, edrc, ias, ibs, na, nb, nbf);
            else if (method == EVANGELISTI)
              calc_hd_block_evangelisti(alplist, betlist, alplist[iac], betlist[ibc],
                 blocks_[block], oei, tei, edrc, ias, ibs, na, nb, nbf);
            else if (method == LEININGER)
              calc_hd_block_mll(alplist[iac], betlist[ibc],
                 blocks_[block], oei, tei, edrc, ias, ibs, na, nb, nbf);
            else if (method == HD_EXACT)
              calc_hd_block(alplist[iac], betlist[ibc], blocks_[block],
                 oei, tei, edrc, ias, ibs, na, nb, nbf);
            else if (method == Z_HD_KAVE)
              calc_hd_block_z_ave(alplist[iac], betlist[ibc], blocks_[block],
              CI_Params_->perturbation_parameter, tei, edrc, ias, ibs, na,
              nb, nbf);
            else {
              throw PsiException("hd_ave option not recognized.",__FILE__,__LINE__);
              }
            }
      } /* end icore==2 */

   else if (icore_ == 0) { /* one subblock at a time */
         block = buf2blk_[buf];
         iac = Ia_code_[block];
         ibc = Ib_code_[block];
         ias = Ia_size_[block];
         ibs = Ib_size_[block];
        if (method == HD_KAVE)
          calc_hd_block_ave(alplist[iac], betlist[ibc], blocks_[block],
             oei, tei, edrc, ias, ibs, na, nb, nbf);
        else if (method == ORB_ENER)
          calc_hd_block_orbenergy(alplist[iac], betlist[ibc], blocks_[block],
             oei, tei, edrc, ias, ibs, na, nb, nbf);
        else if (method == EVANGELISTI)
          calc_hd_block_evangelisti(alplist, betlist, alplist[iac], betlist[ibc], blocks_[block],
             oei, tei, edrc, ias, ibs, na, nb, nbf);
        else if (method == LEININGER)
          calc_hd_block_mll(alplist[iac], betlist[ibc], blocks_[block],
             oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == HD_EXACT)
           calc_hd_block(alplist[iac], betlist[ibc], blocks_[block],
              oei, tei, edrc, ias, ibs, na, nb, nbf);
         else if (method == Z_HD_KAVE)
           calc_hd_block_z_ave(alplist[iac], betlist[ibc], blocks_[block],
              CI_Params_->perturbation_parameter, tei, edrc, ias, ibs, na,
              nb, nbf);
         else {
           throw PsiException("hd_ave option not recognized.",__FILE__,__LINE__);
           }
      } /* end icore==0 */

   else {
     outfile->Printf("(diag_mat_els): Unrecognized icore_ option!\n");
      }
}

void CIvect::print() {
    int block, buf, irrep;

    if (cur_vect_ < 0 || cur_buf_ < 0) {
        outfile->Printf("[Can't print unlocked vector]\n");
    }

    if (vectlen_ > 100000) {
        outfile->Printf("Not printing long (>100000) vector...\n");
        return;
    }

    if (icore_ == 0) {
        for (buf = 0; buf < buf_per_vect_; buf++) {
            read(cur_vect_, buf);
            block = buf2blk_[buf];
            outfile->Printf("\nBlock %2d, codes = (%2d,%2d)\n", block,
                            Ia_code_[block], Ib_code_[block]);
            print_mat(blocks_[block], Ia_size_[block], Ib_size_[block], "outfile");
        }
    }

    else if (icore_ == 1) {
        for (block = 0; block < num_blocks_; block++) {
            outfile->Printf("\nBlock %2d, codes = (%2d,%2d)\n", block,
                            Ia_code_[block], Ib_code_[block]);
            print_mat(blocks_[block], Ia_size_[block], Ib_size_[block], "outfile");
        }
    }

    else if (icore_ == 2) {
        for (buf = 0; buf < buf_per_vect_; buf++) {
            read(cur_vect_, buf);
            irrep = buf2blk_[buf];
            for (block = first_ablk_[irrep]; block <= last_ablk_[irrep];
                 block++) {
                outfile->Printf("\nBlock %2d, codes = (%2d,%2d)\n", block,
                                Ia_code_[block], Ib_code_[block]);
                print_mat(blocks_[block], Ia_size_[block], Ib_size_[block],
                          "outfile");
            }
        }
    }

    else {
        outfile->Printf("(CIvect::print): unrecognized icore option\n");
    }
}

void CIvect::init_vals(int ivect, int nvals, int *alplist, int *alpidx,
      int *betlist, int *betidx, int *blknums, double *value)
{
   int i, j, buf, irrep, blk, ai, bi;

   //ok here it seems safe to set zero blocks
   for (i=0; i<num_blocks_; i++) zero_blocks_[i] = 1;

   /* this used to read >= PARM_GUESS_VEC_H0_BLOCK... but these
      are now gathered from a symnorm so I'll comment this out
      CDS 8/03
   if (CI_Params_->guess_vector == PARM_GUESS_VEC_H0_BLOCK) {
     for (i=0; i<nvals; i++)
        CI_H0block_->c0b[i] = value[i];
     }
   */

   if (icore_ == 1) { /* whole vector in-core */
      zero();
      for (i=0; i<nvals; i++) {
         blk = blknums[i];
         ai = alpidx[i];
         bi = betidx[i];
         blocks_[blk][ai][bi] = value[i];
         zero_blocks_[blk] = 0;
         }
      write(ivect, 0);
      } /* end icore=1 */

   if (icore_ == 2) { /* whole symmetry block in core */
      for (buf=0; buf<buf_per_vect_; buf++) {
         irrep = buf2blk_[buf];
         if (first_ablk_[irrep] < 0) continue;
         zero();
         for (blk=first_ablk_[irrep]; blk<=last_ablk_[irrep]; blk++) {
            for (j=0; j<nvals; j++) {
               if (blknums[j] == blk) {
                  ai = alpidx[j];
                  bi = betidx[j];
                  blocks_[blk][ai][bi] = value[j];
                  zero_blocks_[blk] = 0;
                  }
               }
            } /* end loop over blocks */
         write(ivect, buf);
         } /* end loop over irreps/bufs */

      } /* end icore=2 */

   if (icore_ == 0) { /* one subblock at a time */
      for (buf=0; buf<buf_per_vect_; buf++) {
         zero();
         for (i=0; i<nvals; i++) {
            blk = blknums[i];
            if (blk == buf2blk_[buf]) {
               ai = alpidx[i];
               bi = betidx[i];
               j = ai * Ib_size_[blk] + bi;
               buffer_[j] = value[i];
               zero_blocks_[blk] = 0;
               if (Ms0_)
                  zero_blocks_[decode_[Ib_code_[blk]][Ia_code_[blk]]] = 0;
               }
            }
         write(ivect, buf);
         } /* end loop over buf */
      } /* end icore=0 */

}

void CIvect::set_vals(int ivect, int nvals, int *alplist, int *alpidx,
      int *betlist, int *betidx, int *blknums, double *value)
{
   int i, j, buf, irrep, blk, ai, bi, vec_modified;
   double tval;

   tval = value[0];

   /* this used to read >= PARM_GUESS_VEC_H0_BLOCK... but these
      are now gathered from a symnorm so I'll comment this out
      CDS 8/03
   if (CI_Params_->guess_vector == PARM_GUESS_VEC_H0_BLOCK) {
     for (i=0; i<nvals; i++)
        CI_H0block_->c0b[i] = value[i];
     }
   */

   if (icore_ == 1) { /* whole vector in-core */
      read(ivect, 0);
      vec_modified = 0;
      for (i=0; i<nvals; i++) {
         blk = blknums[i];
         ai = alpidx[i];
         bi = betidx[i];
         blocks_[blk][ai][bi] = value[i];
         zero_blocks_[blk] = 0;
         vec_modified++;
         }
      if (vec_modified) write(ivect, 0);
      } /* end icore=1 */

   if (icore_ == 2) { /* whole symmetry block in core */
      for (buf=0; buf<buf_per_vect_; buf++) {
         vec_modified = 0;
         read(ivect, buf);
         irrep = buf2blk_[buf];
         if (first_ablk_[irrep] < 0) continue;
         for (blk=first_ablk_[irrep]; blk<=last_ablk_[irrep]; blk++) {
            for (j=0; j<nvals; j++) {
               if (blknums[j] == blk) {
                  ai = alpidx[j];
                  bi = betidx[j];
                  blocks_[blk][ai][bi] = value[j];
                  zero_blocks_[blk] = 0;
                  vec_modified++;
                  }
               }
            } /* end loop over blocks */

         if (vec_modified) write(ivect, buf);
         } /* end loop over irreps/bufs */

      } /* end icore=2 */


   if (icore_ == 0) { /* one subblock at a time */
      for (buf=0; buf<buf_per_vect_; buf++) {
         vec_modified = 0;
         read(ivect, buf);
         for (i=0; i<nvals; i++) {
            blk = blknums[i];
            if (blk == buf2blk_[buf]) {
               ai = alpidx[i];
               bi = betidx[i];
               j = ai * Ib_size_[blk] + bi;
               buffer_[j] = value[i];
               zero_blocks_[blk] = 0;
               vec_modified++;
               }
            }
         if (vec_modified) write(ivect, buf);
         } /* end loop over buf */
      } /* end icore=0 */

}


void CIvect::extract_vals(int ivect, int nvals, int *alplist, int *alpidx,
      int *betlist, int *betidx, int *blknums, double *value)
{
   int i, j, buf, irrep, blk, ai, bi, vec_modified;
   double tval;

   tval = value[0];

   if (CI_Params_->guess_vector == PARM_GUESS_VEC_H0_BLOCK) {
     for (i=0; i<nvals; i++)
        CI_H0block_->c0b[i] = value[i];
     }

   if (icore_ == 1) { /* whole vector in-core */
      read(ivect, 0);
      vec_modified = 0;
      for (i=0; i<nvals; i++) {
         blk = blknums[i];
         ai = alpidx[i];
         bi = betidx[i];
         value[i] = blocks_[blk][ai][bi];
         zero_blocks_[blk] = 0;
         vec_modified++;
         }
      if (vec_modified) write(ivect, 0);
      } /* end icore=1 */

   if (icore_ == 2) { /* whole symmetry block in core */
      for (buf=0; buf<buf_per_vect_; buf++) {
         vec_modified = 0;
         read(ivect, buf);
         irrep = buf2blk_[buf];
         if (first_ablk_[irrep] < 0) continue;
         for (blk=first_ablk_[irrep]; blk<=last_ablk_[irrep]; blk++) {
            for (j=0; j<nvals; j++) {
               if (blknums[j] == blk) {
                  ai = alpidx[j];
                  bi = betidx[j];
                  value[j] = blocks_[blk][ai][bi];
                  zero_blocks_[blk] = 0;
                  vec_modified++;
                  }
               }
            } /* end loop over blocks */

         if (vec_modified) write(ivect, buf);
         } /* end loop over irreps/bufs */

      } /* end icore=2 */


   if (icore_ == 0) { /* one subblock at a time */
      for (buf=0; buf<buf_per_vect_; buf++) {
         vec_modified = 0;
         read(ivect, buf);
         for (i=0; i<nvals; i++) {
            blk = blknums[i];
            if (blk == buf2blk_[buf]) {
               ai = alpidx[i];
               bi = betidx[i];
               j = ai * Ib_size_[blk] + bi;
               value[i] = buffer_[j];
               zero_blocks_[blk] = 0;
               vec_modified++;
               }
            }
         if (vec_modified) write(ivect, buf);
         } /* end loop over buf */
      } /* end icore=0 */

}


/*
** CIvect::symnorm()
**
** This function simultaneously symmetrizes and normalizes a CI vector.
** The symmetrization is according to the rule C(Ia,Ib) = (-1)^S C(Ib,Ia)
** and is valid only for Ms=0.
**
** Previously assumed that if Ms!=0, this function is not called.
** Now, re-route such calls to new function CIvect::scale().  Due to
** the old assumption, some required calls to symnorm may be missing
** in the iteration routines.
**
** a = norm
** vecode = 0 for C vector and 1 for Sigma vector
** gather_vec = 0 no gather and 1 for gather
**
*/
void CIvect::symnorm(double a, int vecode, int gather_vec)
{
   int i,j;
   int blk,buf,irrep,ac,bc,len,upper;
   double **mat,*arr, phase, tval;

   if (!Ms0_) {
      scale(a,vecode,gather_vec);
      return;
      }

   if (!CI_Params_->Ms0) phase = 1.0;
   else phase = ((int) CI_Params_->S % 2) ? -1.0 : 1.0;

   if (icore_ == 1) {

      read(cur_vect_, 0);
      for (blk=0; blk<num_blocks_; blk++) {
         ac = Ia_code_[blk];
         bc = Ib_code_[blk];
         mat = blocks_[blk];
         if (ac == bc) { /* diagonal block */
            for (i=0; i<Ia_size_[blk]; i++) {
               mat[i][i] *= a;
               for (j=0; j<i; j++) {
                  mat[i][j] *= a;
                  mat[j][i] = mat[i][j] * phase;
                  }
               }
            }
         if (ac > bc) { /* off-diagonal block */
            xeax(blocks_[blk][0], a, Ia_size_[blk] * Ib_size_[blk]);
            upper = decode_[bc][ac];
            if (upper >= 0) {
               zero_blocks_[upper] = zero_blocks_[blk];
               for (i=0; i<Ia_size_[blk]; i++) {
                  for (j=0; j<Ib_size_[blk]; j++) {
                     blocks_[upper][j][i] = mat[i][j] * phase;
                     }
                  }
               }
            }
         } /* end loop over blocks */

      if (gather_vec) h0block_gather_vec(vecode);
      write(cur_vect_, 0);

      } /* end icore_ == 1 */


   else if (icore_ == 2) { /* irrep at a time */

      for (buf=0; buf<buf_per_vect_; buf++) {
         read(cur_vect_, buf);
         irrep = buf2blk_[buf];
         if (buf_offdiag_[buf]) { /* normalize only. other part never stored */
            for (blk=first_ablk_[irrep]; blk<=last_ablk_[irrep]; blk++) {
                  xeax(blocks_[blk][0], a, Ia_size_[blk] * Ib_size_[blk]);
               }
            }
         else { /* diagonal irrep, symmetrize and normalize */
            for (blk=first_ablk_[irrep]; blk<=last_ablk_[irrep]; blk++) {
               ac = Ia_code_[blk];
               bc = Ib_code_[blk];
               mat = blocks_[blk];
               if (ac == bc) { /* diagonal block */
                  for (i=0; i<Ia_size_[blk]; i++) {
                     mat[i][i] *= a;
                     for (j=0; j<i; j++) {
                        mat[i][j] *= a;
                        mat[j][i] = mat[i][j] * phase;
                        }
                     }
                  }
               if (ac > bc) { /* off-diagonal block in lower triangle */
                  xeax(blocks_[blk][0], a, Ia_size_[blk] * Ib_size_[blk]);
                  upper = decode_[bc][ac];
                  if (upper >= 0) {
                     zero_blocks_[upper] = zero_blocks_[blk];
                     for (i=0; i<Ia_size_[blk]; i++) {
                        for (j=0; j<Ib_size_[blk]; j++) {
                           blocks_[upper][j][i] = mat[i][j] * phase;
                           }
                        }
                     }
                  }
               } /* end loop over blocks */
            } /* end diagonal irrep case */
         if (gather_vec) h0block_gather_vec(vecode);
         write(cur_vect_, buf);
         } /* end loop over buffers */
      } /* end icore==2 */

   else if (icore_ == 0) { /* one RAS block at a time */
      for (buf=0; buf<buf_per_vect_; buf++) {
         blk = buf2blk_[buf];
         read(cur_vect_, buf);
         ac = Ia_code_[blk];
         bc = Ib_code_[blk];
         mat = blocks_[blk];
         if (ac == bc) { /* diagonal block */
            for (i=0; i<Ia_size_[blk]; i++) {
               mat[i][i] *= a;
               for (j=0; j<i; j++) {
                  mat[i][j] *= a;
                  mat[j][i] = mat[i][j] * phase;
                  }
               }
            }
         else { /* off-diagonal block in lower triangle */
            xeax(blocks_[blk][0], a, Ia_size_[blk] * Ib_size_[blk]);
            }

         if (gather_vec) h0block_gather_vec(vecode);
         write(cur_vect_, buf);
         } /* end loop over buffers */

      } /* end case icore==0 */

   else {
     outfile->Printf("(CIvect::symnorm): Unrecognized icore option\n");
      return;
      }

}

/*
   the following subroutine isn't really effective in keeping out some
   lower-energy states that might creep in due to numerical contamination,
   at least not if it is only used after computing the correction vector...
   I added it trying to keep lower-lying sigmas out of delta state
   computations for C2, but it wasn't effective.  ---CDS October 2004
*/

/*
** CIvect::zero_det()
**
** Zero out a specified determinant
** Implement for icore==1 for now... easy to extend.
*/
double CIvect::zero_det(int iac, int ia, int ibc, int ib)
{
  int blk;
  double tval;

  if (icore_ != 1) {
    outfile->Printf( "CIvect::zero_det: Implemented for icore==1 only\n");
    return (0.0);
  }

  blk = decode_[iac][ibc];
  tval = blocks_[blk][ia][ib];
 outfile->Printf("zero_det reports coefficient %12.6lf\n", tval);
  tval = tval*tval;
  blocks_[blk][ia][ib] = 0.0;

  return(tval);
}


/*
** CIvect::scale()
**
** This function scales a CI vector by a given factor.
** Does not concern itself with any possible symmetries or redundancies
** in the CI vector.
**
** Parameters:
**   a = the scale factor
**   vecode = 0 if C vector 1 if Sigma vector
**   gather_vec = 1 if gather 0 otherwise
**
** Returns:
**   none
*/
void CIvect::scale(double a, int vecode, int gather_vec)
{
   int buf;

   for (buf=0; buf<buf_per_vect_; buf++) {
      read(cur_vect_, buf);
      xeax(buffer_, a, buf_size_[buf]);
      if (gather_vec) h0block_gather_vec(vecode);
      write(cur_vect_, buf);
      }
}



/*
** CIvect::buf_lock()
**
** This function "locks in" a memory buffer for the use of a CIvector.
** The appropriate flag is set and pointers are made to point to the
** right regions of the buffer.
**
** Parameters:
**    a  = array to use for the buffer
** Returns: none
*/
void CIvect::buf_lock(double *a)
{
   int i,j,k;

   if (buf_locked_) {
     outfile->Printf("Warning (CIvect::buf_lock): CIvector is already locked!\n");
      }

   if (icore_ == 1) { /* whole vector in-core */
      blocks_[0][0] = a;
      for (j=1; j<Ia_size_[0]; j++) {
         blocks_[0][j] = blocks_[0][0] + Ib_size_[0] * j;
         }
      for (i=1; i<num_blocks_; i++) {
         blocks_[i][0] = blocks_[i-1][0] + Ia_size_[i-1] * Ib_size_[i-1];
         for (j=1; j<Ia_size_[i]; j++) {
            blocks_[i][j] = blocks_[i][0] + Ib_size_[i] * j;
            }
         }
      } /* end icore==1 option */

   if (icore_ == 2) { /* one symmetry block is held at a time */
      blocks_[0][0] = a;
      for (i=0; i<nirreps_; i++) {
         for (j=first_ablk_[i]; j<=last_ablk_[i]; j++) {
            if (j==first_ablk_[i])
               blocks_[j][0] = a;
            else
               blocks_[j][0] = blocks_[j-1][0] + Ia_size_[j-1] * Ib_size_[j-1];
            for (k=1; k<Ia_size_[j]; k++)
               blocks_[j][k] = blocks_[j][0] + Ib_size_[j] * k;
            }
         }
      } /* end icore==2 option */

   if (icore_ == 0) { /* one subblock at a time */
      for (i=0; i<num_blocks_; i++) {
         blocks_[i][0] = a;
         for (j=1; j<Ia_size_[i]; j++) {
            blocks_[i][j] = blocks_[i][0] + Ib_size_[i] * j;
            }
         }
      } /* end icore==0 option */

   buffer_ = a;
   buf_locked_ = 1;
   /* zero(); * commented out 3/13/96: need to eliminate */
}


/*
** CIvect::buf_unlock()
**
** This function "unlocks" a memory buffer from the use of a CIvector.
**
** Parameters: none
** Returns: none
*/
void CIvect::buf_unlock(void)
{
   buf_locked_ = 0;
   blocks_[0][0] = NULL;
   buffer_ = NULL;
   cur_vect_ = -1;
   cur_buf_ = -1;
}



/*
** CIvect::symmetrize(): This function symmetrizes the CI vector
**    to maintain the appropriate spin symmetry.  Symmetrizes only
**    the in-core portion of the vector.
**    Assume that this function is called only if Ms=0
**
** Parameters:
**    phase = the exponent S in C(Ia,Ib) = (-1)^S * C(Ib,Ia)
**
*/
void CIvect::symmetrize(double phase, int iblock)
{
   int i,j;
   int blk,irrep,ac,bc,len,upper;
   double **mat,*arr;

   if (icore_ == 1) {

      for (blk=0; blk<num_blocks_; blk++) {
         ac = Ia_code_[blk];
         bc = Ib_code_[blk];
         mat = blocks_[blk];
         if (ac == bc) { /* diagonal block */
            for (i=0; i<Ia_size_[blk]; i++) {
               for (j=0; j<i; j++) {
                  mat[j][i] = mat[i][j] * phase;
                  }
               }
            }
         if (ac > bc) { /* off-diagonal block */
            upper = decode_[bc][ac];
            if (upper >= 0) {
               zero_blocks_[upper] = zero_blocks_[blk];
               for (i=0; i<Ia_size_[blk]; i++) {
                  for (j=0; j<Ib_size_[blk]; j++) {
                     blocks_[upper][j][i] = mat[i][j] * phase;
                     }
                  }
               }
            }
         } /* end loop over blocks */

      } /* end icore_ == 1 */


   else if (icore_ == 2) { /* irrep at a time */

      irrep = iblock;

      /* do only for diagonal irrep blocks */
      if (CI_CalcInfo_->ref_sym != 0) return;

      for (blk=first_ablk_[irrep]; blk<=last_ablk_[irrep]; blk++) {
         ac = Ia_code_[blk];
         bc = Ib_code_[blk];
         mat = blocks_[blk];
         if (ac == bc) { /* diagonal block */
            for (i=0; i<Ia_size_[blk]; i++) {
               for (j=0; j<i; j++) {
                  mat[j][i] = mat[i][j] * phase;
                  }
               }
            }
         if (ac > bc) { /* off-diagonal block in lower triangle */
            upper = decode_[bc][ac];
            if (upper >= 0) {
               zero_blocks_[upper] = zero_blocks_[blk];
               for (i=0; i<Ia_size_[blk]; i++) {
                  for (j=0; j<Ib_size_[blk]; j++) {
                     blocks_[upper][j][i] = mat[i][j] * phase;
                     }
                  }
               }
            }
         } /* end loop over blocks */
      } /* end icore==2 */

   else if (icore_ == 0) { /* one RAS block at a time */
      ac = Ia_code_[iblock];
      bc = Ib_code_[iblock];
      mat = blocks_[iblock];
      if (ac == bc) { /* diagonal block */
         for (i=0; i<Ia_size_[iblock]; i++) {
            for (j=0; j<i; j++) {
               mat[j][i] = mat[i][j] * phase;
               }
            }
         }

      } /* end case icore==0 */


   else {
     outfile->Printf("(CIvect::symmetrize): Unrecognized icore option\n");
      return;
      }

}



/*
** CIvect::blockptr(): Return a pointer to a block so that it can be
**    accessed outside (yeah, this destroys data protection, but
**    it's better to do this here).
**
*/
double ** CIvect::blockptr(int blknum)
{
   return blocks_[blknum];
}


/*
** CIvect::init_io_files()
**
** Parameters: none
**
** Returns: none
*/
void CIvect::init_io_files(bool open_old) {
    int i;

    for (i = 0; i < nunits_; i++) {
        if (!psio_open_check((ULI)units_[i])) {
            if (open_old) {
                psio_open((ULI)units_[i], PSIO_OPEN_OLD);
            } else {
                psio_open((ULI)units_[i], PSIO_OPEN_NEW);
            }
        }
    }
    fopen_ = true;
}

/*
** CIvect::close_io_files()
**
** Parameters:
**    keep = 1 to keep files, else 0 to delete
**
** Returns: none
*/
void CIvect::close_io_files(int keep) {

    // Nothing to do if its already closed
    if (!fopen_) {
        return;
    }

    for (size_t i = 0; i < nunits_; i++) {
        psio_close(units_[i], keep);
    }
    fopen_ = false;
}

/*
** CIvect::read(): Read in a section of a CI vector from external storage.
**
** Parameters:
**    ivect  = vector number
**    ibuf   = buffer number (ibuf can specify an irrep or a subblock
**                within an irrep, depending on the value of icore.  If
**                icore = 1, then ibuf is ignored.)
**
** Returns: 1 for success, 0 for failure
*/
int CIvect::read(int ivect, int ibuf)
{
   int unit, buf, k, i;
   unsigned long int size;
   int blk;
   char key[20];

   timer_on("CIWave: CIvect read");
   if (nunits_ < 1) {
      cur_vect_ = ivect;
      cur_buf_ = ibuf;
      timer_off("CIWave: CIvect read");
      return(1);
      }

   if (ivect < 0 || ibuf < 0) {
     outfile->Printf("(CIvect::read): Called with negative argument\n");
      timer_off("CIWave: CIvect read");
      return(0);
      }

   if (icore_ == 1) ibuf = 0;
   buf = ivect * buf_per_vect_ + ibuf;

   size = buf_size_[ibuf] * (unsigned long int) sizeof(double);

   /* translate buffer number in case we renumbered after collapse * */
   buf += new_first_buf_;
   if (buf >= buf_total_) buf -= buf_total_;
   sprintf(key, "buffer_ %d", buf);
   unit = file_number_[buf];

   psio_read_entry((ULI) unit, key, (char *) buffer_, size);

   cur_vect_ = ivect;
   cur_buf_ = ibuf;

   timer_off("CIWave: CIvect read");

   return(1);
}


/*
** CIvect::write(): Write a section of a CI vector to external storage.
**
** Parameters:
**    ivect  = vector number
**    ibuf   = buffer number (ibuf can specify an irrep or a subblock
**                within an irrep, depending on the value of icore.  If
**                icore = 1, then ibuf is ignored.)
**
** Returns: 1 for success, 0 for failure
*/
int CIvect::write(int ivect, int ibuf)
{
   int unit, buf, i;
   unsigned long int size;
   int blk;
   char key[20];

   // If we are just an incore buffer
   if (nunits_ < 1) return(1);

   timer_on("CIWave: CIvect write");

   if (ivect >= maxvect_) throw PSIEXCEPTION("(CIvect::write): ivect >= maxvect");
   if (ivect > nvect_) throw PSIEXCEPTION("(CIvect::write): ivect > nvect");
   // {
   //    outfile->Printf( "(CIvect::write): ivect >= maxvect\n");
   //    timer_off("CIWave: CIvect write");
   //    return(0);
   //    }

   if (icore_ == 1) ibuf = 0;
   buf = ivect * buf_per_vect_ + ibuf;
   size = buf_size_[ibuf] * (unsigned long int) sizeof(double);

   /* translate buffer number in case we renumbered after collapse * */
   buf += new_first_buf_;
   if (buf >= buf_total_) buf -= buf_total_;
   sprintf(key, "buffer_ %d", buf);
   unit = file_number_[buf];

   psio_write_entry((ULI) unit, key, (char *) buffer_, size);

   if (ivect >= nvect_) nvect_ = ivect + 1;
   cur_vect_ = ivect;
   cur_buf_ = ibuf;

   timer_off("CIWave: CIvect write");

   return(1);
}


/*
** CIvect::schmidt_add()
**
** This function Gram-Schmidt orthogonalizes a new vector d and adds it to
** the list of vectors in the CIvector c, which must contain room
** for the new vector (i.e. after the new vector is added, nvect <= maxvect).
** Don't add orthogonalized d' if norm(d') < SA_NORM_TOL.
**
** Parameters:
**    L   = number of vectors in CIvect to consider
**
** Returns: 1 if a vector is added, 0 otherwise
**
** Notes: Assumes vectors c,d are same size.  Should account for Ms0 now.
*/
int CIvect::schmidt_add(CIvect &c, int L)
{
   double tval, norm, *dotval;
   int buf, cvect;

   norm = 0.0;

   dotval = init_array(L);

   for (buf=0; buf<buf_per_vect_; buf++) {
      read(cur_vect_, buf);
      for (cvect=0; cvect<L; cvect++) {
         c.read(cvect, buf);
         dot_arr(buffer_, c.buffer_, buf_size_[buf], &tval);
         if (buf_offdiag_[buf]) tval *= 2.0;
         dotval[cvect] += tval;
         }
      }

   for (buf=0; buf<buf_per_vect_; buf++) {
      read(cur_vect_, buf);
      for (cvect=0; cvect<L; cvect++) {
         c.read(cvect, buf);
       /*
         outfile->Printf("dotval[%d] = %2.15lf\n",cvect,dotval[cvect]);
       */
         xpeay(buffer_, -dotval[cvect], c.buffer_, buf_size_[buf]);
         }
      dot_arr(buffer_, buffer_, buf_size_[buf], &tval);
      if (buf_offdiag_[buf]) tval *= 2.0;
      norm += tval;
      write(cur_vect_, buf);
      }

   free(dotval);

   norm = sqrt(norm);
   if (norm < SA_NORM_TOL) return(0);
   norm = 1.0 / norm;

   if (c.nvect_ > c.maxvect_) {
      outfile->Printf( "(CIvect::schmidt_add): no more room to add vectors!\n");
      outfile->Printf( "   c.nvect_ = %d, c.maxvect_ = %d\n", c.nvect_, c.maxvect_);
      return(0);
      }
   else { /* add to c */
      c.cur_vect_ = c.nvect_;
      c.nvect_++;
      for (buf=0; buf<buf_per_vect_; buf++) {
         read(cur_vect_, buf);
         xeay(c.buffer_, norm, buffer_, buf_size_[buf]);
         c.write(c.cur_vect_, buf);
         }
      return(1);
      }

}



/*
** CIvect::schmidt_add2()
**
** This function Gram-Schmidt orthogonalizes a new vector d and adds it to
** the list of vectors in the CIvector c, which must contain room
** for the new vector (i.e. after the new vector is added, nvect <= maxvect).
** Don't add orthogonalized d' if norm(d') < SA_NORM_TOL.  This version
** differs from CIvect::schmidt_add() to allow the user somewhat finer
** control over the numbering of vectors, etc, and allows the return of
** dot products and the normalization factor.
**
** Parameters:
**    L          = number of vectors in CIvect to consider
**    first_vec  = first vec num in C to orthogonalize new vector against
**    last_vec   = last vec num in C to orthogonalize new vector against
**    source_vec = vector number in D file to read
**    target_vec = vector number to write new vector to
**    dotval     = array of dot products of new vector with old vectors
**    nrm        = normalization constant of new vector after orthogonalization
**    ovlpmax    = maximum overlap of current SO vectors to previous vectors
**
** Returns: 1 if a vector is added, 0 otherwise
**
** Notes: Assumes vectors c,d are same size.  Should account for Ms0 now.
*/
int CIvect::schmidt_add2(CIvect &c, int first_vec, int last_vec,
   int source_vec, int target_vec, double *dotval, double *nrm, double *ovlpmax)
{
   double tval, norm, *dotchk, dotck, tmp_norm;
   int buf, cvect, i;

   norm = 0.0;
   tmp_norm = 0.0;
   dotchk = init_array(100);
   *ovlpmax = 0.0;

   for (buf=0; buf<buf_per_vect_; buf++) {
      read(source_vec, buf);
      for (cvect=first_vec; cvect<=last_vec; cvect++) {
         c.read(cvect, buf);
         dot_arr(buffer_, c.buffer_, buf_size_[buf], &tval);
         if (buf_offdiag_[buf]) tval *= 2.0;
         dotval[cvect] += tval;
         }
      }

   for (i=first_vec; i<=last_vec; i++) {
      tval = fabs(dotval[i]);
      if (tval>*ovlpmax) *ovlpmax = tval;
      }

   /* Schmidt orthogonalize and double check orthogonalization */
   for (buf=0; buf<buf_per_vect_; buf++) {
      read(cur_vect_, buf);
      for (cvect=first_vec; cvect<=last_vec; cvect++) {
         c.read(cvect, buf);
         xpeay(buffer_, -dotval[cvect], c.buffer_, buf_size_[buf]);
         }
      dot_arr(buffer_, buffer_, buf_size_[buf], &tval);
      if (buf_offdiag_[buf]) tval *= 2.0;
      norm += tval;
      write(cur_vect_, buf);
      }

   /* outfile->Printf("Norm of %d vec = %20.15f\n",target_vec,norm); */
   norm = sqrt(norm);
   /* outfile->Printf("sqrt Norm of %d vec = %20.15f\n",target_vec,norm); */
   /* There was a dangling else here:
    * if (CI_Params_->mpn_schmidt)
    *   if (norm < MPn_NORM_TOL) return(0);
    *  else if (norm < SA_NORM_TOL) return(0);
    */
   if (CI_Params_->mpn_schmidt) {
     if (norm < MPn_NORM_TOL) return(0);
     else if (norm < SA_NORM_TOL) return(0);
   }
 /*
   if (norm < SA_NORM_TOL && !CI_Params_->mpn) return(0);
 */
   norm = 1.0 / norm;
   /* outfile->Printf("1.0/sqrt(norm) of %d vec = %20.15f\n",target_vec,norm); */
   *nrm = norm;

   if (c.nvect_ > c.maxvect_) {
      outfile->Printf( "(CIvect::schmidt_add2): no more room to add vectors!\n");
      outfile->Printf( "   c.nvect_ = %d, c.maxvect_ = %d\n", c.nvect_, c.maxvect_);
      return(0);
      }
   else { /* add to c */
      c.cur_vect_ = target_vec;
      if (c.cur_vect_ > c.nvect_) c.nvect_++;
      zero_arr(dotchk,100);

      for (buf=0; buf<buf_per_vect_; buf++) {
         read(cur_vect_, buf);
         xeay(c.buffer_, norm, buffer_, buf_size_[buf]);
         c.write(c.cur_vect_, buf);
         }
     /*
      outfile->Printf( "c.cur_vect_ = %d\n",c.cur_vect_);
      outfile->Printf( "dot product of normalized vector in SA2 = %20.10f\n",
              tmp_norm);
      c.print(outfile);
     */

      if (CI_Params_->mpn) {
        zero_arr(dotchk,100);
        for (buf=0; buf<buf_per_vect_; buf++) {
           read(source_vec, buf);
           for (cvect=first_vec; cvect<=last_vec; cvect++) {
              c.read(cvect, buf);
              dot_arr(buffer_, c.buffer_, buf_size_[buf], &tval);
              if (buf_offdiag_[buf]) tval *= 2.0;
              dotchk[cvect] += tval;
              }
           }
        for (i=first_vec; i<=last_vec; i++)
           if (dotchk[i] > *ovlpmax) *ovlpmax = dotchk[i];
        }
      return(1);
      }

}



/*
** CIvect::zero()
**
** Zero out the current memory buffer for a CI vector.
**
** Parameters: none
** Returns: none
**/
void CIvect::zero(void)
{
   zero_arr(buffer_, (int) buffer_size_);
}



/*
** CIvect::sigma_renorm()
**
** Function calculates the numerator part of the Davidson correction vector d
**
** Parameters:
**    nr       = number of roots (=number of d vectors to calculate)
**    L        = number of previous vectors in CI subspace
**    alpha    = subspace CI eigenvector matrix
**    lambda   = array of subspace eigenvalues
**    norm_arr = norm array (hold norm of for each d vector)
**    C        = CIvect for subspace vectors
**    printflag= 1 to print d vector(s), else 0
**    outfile  = where to put any output
**
** Returns: none
*/
void CIvect::sigma_renorm(int nr, int L, double renorm_C, CIvect &S,
             double *buf1, int printflag)
{
   int buf, ivect, root;
   double tval;

      for (buf=0; buf<buf_per_vect_; buf++) {
         for (ivect=0; ivect<L; ivect++) {
            S.buf_lock(buf1);
            S.read(ivect, buf);
            xeay(S.buffer_, renorm_C, S.buffer_, buf_size_[buf]);
            S.buf_unlock();
            } /* end loop over ivect */

         write(nr, buf);
         if (printflag) {
            outfile->Printf( "\nSigma renormalized matrix\n");
            print_buf();
            }
         } /* loop over buffers */
}


/*
** CIvect::dcalc()
**
** Function calculates the numerator part of the Davidson correction vector d
**
** Parameters:
**    nr       = number of roots (=number of d vectors to calculate)
**    L        = number of previous vectors in CI subspace
**    alpha    = subspace CI eigenvector matrix
**    lambda   = array of subspace eigenvalues
**    norm_arr = norm array (hold norm of for each d vector)
**    C        = CIvect for subspace vectors
**    S        = CIvect for sigma vectors
**    root_converged = 1 if root has converged, else 0
**    printflag= 1 to print d vector(s), else 0
**    outfile  = where to put any output
**    E_est     = Intermediate is OLSEN update
**
** Returns: none
*/
void CIvect::dcalc(int nr, int L, double **alpha, double *lambda,
      double *norm_arr, CIvect &C, CIvect &S, double *buf1, double *buf2,
      int *root_converged, int printflag, double *E_est)
{
   int buf, ivect, root, tmproot, converged=0, i;
   double tval;


   buf_lock(buf2);

   /* Calculate the d vector for each converged root but do
   ** not form the complete correction (f) vector
   */

   /*
   for (root=0; root<nr; root++)
      if (root_converged[root]) {
        converged = root;
        break;
        }

   if (converged) nr++;
   */

   for (root=0; root<nr; root++) {
      norm_arr[root] = 0.0;
      /* if (converged && root==nr-1 && L<nr) break; */
      for (buf=0; buf<buf_per_vect_; buf++) {
         zero();
         if (CI_Params_->update==UPDATE_OLSEN) {
           read(root,buf);
           xeax(buffer_, -E_est[root], buf_size_[buf]);
           /* buffer is know E_est*C^k */
           }
         for (ivect=0; ivect<L; ivect++) {
            if (CI_Params_->update == UPDATE_DAVIDSON) { /* DAVIDSON update formula */
              C.buf_lock(buf1);
              C.read(ivect, buf);
              tval = -alpha[ivect][root] * lambda[root];
              xpeay(buffer_, tval, C.buffer_, buf_size_[buf]);
              C.buf_unlock();
              }
            S.buf_lock(buf1);
            S.read(ivect, buf);
            xpeay(buffer_, alpha[ivect][root], S.buffer_, buf_size_[buf]);
            S.buf_unlock();
            } /* end loop over ivect */
         dot_arr(buffer_, buffer_, buf_size_[buf], &tval);
         if (buf_offdiag_[buf]) tval *= 2.0;
         norm_arr[root] += tval;
         write(root, buf);
         /*
         if (root==nr-1 && converged) write(converged, buf);
         else write(root, buf);
         */

         if (printflag) {
            outfile->Printf( "\nfirst D matrix\n");
            print_buf();
            }
         } /* loop over buffers */

      norm_arr[root] = sqrt(norm_arr[root]);

      } /* end loop over roots */

   buf_unlock();
}



/*
**
** CIvect::dcalc2()
**
** Function calculates the denominator part of the Davidson correction
** vector 'd'.
**
** Parameters:
**    rootnum   = number of current root
**    lambda    = current iteration's energy eigenvalue for the current root
**    Hd        = CIvector for the Hamiltonian diagonal
**
** Returns: sum of squares of coefficients in d.
*/
double CIvect::dcalc2(int rootnum, double lambda, CIvect &Hd,
        int precon, struct stringwr **alplist, struct stringwr **betlist)
{
   int buf, errcod, i;
   double tval, norm = 0.0;

   for (buf=0; buf<buf_per_vect_; buf++) {
      read(rootnum, buf);
      if (CI_Params_->hd_otf == FALSE) Hd.read(0, buf);
      else if (CI_Params_->hd_otf == TRUE) {
        if (CI_Params_->mpn)
          Hd.diag_mat_els_otf(alplist, betlist, CI_CalcInfo_->onel_ints->pointer(),
             CI_CalcInfo_->twoel_ints->pointer(), CI_CalcInfo_->e0_drc, CI_CalcInfo_->num_alp_expl,
             CI_CalcInfo_->num_bet_expl, CI_CalcInfo_->nmo, buf, CI_Params_->hd_ave);
        else
          Hd.diag_mat_els_otf(alplist, betlist, CI_CalcInfo_->onel_ints->pointer(),
             CI_CalcInfo_->twoel_ints->pointer(), CI_CalcInfo_->edrc, CI_CalcInfo_->num_alp_expl,
             CI_CalcInfo_->num_bet_expl, CI_CalcInfo_->nmo, buf, CI_Params_->hd_ave);
        }

      if (CI_Params_->mpn) norm = calc_mpn_vec(buffer_, lambda, Hd.buffer_,
                           buf_size_[buf], 1.0, -1.0, DIV);
      else {
        if (CI_Params_->precon >= PRECON_GEN_DAVIDSON)
          h0block_gather_vec(CI_VEC);
        tval = calc_d2(buffer_, lambda, Hd.buffer_, buf_size_[buf], precon);
       }

      if (buf_offdiag_[buf]) tval *= 2.0;
      norm += tval;
      write(rootnum, buf);
      }
   //if (!CI_Params_->mpn) errcod = H0block_calc(lambda); /* MLL */
   return(norm);
}

/*
**
** CIvect::construct_kth_order_wf()
**
** Function constructs the kth order wavefunction from all
** other nth order wavefunction (n<k) and the current Sigma vector
** Hc_k-1.  Uses only the two buffers.
**
** Parameters:
**    Hd         = CIvect for H0
**    S          = CIvect for Sigma vector
**    C          = CIvect for Cvec vector
**    alplist    = alpha string list
**    betlist    = beta  string list
**    buf1       = first buffer for a CIvect
**    buf2       = second buffer for a CIvect
**
** Returns: none
*/
void CIvect::construct_kth_order_wf(CIvect &Hd, CIvect &S, CIvect &C,
       struct stringwr **alplist, struct stringwr **betlist, double *buf1,
       double *buf2, int k, double *mp_energy, double **cvec_coeff,
       double *cvec_norm)
{
   int i, j, r, order, buf, block;
   double tval, norm;

   //outfile->Printf("\nCVEC_COEFF and CVEC_NORMS in CONSTRUCT\n");
   //print_mat(cvec_coeff, k-2, k-2, outfile);
   //for (i=0; i<k-2; i++)
   //   outfile->Printf("cvec_norm[%d] = %lf\n",i, cvec_norm[i]);
   //

     for (buf=0; buf<buf_per_vect_; buf++) {
        Hd.buf_lock(buf2);
        Hd.diag_mat_els_otf(alplist, betlist, CI_CalcInfo_->onel_ints->pointer(),
             CI_CalcInfo_->twoel_ints->pointer(), CI_CalcInfo_->e0_drc, CI_CalcInfo_->num_alp_expl,
             CI_CalcInfo_->num_bet_expl, CI_CalcInfo_->nmo, buf, CI_Params_->hd_ave);
        read(k-1, buf);
        norm = calc_mpn_vec(buffer_, (mp_energy[1]-CI_CalcInfo_->edrc),
                Hd.buffer_, buf_size_[buf], 1.0, 1.0, MULT);
        Hd.buf_unlock();

        C.buf_lock(buf2);
        if (CI_Params_->mpn_schmidt) {
          for (i=0; i<=k-2; i++) {
             C.read(i, buf);
             tval = 0.0;
             for (r=2; r<=k; r++) {
               if ((k-r)==i) tval+= mp_energy[r]*cvec_coeff[k-r][i]
                                    *(1.0/cvec_norm[k-r]);
               else tval+= mp_energy[r]*cvec_coeff[k-r][i];
               }
             xpeay(buffer_, tval, C.buffer_, buf_size_[buf]);
             }
          }
        else {
          for (i=2; i<=k; i++) {
             C.read(k-i, buf);
             xpeay(buffer_, mp_energy[i], C.buffer_, buf_size_[buf]);
             }
          }
        C.buf_unlock();

        S.buf_lock(buf2);
        S.read(0, buf);
        xeaxmy(buffer_, S.buffer_, 1.0, S.buf_size_[buf]);
        S.buf_unlock();

        Hd.buf_lock(buf2);
        Hd.diag_mat_els_otf(alplist, betlist, CI_CalcInfo_->onel_ints->pointer(),
             CI_CalcInfo_->twoel_ints->pointer(), CI_CalcInfo_->e0_drc, CI_CalcInfo_->num_alp_expl,
             CI_CalcInfo_->num_bet_expl, CI_CalcInfo_->nmo, buf, CI_Params_->hd_ave);
        norm = calc_mpn_vec(buffer_, CI_CalcInfo_->e0, Hd.buffer_, buf_size_[buf],
                -1.0, 1.0, DIV);

        if (Ms0_) {
          block = buf2blk_[buf];
          if ((int) CI_Params_->S % 2) symmetrize(-1.0, block);
          else symmetrize(1.0, block);
         }
        copy_zero_blocks(S);
        write(k, buf);
        Hd.buf_unlock();
        }

}

/*
**
** CIvect::wigner_E2k_formula()
**
** Uses the kth order wavefunction and the Wigner formulas
**   to compute the 2k and 2k+1 th order energies
**
** Parameters:
**    Hd         = CIvect for H0
**    S          = CIvect for Sigma vector
**    C          = CIvect for Cvec vector
**    alplist    = alpha string list
**    betlist    = beta  string list
**    buf1       = first buffer for a CIvect
**    buf2       = second buffer for a CIvect
**    k          = kth order
**    mp2k_energy= array for storing the energy values
**
** Returns: none
*/
void CIvect::wigner_E2k_formula(CIvect &Hd, CIvect &S, CIvect &C,
       struct stringwr **alplist, struct stringwr **betlist, double *buf1,
       double *buf2, int k, double *mp2k_energy, double **wfn_overlap,
       double **cvec_coeff, double *cvec_norm, int kvec_offset)
{
   int i, j, buf, I, J;
   double tval, E2k, E2kp1, tval2;

   E2k = E2kp1 = 0.0;

   /* First determine the overlap of kth order wavefunction with
   ** all previous order wavefunctions
   */
   if (CI_Params_->mpn_schmidt) {
     zero_mat(wfn_overlap, CI_Params_->maxnvect+1, CI_Params_->maxnvect+1);
     /* for (i=0; i<k-1; i++) wfn_overlap[i][i] = 1.0; */
     }

   C.buf_lock(buf2);
   for (buf=0; buf<buf_per_vect_; buf++) {
      if (CI_Params_->mpn_schmidt) {
        for (i=0; i<=(k-kvec_offset); i++) {
           read(i, buf);
           for (j=i; j<=(k-kvec_offset); j++) {
              C.read(j,buf);
              dot_arr(buffer_, C.buffer_, C.buf_size_[buf], &tval);
              if (buf_offdiag_[buf]) tval *= 2.0;
              wfn_overlap[i+kvec_offset][j+kvec_offset] += tval;
              if (i!=j) wfn_overlap[j+kvec_offset][i+kvec_offset] += tval;
              }
           }
      /*
        read(k-1, buf);
        for (i=1; i<=k-1; i++) {
           C.read(i, buf);
           dot_arr(buffer_, C.buffer_, C.buf_size_[buf], &tval);
           if (buf_offdiag_[buf]) tval *= 2.0;
           wfn_overlap[k-1][i] += tval;
           if (i!=(k-1)) wfn_overlap[i][k-1] += tval;
           }
       */
        }
      else {
        read(k, buf);
        for (i=(1-kvec_offset); i<=(k-kvec_offset); i++) {
           C.read(i, buf);
           dot_arr(buffer_, C.buffer_, C.buf_size_[buf], &tval);
           if (buf_offdiag_[buf]) tval *= 2.0;
           wfn_overlap[k][i+kvec_offset] += tval;
           if ((i+kvec_offset)!=k) wfn_overlap[i+kvec_offset][k] += tval;
           }
        }
      }
   C.buf_unlock();


   if (print_lvl_ > 3) {
     outfile->Printf("\nwfn_overlap = \n");
     print_mat(wfn_overlap, k+1, k+1, "outfile");
     outfile->Printf("\t\t\t\t");
     }

   /* Compute E_2k and E_2k+1 */

   for (buf=0; buf<buf_per_vect_; buf++) {
      S.buf_lock(buf2);
      S.read(0, buf);
      read(k-1-kvec_offset, buf);
      dot_arr(buffer_, S.buffer_, buf_size_[buf], &tval);
      if (buf_offdiag_[buf]) tval *= 2.0;
      E2k += tval;
      read(k-kvec_offset, buf);
      dot_arr(buffer_, S.buffer_, buf_size_[buf], &tval);
      if (buf_offdiag_[buf]) tval *= 2.0;
      E2kp1 += tval;
      S.buf_unlock();
      Hd.buf_lock(buf2);
      Hd.diag_mat_els_otf(alplist, betlist, CI_CalcInfo_->onel_ints->pointer(),
           CI_CalcInfo_->twoel_ints->pointer(), CI_CalcInfo_->e0_drc, CI_CalcInfo_->num_alp_expl,
           CI_CalcInfo_->num_bet_expl, CI_CalcInfo_->nmo, buf, CI_Params_->hd_ave);
      xexy(Hd.buffer_, buffer_, buf_size_[buf]);
      dot_arr(buffer_, Hd.buffer_, buf_size_[buf], &tval);
      if (buf_offdiag_[buf]) tval *= 2.0;
      E2kp1 -= tval;
      read(k-1-kvec_offset, buf);
      dot_arr(buffer_, Hd.buffer_, buf_size_[buf], &tval);
      if (buf_offdiag_[buf]) tval *= 2.0;
      E2k -= tval;
      Hd.buf_unlock();
      }

   if (CI_Params_->mpn_schmidt) {
   /*
     C.buf_lock(buf2);
     for (i=1; i<=k-2; i++) {
        zero_arr(buffer_, buf_size_[0]);
        for (I=1; I<=k-2; I++) {
           C.read(I,0);
           outfile->Printf( " prescaled Bvec %d = \n", I);
           C.print(outfile);
           if (i==I) tval = cvec_coeff[i][I]*(1.0/cvec_norm[i]);
           else tval = cvec_coeff[i][I];
           xpeay(buffer_, tval, C.buffer_, buf_size_[0]);
           }
        if (i==1) {
          outfile->Printf( " Cvec %d = \n", i);
          print(outfile);
          }
        }
     C.buf_unlock();
    */
     for (i=1-kvec_offset; i<=k-2-kvec_offset; i++) {
        for (j=1-kvec_offset; j<=k-2-kvec_offset; j++) {
           tval = 0.0;
           for (I=1-kvec_offset; I<=k-2-kvec_offset; I++) {
              if (I==j && I==i) tval += cvec_coeff[i+kvec_offset][I+kvec_offset]*
                                        cvec_coeff[j+kvec_offset][I+kvec_offset]*
                                        (1.0/cvec_norm[i+kvec_offset])*
                                        (1.0/cvec_norm[j+kvec_offset]);
              else if (I==i) tval += cvec_coeff[i+kvec_offset][I+kvec_offset]*
                                     cvec_coeff[j+kvec_offset][I+kvec_offset]*
                                     (1.0/cvec_norm[i+kvec_offset]);
              else if (I==j) tval += cvec_coeff[i+kvec_offset][I+kvec_offset]*
                                     cvec_coeff[j+kvec_offset][I+kvec_offset]*
                                     (1.0/cvec_norm[j+kvec_offset]);
              else tval += cvec_coeff[i+kvec_offset][I+kvec_offset]*
                           cvec_coeff[j+kvec_offset][I+kvec_offset];
              }
           E2k -= tval*mp2k_energy[2*k+kvec_offset-i+kvec_offset-j+kvec_offset];
           E2kp1 -= tval*mp2k_energy[2*k+1+kvec_offset-i+kvec_offset-j+kvec_offset];
           }
        tval = tval2 = 0.0;
        for (I=1-kvec_offset; I<=k-kvec_offset-2; I++) {
           if (I==i) {
             tval += wfn_overlap[k][I]*cvec_coeff[i][I]*(1.0/cvec_norm[i]);
             tval2 += wfn_overlap[k-1][I]*cvec_coeff[i][I]*(1.0/cvec_norm[i]);
             }
           else {
             tval += wfn_overlap[k][I]*cvec_coeff[i][I];
             tval2 += wfn_overlap[k-1][I]*cvec_coeff[i][I];
            }
           }
        E2kp1 -= tval*2.0*mp2k_energy[k+1-i];
        E2kp1 -= tval2*2.0*mp2k_energy[k+2-i];
        E2k -= tval*mp2k_energy[k-i];
        E2k -= tval2*2.0*mp2k_energy[k+1-i];
       /*
        outfile->Printf( "E2kp1 -> - tval*2.0*mp2k_energy[k+1-i] = %20.10f\n",
                tval*2.0*mp2k_energy[k+1-i]);
        outfile->Printf( "E2kp1 -> - tval2*2.0*mp2k_energy[k+2-i] = %20.10f\n",
                tval*2.0*mp2k_energy[k+2-i]);
        outfile->Printf( "E2k -> - tval*mp2k_energy[k-i] = %20.10f\n",
                tval*mp2k_energy[k-i]);
        outfile->Printf( "E2k -> - tval*2.0*mp2k_energy[k+1-i] = %20.10f\n",
                tval*2.0*mp2k_energy[k+1-i]);
        */
        }
     E2kp1 += (CI_CalcInfo_->edrc-mp2k_energy[1])*wfn_overlap[k][k];
     E2kp1 -= 2.0*mp2k_energy[2]*wfn_overlap[k-1][k];
     E2kp1 -= mp2k_energy[3]*wfn_overlap[k-1][k-1];
     E2k += (CI_CalcInfo_->edrc-mp2k_energy[1])*wfn_overlap[k][k-1];
     E2k -= mp2k_energy[2]*wfn_overlap[k-1][k-1];
     }

  /*
   else {
     C.buf_lock(buf2);
     for (I=1; I<=k; I++) {
        C.read(I,0);
        if (I==1) {
          outfile->Printf("Cvec %d = \n", I);
          C.print(outfile);
          }
        }
     C.buf_unlock();
     for (i=1; i<=k-2; i++) {
        E2kp1 -= 2.0 * mp2k_energy[k+1-i] * wfn_overlap[i][k];
        E2kp1 -= 2.0 * mp2k_energy[k+2-i] * wfn_overlap[i][k-1];
        for (j=1; j<=k-2; j++)
           E2kp1 -= mp2k_energy[2*k+1-i-j] * wfn_overlap[i][j];
        }

     for (i=1; i<=k-2; i++) {
        E2k -= mp2k_energy[k-i] * wfn_overlap[k][i];
        E2k -= 2.0 * mp2k_energy[k+1-i] * wfn_overlap[k-1][i];
        for (j=1; j<=k-2; j++)
           E2k -= mp2k_energy[2*k-i-j] * wfn_overlap[i][j];
        }
     E2k += (CI_CalcInfo_->edrc-mp2k_energy[1]) * wfn_overlap[k][k-1];
     E2k -= mp2k_energy[2] * wfn_overlap[k-1][k-1];
     E2kp1 += (CI_CalcInfo_->edrc-mp2k_energy[1])*wfn_overlap[k][k];
     E2kp1 -= mp2k_energy[3]*wfn_overlap[k-1][k-1];
     E2kp1 -= 2.0 * mp2k_energy[2] * wfn_overlap[k-1][k];
     }
  */

   else {
     for (i=1; i<=k; i++)
        for (j=1; j<=k; j++) {
           E2kp1 -= mp2k_energy[2*k+1-i-j] * wfn_overlap[i][j];
           if ((i==k) && (j==k)) E2kp1 += CI_CalcInfo_->edrc * wfn_overlap[k][k];
           }

     for (i=1; i<=k; i++)
        for (j=1; j<k; j++) {
           E2k -= mp2k_energy[2*k-i-j] * wfn_overlap[i][j];
           if ((i==k) && (j==k-1)) E2k += CI_CalcInfo_->edrc * wfn_overlap[k][k-1];
           }
     }

  /*  outfile->Printf( "final E2k = %lf\n", E2k);
   outfile->Printf( "final E2kp1 = %lf\n", E2kp1);
  */

   mp2k_energy[2*k] = E2k;
   mp2k_energy[2*k+1] = E2kp1;

}


/*
** CIvect::print_buf()
**
** Function prints the in-core section of a CI vector
**
** Parameters:
**    outfile  = file pointer for output file
**
** Returns: none
*/
void CIvect::print_buf()
{
   int blk;
   int irrep;

   if (icore_ == 1) {
      for (blk=0; blk<num_blocks_; blk++) {
         outfile->Printf( "\nBlock %2d, codes = (%2d,%2d)\n", blk,
            Ia_code_[blk], Ib_code_[blk]);
         print_mat(blocks_[blk], Ia_size_[blk], Ib_size_[blk], "outfile");
         }
      }

   if (icore_ == 2) { /* symmetry block in-core */
      irrep = buf2blk_[cur_buf_];
      if (first_ablk_[irrep] < 0) {
         outfile->Printf( "(CIvect::print_blk): No blks for irrep %d\n",irrep);
         return;
         }
      else {
         for (blk=first_ablk_[irrep]; blk <= last_ablk_[irrep]; blk++) {
            outfile->Printf( "\nBlock %2d, codes = (%2d,%2d)\n", blk,
               Ia_code_[blk], Ib_code_[blk]);
            print_mat(blocks_[blk], Ia_size_[blk], Ib_size_[blk], "outfile");
            }
         }
      }

   if (icore_ == 0) { /* one subblock in-core */
      blk = buf2blk_[cur_buf_];

      outfile->Printf( "\nBlock %2d, codes = (%2d,%2d)\n", blk,
         Ia_code_[blk], Ib_code_[blk]);
      print_mat(blocks_[blk], Ia_size_[blk], Ib_size_[blk], "outfile");
      }
}



/*
** CIvect::civ_xeay()
**
** Function does the operation X = a * Y for two CI vectors X and Y and
** some constant a.
**
** Parameters:
**    a     =  constant in X = a * Y
**    Y     =  vector in X = a * Y
**    xvect = vector number in X
**    yvect = vector number in Y
**
** Returns: none
*/
void CIvect::civ_xeay(double a, CIvect &Y, int xvect, int yvect)
{
   int buf;

   for (buf=0; buf<buf_per_vect_; buf++) {
      Y.read(yvect, buf);
      xeay(buffer_, a, Y.buffer_, buf_size_[buf]);
      write(xvect, buf);
      }
}



/*
** CIvect::civ_xpeay()
**
** Function does the operation X += a * Y for two CI vectors X and Y and
** some constant a.
**
** Parameters:
**    a     =  constant in X = a * Y
**    Y     =  vector in X = a * Y
**    xvect = vector number in X
**    yvect = vector number in Y
**
** Returns: none
*/
void CIvect::civ_xpeay(double a, CIvect &Y, int xvect, int yvect)
{
   int buf;

   for (buf=0; buf<buf_per_vect_; buf++) {
      Y.read(yvect, buf);
      read(xvect, buf);
      xpeay(buffer_, a, Y.buffer_, buf_size_[buf]);
      write(xvect, buf);
      }
}


/*
** CIvect::transp_block()
**
** Transpose a CI vector (or a piece of a CI vector)
**
** Parameters:
**    iblock = RAS subblock number
**    tmparr = scratch array to use as intermediate
**
** Returns: none
*/
void CIvect::transp_block(int iblock, double **tmparr)
{
   int i,j,nrows,ncols;
   double **src, *dest;

   src = blocks_[iblock];
   dest = tmparr[0];

   /* bind pointers to subbblock topology */
   nrows = Ib_size_[iblock];
   ncols = Ia_size_[iblock];
   for (i=1; i<nrows; i++) {
      tmparr[i] = dest + i * ncols;
      }

   /* copy data */
   for (i=0; i<Ib_size_[iblock]; i++) {
      for (j=0; j<Ia_size_[iblock]; j++) {
         *dest++ = src[j][i];
         }
      }

}


/*
** CIvect::get_max_blk_size()
**
** Return the maximum RAS subblock size as a long unsigned integer
**
*/
unsigned long CIvect::get_max_blk_size(void)
{
   int i;
   unsigned long blksize, maxblksize=0;

   for (i=0; i<num_blocks_; i++) {
      blksize = (unsigned long) Ia_size_[i] * (unsigned long) Ib_size_[i];
      if (blksize > maxblksize) maxblksize = blksize;
      }

   return(maxblksize);
}


/*
** CIvect::checknorm()
**
** Check the norm of a CI vector
*/
double CIvect::checknorm(void) {
    double tval, dotprod = 0.0;

    for (int buf = 0; buf < buf_per_vect_; buf++) {
        read(cur_vect_, buf);
        dot_arr(buffer_, buffer_, buf_size_[buf], &tval);
        if (buf_offdiag_[buf]) tval *= 2.0;
        dotprod += tval;
    }

    return (dotprod);
}

/*
** CIvect::copy()
**
** This copies one CI vector to another
**
*/
void CIvect::copy(CIvect &Src, int targetvec, int srcvec) {

    for (int buf = 0; buf < buf_per_vect_; buf++) {
        Src.read(srcvec, buf);
        xey(buffer_, Src.buffer_, buf_size_[buf]);
        int blk = buf2blk_[buf];
        if ((blk >= 0) && ((zero_blocks_[blk] == 0) || (Src.zero_blocks_[blk] == 0)))
            zero_blocks_[blk] = 0;
        write(targetvec, buf);
    }
}

/*
** CIvect::restart_gather()
**
** This function takes a linear combination of previous 'b' vectors to
** make a current 'c' vector during a Davidson iteration procedure.
** The coefficients are given by array alpha, the current vector number
** is given by ivec, and the number of previous vectors is nvec.
**
*/
void CIvect::restart_gather(int ivec, int nvec, int nroot, double **alpha,
      double *buffer_1, double *buffer_2)
{
   int buf, oldvec;

   for (buf=0; buf<buf_per_vect_; buf++) {
      zero_arr(buffer_2, buf_size_[buf]);
      buf_lock(buffer_1);
      for (oldvec=0; oldvec<nvec; oldvec++) {
         read(oldvec, buf);
         xpeay(buffer_2, alpha[oldvec][nroot], buffer_1, buf_size_[buf]);
         }
      buf_unlock();
      buf_lock(buffer_2);
      write(ivec, buf);
      buf_unlock();
      }
}



/*
** CIvect::gather()
**
** This function takes a linear combination of previous 'b' vectors to
** make a current 'c' vector during a Davidson iteration procedure.
** The coefficients are given by array alpha, the current vector number
** is given by ivec, and the number of previous vectors is nvec.
** Similar to CIvect::restart_gather() except that we no longer assume
** two different CIvects.  Assume buffer locks are already done.
**
*/
void CIvect::gather(int ivec, int nvec, int nroot, double **alpha,
      CIvect &C)
{
   int buf, oldvec;

   /* outfile->Printf("In CIvect::gather\n"); */
   for (buf=0; buf<buf_per_vect_; buf++) {
      zero_arr(buffer_, buf_size_[buf]);
      for (oldvec=0; oldvec<nvec; oldvec++) {
         C.read(oldvec, buf);
         xpeay(buffer_, alpha[oldvec][nroot], C.buffer_, buf_size_[buf]);
         /* outfile->Printf("coef[%d][%d] = %10.7f\n",oldvec,nroot,alpha[oldvec][nroot]); */
         }
      write(ivec, buf);
      }
}


/*
** CIvect::restart_reord_fp()
**
** This function reorders the file pointers to do a restart of the iterative
** diagonalization method.  The idea is that, for nroot > 1, it is not
** possible to restart and make CI vectors 0...nroot equal to the restarted
** approximate eigenvectors, since this overwrites info needed to construct
** them.  Thus, the new CI vectors must occupy the LAST nroot positions.
** However, it is nevertheless useful to index them as 0...nroot.  This
** is most easily accomplished by a remapping (rotation) of the file
** pointer info.  The one parameter is L, the new "0" vector.
**
** Actually, it's slightly more complex.  For multiple restarts in a given
** calc, the 0 position rotates around.  This routine should still work.
**
** In the latest version, I am phasing out the "offset" array in favor
** of a new_first_buf array which basically gives the buffer number of
** the new "0" vector.  This is more natural for the libpsio implementation.
*/
void CIvect::restart_reord_fp(int L)
{
   int buf, newbuf;
   int *tmp_file_number_;

   new_first_buf_ = L*buf_per_vect_ + new_first_buf_;
   if (new_first_buf_ >= buf_total_) new_first_buf_ -= buf_total_;

   /*
   tmp_file_offset_ = (unsigned long *) malloc (buf_total_ *
                      sizeof(unsigned long));
   tmp_file_number_ = init_int_array(buf_total_);


   for (buf=L*buf_per_vect_,newbuf=0; buf<buf_total_; buf++,newbuf++) {
     tmp_file_offset_[newbuf] = file_offset_[buf];
     tmp_file_number_[newbuf] = file_number_[buf];
   }

   for (buf=0; buf<L*buf_per_vect_; buf++,newbuf++) {
     tmp_file_offset_[newbuf] = file_offset_[buf];
     tmp_file_number_[newbuf] = file_number_[buf];
   }

   for (buf=0; buf<buf_total_; buf++) {
     file_offset_[buf] = tmp_file_offset_[buf];
     file_number_[buf] = tmp_file_number_[buf];
   }

   free(tmp_file_offset_);
   free(tmp_file_number_);
   */
}


void CIvect::print_fptrs()
{
   int buf;

   outfile->Printf("Printing file pointer information\n");
   for (buf=0; buf<buf_total_; buf++)
      outfile->Printf("%d %d\n",buf, file_number_[buf]);


}



/*
** CIvect::h0block_buf_init()
**
** Initialize H0block stuff pertaining to buffers
**
*/
void CIvect::h0block_buf_init(void)
{
   int i, cnt, irrep, buf, blk;
   int *tmparr;

   CI_H0block_->nbuf = buf_per_vect_;
   CI_H0block_->buf_num = init_int_array(buf_per_vect_);
   if (CI_H0block_->size < 1) return;

   tmparr = init_int_array(CI_H0block_->size+CI_H0block_->coupling_size);

   if (icore_ == 1) {
      CI_H0block_->buf_member =
        init_int_matrix(1, CI_H0block_->size+CI_H0block_->coupling_size);
      for (i=0; i<(CI_H0block_->size+CI_H0block_->coupling_size); i++) {
         CI_H0block_->buf_member[0][i] = i;
         }
      CI_H0block_->buf_num[0] = CI_H0block_->size+CI_H0block_->coupling_size;
      }
   else if (icore_ == 2) {
      CI_H0block_->buf_member = (int **) malloc (buf_per_vect_ * sizeof(int *));
      for (buf=0; buf<buf_per_vect_; buf++) {
         cnt = 0;
         irrep = buf2blk_[buf];
         for (blk=first_ablk_[irrep]; blk<=last_ablk_[irrep]; blk++) {
            for (i=0; i<CI_H0block_->size+CI_H0block_->coupling_size; i++) {
               if (CI_H0block_->blknum[i] == blk) {
                  tmparr[cnt++] = i;
                  }
               }
            }
         CI_H0block_->buf_num[buf] = cnt;
         if (cnt) CI_H0block_->buf_member[buf] = init_int_array(cnt);
         for (i=0; i<cnt; i++) {
            CI_H0block_->buf_member[buf][i] = tmparr[i];
            }
         } /* end loop over bufs */
      } /* end icore==2 */
   else {
      CI_H0block_->buf_member = (int **) malloc (buf_per_vect_ * sizeof(int *));
      for (buf=0; buf<buf_per_vect_; buf++) {
         cnt = 0;
         blk = buf2blk_[buf];
         for (i=0; i<CI_H0block_->size+CI_H0block_->coupling_size; i++) {
            if (CI_H0block_->blknum[i] == blk) {
               tmparr[cnt++] = i;
               }
            }
         CI_H0block_->buf_num[buf] = cnt;
         if (cnt) CI_H0block_->buf_member[buf] = init_int_array(cnt);
         for (i=0; i<cnt; i++) {
            CI_H0block_->buf_member[buf][i] = tmparr[i];
            }
         } /* end loop over bufs */
      } /* end icore==0 */

   free(tmparr);
}



/*
** CIvect::h0block_buf_ols()
**
** Parameters:
**    norm   = block's norm
**    ovrlap = blocks' overlap
*/
void CIvect::h0block_buf_ols(double *nx, double *ox, double *c1norm,double E_est)
{
   int i, j, k, blk, al, bl;
   double c, cn, tval, c1;

   for (i=0; i<CI_H0block_->buf_num[cur_buf_]; i++) {
      j = CI_H0block_->buf_member[cur_buf_][i];
      blk = CI_H0block_->blknum[j];
      al = CI_H0block_->alpidx[j];
      bl = CI_H0block_->betidx[j];
      c = CI_H0block_->c0b[j];
      cn = blocks_[blk][al][bl];
      c1 = cn - c;
      *nx -= cn * cn;
      *ox -= cn * c;
      *c1norm -= c1 * c1;
      tval = c + E_est * CI_H0block_->c0bp[j];
      tval -= CI_H0block_->s0bp[j];
      blocks_[blk][al][bl] = tval;
      /* CI_H0block_->c0b[j] = tval; */ /* this should gather all of c0b. Norm later */
                              /* What about the effect of symmetrization? */
                              /* symmetrization was the bug MLL */
     /*
      if (buf_offdiag_[cur_buf_]) {
         k = CI_H0block_->pair[j];
         if (k >= 0 && k != j) {
            CI_H0block_->c0b[k] = tval * phase;
            }
         }
      */
      *nx += tval * tval;
      *ox += tval * c;
      *c1norm += (tval - c) * (tval - c);
      }
}

/*
** CIvect::h0block_gather_vec(int vecode)
**
** Parameters:
**    curr = current vector number
**    vecode = 0 for C vector and 1 for Sigma vector
*/
void CIvect::h0block_gather_vec(int vecode)
{
   int buf, i, j, k, blk, al, bl;
   double c, cn, tval, phase, norm = 0.0;

   if (!CI_Params_->Ms0) phase = 1.0;
   else phase = ((int) CI_Params_->S % 2) ? -1.0 : 1.0;

   for (i=0; i<CI_H0block_->buf_num[cur_buf_]; i++) {
      j = CI_H0block_->buf_member[cur_buf_][i];
      blk = CI_H0block_->blknum[j];
      al = CI_H0block_->alpidx[j];
      bl = CI_H0block_->betidx[j];
      tval = blocks_[blk][al][bl];
      if (vecode) CI_H0block_->s0b[j] = tval;
      else CI_H0block_->c0b[j] = tval;
      if (buf_offdiag_[cur_buf_]) {
        k = CI_H0block_->pair[j];
        if (k >= 0 && k != j) {
        /* if (k >= 0 && k != j && CI_Params_->Ms0)  */
          if (vecode) CI_H0block_->s0b[k] = tval * phase;
          else CI_H0block_->c0b[k] = tval * phase;
           }
        }
      }
    /*
     if (!vecode) {
       outfile->Printf("c0b in h0block_gather_vec = \n");
       print_mat(&(CI_H0block_->c0b), 1, CI_H0block_->size, outfile);
       }
    */
}


/*
** CIvect::h0block_gather_multivec(double *vec)
**
** Parameters:
**    curr = current vector number
**    vecode = 0 for C vector and 1 for Sigma vector
*/
void CIvect::h0block_gather_multivec(double *vec)
{
   int buf, i, j, k, blk, al, bl;
   double c, cn, tval, phase, norm = 0.0;

   if (!CI_Params_->Ms0) phase = 1.0;
   else phase = ((int) CI_Params_->S % 2) ? -1.0 : 1.0;

   for (i=0; i<CI_H0block_->buf_num[cur_buf_]; i++) {
      j = CI_H0block_->buf_member[cur_buf_][i];
      blk = CI_H0block_->blknum[j];
      al = CI_H0block_->alpidx[j];
      bl = CI_H0block_->betidx[j];
      tval = blocks_[blk][al][bl];
      vec[j] = tval;
      if (buf_offdiag_[cur_buf_]) {
        k = CI_H0block_->pair[j];
        if (k >= 0 && k != j) {
        /* if (k >= 0 && k != j && CI_Params_->Ms0)  */
          vec[k] = tval * phase;
           }
        }
      }
}


/*
** CIvect::h0block_buf_precon(double *nx, int root)
**
** Routine used by sem
** Parameters:
**    norm   = block's norm
*/
void CIvect::h0block_buf_precon(double *nx, int root)
{
   int i, j, k, blk, al, bl, buf;
   double c, cn, tval, phase, norm = 0.0;

   if (!CI_Params_->Ms0) phase = 1.0;
   else phase = ((int) CI_Params_->S % 2) ? -1.0 : 1.0;

   for (buf=0; buf<buf_per_vect_; buf++) {
      read(root,buf);
      for (i=0; i<CI_H0block_->buf_num[buf]; i++) {
         j = CI_H0block_->buf_member[buf][i];
         blk = CI_H0block_->blknum[j];
         al = CI_H0block_->alpidx[j];
         bl = CI_H0block_->betidx[j];
         tval = blocks_[blk][al][bl] * blocks_[blk][al][bl];
         *nx -= tval;
         if (buf_offdiag_[buf]) {
           k = CI_H0block_->pair[j];
           if (k >= 0 && k!=j) *nx -= tval * phase;
           }
         tval = CI_H0block_->c0bp[j] * CI_H0block_->c0bp[j];
         *nx += tval;
         if (buf_offdiag_[buf]) {
           k = CI_H0block_->pair[j];
           if (k>= 0 && k!=j) *nx += tval * phase;
           }
         blocks_[blk][al][bl] = -CI_H0block_->c0bp[j];
         }
      write(root,buf);
      }
}


/*
** CALC_SSQ: Computes S^2
**
*/

double CIvect::calc_ssq(double *buffer1, double *buffer2,
             struct stringwr **alplist, struct stringwr **betlist, int vec_num)
{
   int bra_block, ket_block;
   int ket_ac, ket_bc, ket_nas, ket_nbs;
   int bra_ac, bra_bc, bra_nas, bra_nbs;
   int ket_birr, bra_birr;
   double tval = 0.0;
   double tval2 = 0.0;
   double S2, Ms;
   int i;

   buf_lock(buffer1);
   read(vec_num, 0);

   if (print_lvl_ > 4) {
     for (i=0; i<num_blocks_; i++) {
        ket_nas = Ia_size_[i];
        ket_nbs = Ib_size_[i];
        if (ket_nas==0 || ket_nbs==0) continue;
        print_mat(blocks_[i], ket_nas, ket_nbs, "outfile");
     }
   }

   /* loop over ket blocks of c */
   for (ket_block=0; ket_block<num_blocks_; ket_block++) {
      ket_ac = Ia_code_[ket_block];
      ket_bc = Ib_code_[ket_block];
      ket_nas = Ia_size_[ket_block];
      ket_nbs = Ib_size_[ket_block];
      if (ket_nas==0 || ket_nbs==0) continue;

      // CDS help: Is this correct?
      //ket_birr = ket_block / BetaG->subgr_per_irrep;
      ket_birr = ket_block / codes_per_irrep_;

      for (bra_block=0; bra_block<num_blocks_; bra_block++) {
         bra_ac = Ia_code_[bra_block];
         bra_bc = Ib_code_[bra_block];
         bra_nas = Ia_size_[bra_block];
         bra_nbs = Ib_size_[bra_block];
         if (bra_nas==0 || bra_nbs==0) continue;

         // CDS help: Is this correct?
         //bra_birr = bra_bc / BetaG->subgr_per_irrep;
         bra_birr = bra_bc / codes_per_irrep_;

         tval2 = ssq(alplist[ket_ac], betlist[ket_bc], blocks_[bra_block],
                   blocks_[ket_block], ket_nas, ket_nbs, bra_ac, bra_bc);
         tval += tval2;
         if (print_lvl_ > 4) {
           outfile->Printf("\nbra_block = %d\n",bra_block);
           outfile->Printf("ket_block = %d\n",ket_block);
           outfile->Printf("Contribution to <S_S+> = %lf\n",tval2);
         }
       } /* end loop over bra_blocks */

    } /* end loop over ket_block */

    Ms = 0.5 * (CI_CalcInfo_->num_alp_expl - CI_CalcInfo_->num_bet_expl);
    if (print_lvl_ > 1) {
      outfile->Printf("\n\n<S_z> = %lf\n", Ms);
      outfile->Printf("<S_z>^2 = %lf\n", Ms*Ms);
      outfile->Printf("<S_S+> = %lf\n", tval);
    }
    S2 = CI_CalcInfo_->num_bet_expl + tval + Ms + Ms*Ms;

    if (print_lvl_) outfile->Printf("Computed <S^2> vector %d = %20.15f\n\n", vec_num, S2);

  buf_unlock();
  return(S2);
}

int CIvect::check_zero_block(int blocknum)
{
   if (blocknum < 0 || blocknum > num_blocks_) {
      outfile->Printf( "CIvect::check_zero_block(): Block %d out of range\n",
              blocknum);
   }

   return(zero_blocks_[blocknum]);
}

void CIvect::set_zero_block(int blocknum, int value)
{
   if (blocknum < 0 || blocknum > num_blocks_) {
      outfile->Printf( "CIvect::set_zero_block(): Block %d out of range\n",
              blocknum);
   }

   if (value != 0 && value != 1) {
      outfile->Printf( "CIvect::set_zero_block(): Value %d out of range\n",
              value);
   }

   zero_blocks_[blocknum] = value;
}

void CIvect::set_zero_blocks_all(void)
{
   int i;

   for (i=0; i<num_blocks_; i++) zero_blocks_[i] = 1;

}

void CIvect::copy_zero_blocks(CIvect &src)
{
   int i;

   for (i=0; i<num_blocks_; i++) {
      zero_blocks_[i] = src.zero_blocks_[i];
      // outfile->Printf( "zero_block[%d] = %d\n", i, zero_blocks[i]);
   }
}

void CIvect::print_zero_blocks(void)
{
   int i;

   for (i=0; i<num_blocks_; i++) {
      outfile->Printf( "zero_block[%d] = %d\n", i, zero_blocks_[i]);
   }

}


/*
**
** scale_sigma()
**
** Modifies a sigma vector to account for scaling in
** perturbed Hamiltonian.
**
** Parameters:
**    Hd         = CIvect for H0
**    S          = CIvect for Sigma vector
**    C          = CIvect for Cvec vector
**    alplist    = alpha string list
**    betlist    = beta  string list
**    buf1       = first buffer for a CIvect
**    buf2       = second buffer for a CIvect
**
** Returns: none
*/
void CIvect::scale_sigma(CIvect &Hd, CIvect &C,
       struct stringwr **alplist, struct stringwr **betlist, int i,
       double *buf1, double *buf2)
{
   int buf;

   for (buf=0; buf<buf_per_vect_; buf++) {
      /* outfile->Printf(" i = %d\n", i);
      outfile->Printf("In scale_sigma\n"); */
      Hd.buf_lock(buf1);
      Hd.diag_mat_els_otf(alplist, betlist, CI_CalcInfo_->onel_ints->pointer(),
         CI_CalcInfo_->twoel_ints->pointer(), CI_CalcInfo_->e0_drc, CI_CalcInfo_->num_alp_expl,
         CI_CalcInfo_->num_bet_expl, CI_CalcInfo_->nmo, buf, ORB_ENER);
      C.buf_lock(buf2);
      C.read(i, buf);
      xexy(buf1, buf2, C.buf_size_[buf]);
      C.buf_unlock();
      buf_lock(buf2);
      read(i, buf);
      xexmy(buf2, buf1, buf_size_[buf]);
      xpeay(buf1, CI_Params_->perturbation_parameter, buf2, buf_size_[buf]);
      buf_unlock();
      Hd.buf_unlock();
      buf_lock(buf1);
      write(i, buf);
      buf_unlock();
      }

}

/*
**
** CIvect::dcalc_evangelisti()
**
** Function calculates the denominator part of the Davidson correction
** vector 'd'.
**
** Parameters:
**    rootnum   = number of current root
**    lambda    = current iteration's energy eigenvalue for the current root
**    Hd        = CIvector for the Hamiltonian diagonal
**    C         = CI vector
**
** Returns: sum of squares of coefficients in d.
*/
double CIvect::dcalc_evangelisti(int rootnum, int num_vecs, double lambda,
        CIvect &Hd, CIvect &C, double *buf1, double *buf2, int precon, int L,
        struct stringwr **alplist, struct stringwr **betlist, double **alpha)
{
   int buf, errcod, i;
   double tval, norm = 0.0;

   for (buf=0; buf<buf_per_vect_; buf++) {
      Hd.buf_unlock();

      buf_unlock();
      zero_arr(buf1, buf_size_[buf]);
      C.buf_lock(buf2);
      for (i=0; i<L; i++) {
         C.read(i, buf);
         xpeay(buf1, alpha[rootnum][i], buf2, C.buf_size_[buf]);
         }
      C.buf_unlock();
      buf_lock(buf2);
      read(rootnum, buf);

      xexy(buf2, buf1, buf_size_[buf]); /* r_I*c_I */
      xeax(buf2, -2.0, buf_size_[buf]); /* -2*r_I*c_I */
      xexy(buf1, buf1, buf_size_[buf]); /* c_I*c_I */
      xpey(buf1, buf2, buf_size_[buf]); /* -2*r_I*c_I + c_I*c_I */
      buf_unlock();
      Hd.buf_lock(buf2);
      if (CI_Params_->hd_otf == FALSE) Hd.read(0, buf);
      else if (CI_Params_->hd_otf == TRUE) {
          Hd.diag_mat_els_otf(alplist, betlist, CI_CalcInfo_->onel_ints->pointer(),
             CI_CalcInfo_->twoel_ints->pointer(), CI_CalcInfo_->edrc, CI_CalcInfo_->num_alp_expl,
             CI_CalcInfo_->num_bet_expl, CI_CalcInfo_->nmo, buf, CI_Params_->hd_ave);
        }
      xpey(buf2, buf1, buf_size_[buf]); /* Hd -2*r_I*c_I + c_I*c_I */
      buf_lock(buf1);
      read(rootnum, buf);
      tval = calc_d2(buf1, lambda, buf2, buf_size_[buf], precon);
      if (buf_offdiag_[buf]) tval *= 2.0;
      norm += tval;
      write(rootnum, buf);
      }

  return(norm);
}

/*
** Write the number of the new first buffer to disk.
** The new first buffer is the buffer which is renumbered as "zero" after
** a collapse of the subspace.  The labels on disk are not actually
** changed so that it is a little easier to deal with two logical CIvectors
** which point to the same physical CIvector (as happens if nodfile).
*/
void CIvect::write_new_first_buf(void)
{
  int unit;

  unit = first_unit_;
  psio_write_entry((ULI) unit, "New First Buffer", (char *) &new_first_buf_,
    sizeof(int));
}

/*
** Read the number of the new first buffer from disk.
** The new first buffer is the buffer which is renumbered as "zero" after
** a collapse of the subspace.  The labels on disk are not actually
** changed so that it is a little easier to deal with two logical CIvectors
** which point to the same physical CIvector (as happens if nodfile).
** Return -1 if "New First Buffer" is not stored in the file yet.
*/
int CIvect::read_new_first_buf(void)
{
  int unit;
  int nfb;

  unit = first_unit_;
  if (psio_tocscan((ULI) unit, "New First Buffer") == NULL) return(-1);
  psio_read_entry((ULI) unit, "New First Buffer", (char *) &nfb,
    sizeof(int));
  return(nfb);

}

/*
** Set the number of the new first buffer.
** The new first buffer is the buffer which is renumbered as "zero" after
** a collapse of the subspace.  The labels on disk are not actually
** changed so that it is a little easier to deal with two logical CIvectors
** which point to the same physical CIvector (as happens if nodfile).
*/
void CIvect::set_new_first_buf(int nfb)
{
  new_first_buf_ = nfb;
}


/*
** Read the number of valid vectors in this object.  That will be stored
** in the first unit.
*/
int CIvect::read_num_vecs(void)
{
  int unit;
  int nv;

  unit = first_unit_;
  if (psio_tocscan((ULI) unit, "Num Vectors") == NULL) return(-1);
  psio_read_entry((ULI) unit, "Num Vectors", (char *) &nv, sizeof(int));
  return(nv);
}


/*
** Write the number of valid vectors in this object.  That will be stored
** in the first unit.
*/
void CIvect::write_num_vecs(int nv)
{
  int unit;

  unit = first_unit_;
  psio_write_entry((ULI) unit, "Num Vectors", (char *) &nv, sizeof(int));
  write_toc();
  //civect_psio_debug();
}


/*
** Write the libpsio table of contents to disk in case we crash before
** we're done.  The TOC is written to the end of the file.  If we aren't
** done filling it up, it will be wiped out by the next write but written
** again at the new end of file next time we call this function.
*/
void CIvect::write_toc(void)
{
  int i,unit;

  for (i=0; i<nunits_; i++) {
    psio_tocwrite(units_[i]);
  }

}


/*
** Print libpsio debug info
*/
void CIvect::civect_psio_debug(void)
{
  int i, unit;

  /* psio_tocprint not available right now in PSI4; re-enable it if you
     need this functionality.
  */
  /*
  for (i=0; i<nunits_; i++)
    psio_tocprint(units_[i], outfile);
  */
  outfile->Printf( "Number of vectors = %d\n", read_num_vecs());
  outfile->Printf( "New first buffer_ = %d\n", read_new_first_buf());
  outfile->Printf( "Internal new first buffer_ = %d\n", new_first_buf_);
}


/*
** perturbation theory correction
** CDS 2/04
*/
void CIvect::pt_correction(struct stringwr **alplist, struct stringwr
      **betlist)
{
  int block, iac, ibc;
  int nas, nbs;

  if (icore_ == 1) { /* whole vector at once */
    for (block=0; block<num_blocks_; block++) {
      iac = Ia_code_[block];  nas = Ia_size_[block];
      ibc = Ib_code_[block];  nbs = Ib_size_[block];
    }
  }
  else {
    outfile->Printf( "only icore_=1 works for now\n");
  }

}

/*
** CIvect::compute_follow_overlap
**
** Computes the overlap with some user-supplied vector
** Only works for icore==1 (whole vector) for now
**
** Returns: the overlap
*/
double CIvect::compute_follow_overlap(int troot, int ncoef, double *coef,
  int *Iac, int *Iaridx, int *Ibc, int *Ibridx)
{
  int i, a, b, blk;
  double tval;

  if (icore_ != 1) {
    outfile->Printf( "CIvect::compute_follow_overlap: can't use icore != 1\n");
    return(0.0);
  }

  read(troot,0);

  tval = 0.0;

  for (i=0; i<ncoef; i++) {
    blk = decode_[Iac[i]][Ibc[i]];
    a = Iaridx[i];
    b = Ibridx[i];
    tval += blocks_[blk][a][b] * coef[i];
  }

  tval = fabs(tval);
  return(tval);

}




/*
** calc_hd_block(): Function calculates a block of H0, the diagonal elements of
**    the Hamiltonian matrix.
**
** Parameters:
**    alplist_local = list of alpha strings with replacements (used to get occs)
**    betlist_local = list of beta strings with replacements
**    nas     = number of alpha strings in list
**    nbs     = number of beta strings in list
**    H0      = matrix to hold results (stored as H0[alpidx][betidx])
**    oei     = one-electron integrals
**    tei     = two-electron integrals
**    na      = number of explicit alpha electrons
**    nb      = number of explicit beta electrons
**    nbf     = number of orbitals in CI
**    edrc    = energy of the dropped core orbitals
**
*/
void CIvect::calc_hd_block(struct stringwr *alplist_local, struct stringwr *betlist_local,
      double **H0, double *oei, double *tei, double edrc,
      int nas, int nbs, int na, int nb, int nbf)
{
   int acnt, bcnt;
   int a1, a2, b1, b2;
   int i,j, ii, iii, jj, ij, iijj, ijij;
   double value;
   struct stringwr *betlist0;

   betlist0 = betlist_local;

   for (acnt=0; acnt<nas; acnt++) {

      for (bcnt=0, betlist_local=betlist0; bcnt<nbs; bcnt++) {

         /* add dropped core energy first */
/************************************************/
         value = edrc;

         /* loop over alpha occs */
         for (a1=0; a1<na; a1++) {
            i = (int) alplist_local->occs[a1];
            ii = ioff[i] + i;
            value += oei[ii];
            /* outfile->Printf("oei[%d] = %lf\n",ii,oei[ii]); */
            iii = ioff[ii];

            for (a2=0; a2<a1; a2++) {
               j = (int) alplist_local->occs[a2];
               jj = ioff[j] + j;
               iijj = iii + jj;
               ij = ioff[i] + j;
               ijij = ioff[ij] + ij;
               value += tei[iijj] - tei[ijij];
               }

            for (b1=0; b1<nb; b1++) {
               j = (int) betlist_local->occs[b1];
               jj = ioff[j] + j;
               iijj = ioff[MAX0(ii,jj)] + MIN0(ii,jj);
               value += tei[iijj];
               }
           }

         for (b1=0; b1<nb; b1++) {
            i = (int) betlist_local->occs[b1];
            ii = ioff[i] + i;
            value += oei[ii];
            iii = ioff[ii];

            for (b2=0; b2<b1; b2++) {
               j = (int) betlist_local->occs[b2];
               jj = ioff[j] + j;
               iijj = iii + jj;
               ij = ioff[i] + j;
               ijij = ioff[ij] + ij;
               value += tei[iijj] - tei[ijij];
               }
            }

         H0[acnt][bcnt] = value;
      /*
         outfile->Printf("H0[%d][%d] = %lf\n",acnt,bcnt,value);
      */
         betlist_local++;
         } /* end loop over bcnt */

      alplist_local++;
      }

}
/*
** calc_hd_block_ave(): Function calculates a block of H0 and the diagonal elements
** of the Hamiltonian matrix averaged over spin-coupling sets to correct any
** spin contamination of the c and sigma vectors.
**
** Parameters:
**    alplist_local = list of alpha strings with replacements (used to get occs)
**    betlist_local = list of beta strings with replacements
**    nas     = number of alpha strings in list
**    nbs     = number of beta strings in list
**    H0      = matrix to hold results (stored as H0[alpidx][betidx])
**    oei     = one-electron integrals
**    tei     = two-electron integrals
**    na      = number of explicit alpha electrons
**    nb      = number of explicit beta electrons
**    nbf     = number of orbitals in CI
**    edrc    = energy of the dropped core orbitals
**
*/
void CIvect::calc_hd_block_ave(struct stringwr *alplist_local, struct stringwr *betlist_local,
      double **H0, double *tf_oei, double *tei, double edrc,
      int nas, int nbs, int na, int nb, int nbf)
{
   int acnt, bcnt;
   int a1, a2, a3, b1, b2, b3;
   int i,j, ii, iii, jj, ij, iijj, ijij;
   double value, tval = 0.0, tval2, Kave;
   struct stringwr *betlist0;
   double k_total; /* total number of K ints in energy expression */
   int k_combo; /* total combination of unique K ints over spin-coupling set */
   int *unique_occs; /* the uniquely occupied orbitals for a given determinant */
   int num_el;  /* total number of electrons explicitly treated */
   int num_unique; /* number of unique orbitals */
   betlist0 = betlist_local;


   k_total = combinations(na,2) + combinations(nb,2);

   num_el = na + nb;
   unique_occs = init_int_array(num_el);

   for (acnt=0; acnt<nas; acnt++) {

      for (bcnt=0, betlist_local=betlist0; bcnt<nbs; bcnt++) {

         /* add dropped core energy first */
         value = edrc;

         /* loop over alpha occs */
         for (a1=0; a1<na; a1++) {
            i = (int) alplist_local->occs[a1];
            ii = ioff[i] + i;
            /* h_ii bar alpha alpha */
            value += tf_oei[ii];
            /* outfile->Printf("tf_oei[%d] = %lf\n",ii,tf_oei[ii]); */
            iii = ioff[ii];

            /* loop over alpha occs */
            for (a2=0; a2<a1; a2++) {
               j = (int) alplist_local->occs[a2];
               jj = ioff[j] + j;
               iijj = iii + jj;
               /* J alpha alpha */
               value += tei[iijj];
               }

            /* loop over beta occs */
            for (b1=0; b1<nb; b1++) {
               j = (int) betlist_local->occs[b1];
               jj = ioff[j] + j;
               iijj = ioff[MAX0(ii,jj)] + MIN0(ii,jj);
               value += tei[iijj];
               }
           }

         /* loop over beta occs */
         for (b1=0; b1<nb; b1++) {
            i = (int) betlist_local->occs[b1];
            ii = ioff[i] + i;
            value += tf_oei[ii];
            /* outfile->Printf("tf_oei[%d] = %lf\n",ii,tf_oei[ii]); */
            iii = ioff[ii];

            /* loop over beta occs */
            for (b2=0; b2<b1; b2++) {
               j = (int) betlist_local->occs[b2];
               jj = ioff[j] + j;
               iijj = iii + jj;
               ij = ioff[i] + j;
               ijij = ioff[ij] + ij;
               value += tei[iijj];
               }
            }

         /* determine average K over spin-coupling set */
         num_unique = 0;
         for (a1=0; a1<na; a1++) unique_occs[num_unique++] = (int) alplist_local->occs[a1];
         /* for (j=0; j<num_unique; j++)
            outfile->Printf("unique_occs[%d] = %d\n",j,unique_occs[j]); */
            for (b1=0; b1<nb; b1++) {
               j = (int) betlist_local->occs[b1];
               for (a1=0; a1<na; a1++) {
                  if (j==unique_occs[a1]) break;
                  if (a1==(na-1)) unique_occs[num_unique++] = j;
                  }
               }
         /* outfile->Printf("num_unique = %d\n",num_unique);
         outfile->Printf("num_el = %d\n",num_el);
         */
         if (num_unique>num_el) outfile->Printf("WARNING: The number of explicit electrons" \
                             "!= num_el\n");

       /*
         for (j=0; j<na; j++)
            outfile->Printf("alp_occs[%d] = %d\n",j,(int)alplist_local->occs[j]);
         for (j=0; j<nb; j++)
            outfile->Printf("bet_occs[%d] = %d\n",j,(int)betlist_local->occs[j]);
         for (j=0; j<num_unique; j++)
            outfile->Printf("unique_occs[%d] = %d\n",j,unique_occs[j]);
       */

         Kave = 0.0;
         for (a1=0; a1<num_unique; a1++) {
            i = unique_occs[a1];
            for (b1=0; b1<a1; b1++) {
               j = unique_occs[b1];
               ij = ioff[MAX0(i,j)] + MIN0(i,j);
               ijij = ioff[ij] + ij;
               Kave += tei[ijij];
               /* outfile->Printf("tei[%d] = %lf\n",ijij,tei[ijij]); */
               }
            }

         /* outfile->Printf("num_unique = %d\n",num_unique);
         outfile->Printf("ioff[num_unique-1] = %d\n",ioff[num_unique]);
         outfile->Printf("k_total = %d\n",k_total);
         */

         if (num_unique > 1) Kave /= ioff[num_unique-1];
         value -= 0.5 * Kave * k_total;
         /* outfile->Printf("Kave = %lf\n",Kave); */

         if (print_lvl_ > 5) {
           outfile->Printf("acnt = %d\t bcnt = %d\n",acnt,bcnt);
           outfile->Printf("tval = %lf\n",tval);
           for(a1=0; a1<na; a1++)
             outfile->Printf(" %d",alplist_local->occs[a1]);
           outfile->Printf(" \n");
           for(b1=0; b1<nb; b1++)
             outfile->Printf(" %d",betlist_local->occs[b1]);
           outfile->Printf(" \n");
           }

         H0[acnt][bcnt] = value;
         /* outfile->Printf("H0[%d][%d] = %lf\n",acnt,bcnt,value); */
         betlist_local++;
         } /* end loop over bcnt */

      alplist_local++;
      }

}

/*
** calc_hd_block_orbenergy(): Function calculates a block of H0 and the diagonal elements
** of the Hamiltonian matrix as the sum of orbital energies.
**
** Parameters:
**    alplist_local = list of alpha strings with replacements (used to get occs)
**    betlist_local = list of beta strings with replacements
**    nas     = number of alpha strings in list
**    nbs     = number of beta strings in list
**    H0      = matrix to hold results (stored as H0[alpidx][betidx])
**    oei     = one-electron integrals
**    tei     = two-electron integrals
**    na      = number of explicit alpha electrons
**    nb      = number of explicit beta electrons
**    nbf     = number of orbitals in CI
**    edrc    = energy of the dropped core orbitals
**
*/
void CIvect::calc_hd_block_orbenergy(struct stringwr *alplist_local,
      struct stringwr *betlist_local, double **H0, double *oei,
      double *tei, double edrc, int nas, int nbs, int na, int nb, int nbf)
{
   int acnt, bcnt;
   int a1, b1, i,j;
   double value, tval;
   struct stringwr *betlist0, *alplist0;
   double *orb_e_diff_alp, *orb_e_diff_bet;
   double sum_orb_energies = 0.0;

   betlist0 = betlist_local;
   alplist0 = alplist_local;

   orb_e_diff_alp = init_array(nas);
   orb_e_diff_bet = init_array(nbs);
  /* if (CI_Params_->Ms0) orb_e_diff_bet = &orb_e_diff_alp;
   else orb_e_diff_bet = init_array(CI_CalcInfo_->num_bet_str);
  */

   for (acnt=0; acnt<nas; acnt++) {
      orb_e_diff_alp[acnt] = 0.0;
      for (a1=0; a1<na; a1++) {
         i = (int) alplist_local->occs[a1];
         i += CI_CalcInfo_->num_drc_orbs;
         if(CI_Params_->zaptn)
           orb_e_diff_alp[acnt] += CI_CalcInfo_->scfeigvala[i];
         else
           orb_e_diff_alp[acnt] += CI_CalcInfo_->scfeigval[i];
         }
      alplist_local++;
      }

   for (bcnt=0; bcnt<nbs; bcnt++) {
      orb_e_diff_bet[bcnt] = 0.0;
      for (b1=0; b1<nb; b1++) {
         j = (int) betlist_local->occs[b1];
         j += CI_CalcInfo_->num_drc_orbs;
         if(CI_Params_->zaptn)
           orb_e_diff_bet[bcnt] += CI_CalcInfo_->scfeigvalb[j];
         else
           orb_e_diff_bet[bcnt] += CI_CalcInfo_->scfeigval[j];
         }
      betlist_local++;
      }

   alplist_local = alplist0;
   betlist_local = betlist0;

   for (acnt=0; acnt<nas; acnt++) {
         tval = edrc + orb_e_diff_alp[acnt];
      for (bcnt=0; bcnt<nbs; bcnt++) {
         value = orb_e_diff_bet[bcnt] + tval;
         H0[acnt][bcnt] = value;
        /*
         outfile->Printf("H0[%d][%d] = %lf\n",acnt,bcnt,value);
        */
         betlist_local++;
         } /* end loop over bcnt */
      alplist_local++;
      }

/* Free up memory */
free(orb_e_diff_alp);
free(orb_e_diff_bet);


}

/*
** calc_hd_block_evangelisti(): Function calculates a block of H0 and the diagonal elements
** of the Hamiltonian matrix averaged over spin-coupling sets to correct any
** spin contamination of the c and sigma vectors.
**
** Parameters:
**    alplist_local = list of alpha strings with replacements (used to get occs)
**    betlist_local = list of beta strings with replacements
**    nas     = number of alpha strings in list
**    nbs     = number of beta strings in list
**    H0      = matrix to hold results (stored as H0[alpidx][betidx])
**    oei     = one-electron integrals
**    tei     = two-electron integrals
**    na      = number of explicit alpha electrons
**    nb      = number of explicit beta electrons
**    nbf     = number of orbitals in CI
**    edrc    = energy of the dropped core orbitals
**
*/
void CIvect::calc_hd_block_evangelisti(struct stringwr **alplist, struct stringwr **betlist,
     struct stringwr *alplist_local, struct stringwr *betlist_local,
      double **H0, double *tf_oei, double *tei, double edrc,
      int nas, int nbs, int na, int nb, int nbf)
{
   int acnt, bcnt;
   int a1, b1, i,j;
   double value, tval;
   struct stringwr *betlist0, *alplist0;
   double *orb_e_diff_alp, *orb_e_diff_bet;
   int num_alp_diff, num_bet_diff;
   int **orb_diff, *jnk;
   int sign;

   betlist0 = betlist_local;
   alplist0 = alplist_local;

   orb_diff = init_int_matrix(2,na);
   jnk = init_int_array(na);
   orb_e_diff_alp = init_array(nas);
   orb_e_diff_bet = init_array(nbs);

   for (acnt=0; acnt<nas; acnt++) {
      orb_e_diff_alp[acnt] = 0.0;
      num_alp_diff = calc_orb_diff(na,
                     alplist[CI_CalcInfo_->ref_alp_list][CI_CalcInfo_->ref_alp_rel].occs,
                     alplist_local->occs, orb_diff[0], orb_diff[1], &sign,
                     jnk, 1);
      for (a1=0; a1<num_alp_diff; a1++) {
         i = orb_diff[0][a1];
         j = orb_diff[1][a1];
         i += CI_CalcInfo_->num_drc_orbs;
         j += CI_CalcInfo_->num_drc_orbs;
         orb_e_diff_alp[acnt] += CI_CalcInfo_->scfeigval[j]
                                 - CI_CalcInfo_->scfeigval[i];
         }
      alplist_local++;
      }

   for (bcnt=0; bcnt<nbs; bcnt++) {
      orb_e_diff_bet[bcnt] = 0.0;
      num_bet_diff = calc_orb_diff(nb,
                     betlist[CI_CalcInfo_->ref_bet_list][CI_CalcInfo_->ref_bet_rel].occs,
                     betlist_local->occs, orb_diff[0], orb_diff[1], &sign,
                     jnk, 1);
      for (b1=0; b1<num_bet_diff; b1++) {
         i = orb_diff[0][b1];
         j = orb_diff[1][b1];
         i += CI_CalcInfo_->num_drc_orbs;
         j += CI_CalcInfo_->num_drc_orbs;
         orb_e_diff_bet[bcnt] += CI_CalcInfo_->scfeigval[j]
                                 - CI_CalcInfo_->scfeigval[i];
         }
      betlist_local++;
      }

   alplist_local = alplist0;
   betlist_local = betlist0;

   for (acnt=0; acnt<nas; acnt++) {
         /* add dropped core energy first */
         tval = CI_CalcInfo_->escf - CI_CalcInfo_->enuc;
         tval += orb_e_diff_alp[acnt];
      for (bcnt=0; bcnt<nbs; bcnt++) {
         value = 0.0;
         value = orb_e_diff_bet[bcnt] + tval;
         H0[acnt][bcnt] = value;
         /* outfile->Printf("H0[%d][%d] = %lf\n",acnt,bcnt,value); */
         betlist_local++;
         } /* end loop over bcnt */
      alplist_local++;
      }

/* Free memory */
/*
free(jnk);
free(orb_e_diff_alp);
free(orb_e_diff_bet);
free(orb_diff);
*/


}


/*
** calc_hd_block_mll(): Function calculates a block of H0 and the diagonal elements
** of the Hamiltonian matrix as the sum of orbital energies.
**
** Parameters:
**    alplist_local = list of alpha strings with replacements (used to get occs)
**    betlist_local = list of beta strings with replacements
**    nas     = number of alpha strings in list
**    nbs     = number of beta strings in list
**    H0      = matrix to hold results (stored as H0[alpidx][betidx])
**    oei     = one-electron integrals
**    tei     = two-electron integrals
**    na      = number of explicit alpha electrons
**    nb      = number of explicit beta electrons
**    nbf     = number of orbitals in CI
**    edrc    = energy of the dropped core orbitals
**
*/
void CIvect::calc_hd_block_mll(struct stringwr *alplist_local,
      struct stringwr *betlist_local, double **H0, double *oei,
      double *tei, double edrc, int nas, int nbs, int na, int nb, int nbf)
{
   int acnt, bcnt;
   int a1, b1, i,j, i_offset, j_offset, ii, jj;
   double value, tval;
   struct stringwr *betlist0, *alplist0;
   double *orb_e_diff_alp, *orb_e_diff_bet;
   double *oei_alp, *oei_bet, *eigval;

   betlist0 = betlist_local;
   alplist0 = alplist_local;

   oei_alp = init_array(nas);
   oei_bet = init_array(nbs);
   orb_e_diff_alp = init_array(nas);
   orb_e_diff_bet = init_array(nbs);
  /* if (CI_Params_->Ms0) orb_e_diff_bet = &orb_e_diff_alp;
   else orb_e_diff_bet = init_array(nbs);
  */

   for (acnt=0; acnt<nas; acnt++) {
      orb_e_diff_alp[acnt] = oei_alp[acnt] = 0.0;
      for (a1=0; a1<na; a1++) {
         i = (int) alplist_local->occs[a1];
         ii = ioff[i] + i;
         i_offset = i + CI_CalcInfo_->num_drc_orbs;
         oei_alp[acnt] += oei[ii];
         orb_e_diff_alp[acnt] += CI_CalcInfo_->scfeigval[i_offset] - oei[ii];
         }
      alplist_local++;
      }

   for (bcnt=0; bcnt<nbs; bcnt++) {
      orb_e_diff_bet[bcnt] = oei_bet[bcnt] = 0.0;
      for (b1=0; b1<nb; b1++) {
         j = (int) betlist_local->occs[b1];
         jj = ioff[j] + j;
         j_offset = j + CI_CalcInfo_->num_drc_orbs;
         oei_bet[bcnt] += oei[jj];
         orb_e_diff_bet[bcnt] += CI_CalcInfo_->scfeigval[j_offset] - oei[jj];
         }
      betlist_local++;
      }

   alplist_local = alplist0;
   betlist_local = betlist0;

   for (acnt=0; acnt<nas; acnt++) {
         tval = edrc + 0.5 * orb_e_diff_alp[acnt] + oei_alp[acnt];
      for (bcnt=0; bcnt<nbs; bcnt++) {
         value = 0.5 * orb_e_diff_bet[bcnt] + oei_bet[bcnt] + tval;
         H0[acnt][bcnt] = value;
         betlist_local++;
         } /* end loop over bcnt */
      alplist_local++;
      }

 free(oei_alp);
 free(oei_bet);
 free(orb_e_diff_alp);
 free(orb_e_diff_bet);
}

/*
** calc_hd_block_z_ave(): Function calculates a block of H0 and the diagonal
** elements of the Hamiltonian matrix averaged over spin-coupling sets to
** correct any spin contamination of the c and sigma vectors.
**
** Parameters:
**    alplist_local = list of alpha strings with replacements (used to get occs)
**    betlist_local = list of beta strings with replacements
**    nas     = number of alpha strings in list
**    nbs     = number of beta strings in list
**    H0      = matrix to hold results (stored as H0[alpidx][betidx])
**    oei     = one-electron integrals
**    tei     = two-electron integrals
**    na      = number of explicit alpha electrons
**    nb      = number of explicit beta electrons
**    nbf     = number of orbitals in CI
**    edrc    = energy of the dropped core orbitals
**
*/
void CIvect::calc_hd_block_z_ave(struct stringwr *alplist_local,
  struct stringwr *betlist_local, double **H0, double pert_param,
  double *tei, double edrc, int nas, int nbs, int na, int nb, int nbf)
{
   int acnt, bcnt;
   int a1, a2, a3, b1, b2, b3;
   int i,j, ii, iii, jj, ij, iijj, ijij;
   double value, tval = 0.0, tval2, Kave;
   struct stringwr *betlist0;
   double k_total; /* total number of K ints in energy expression */
   int k_combo; /* total combination of unique K ints over spin-coupling set */
   int *unique_occs; /* the uniquely occupied orbitals for a given determinant */
   int num_el;  /* total number of electrons explicitly treated */
   int num_unique; /* number of unique orbitals */
   betlist0 = betlist_local;


   k_total = combinations(na,2) + combinations(nb,2);
   num_el = na + nb;
   unique_occs = init_int_array(num_el);

   for (acnt=0; acnt<nas; acnt++) {

      for (bcnt=0, betlist_local=betlist0; bcnt<nbs; bcnt++) {

         /* add dropped core energy first */
         value = edrc;

         /* loop over alpha occs */
         for (a1=0; a1<na; a1++) {
            i = (int) alplist_local->occs[a1];
            value += CI_CalcInfo_->scfeigval[i+CI_CalcInfo_->num_drc_orbs];
            ii = ioff[i] + i;
            /* h_ii bar alpha alpha */
            iii = ioff[ii];

            /* loop over alpha occs */
            for (a2=0; a2<a1; a2++) {
               j = (int) alplist_local->occs[a2];
               jj = ioff[j] + j;
               iijj = iii + jj;
               /* J alpha alpha */
               value -= pert_param * tei[iijj];
               }

            /* loop over beta occs */
            for (b1=0; b1<nb; b1++) {
               j = (int) betlist_local->occs[b1];
               jj = ioff[j] + j;
               iijj = ioff[MAX0(ii,jj)] + MIN0(ii,jj);
               value -= pert_param * tei[iijj];
               }
           }

         /* loop over beta occs */
         for (b1=0; b1<nb; b1++) {
            i = (int) betlist_local->occs[b1];
            value += CI_CalcInfo_->scfeigval[i+CI_CalcInfo_->num_drc_orbs];
            ii = ioff[i] + i;
            iii = ioff[ii];

            /* loop over beta occs */
            for (b2=0; b2<b1; b2++) {
               j = (int) betlist_local->occs[b2];
               jj = ioff[j] + j;
               iijj = iii + jj;
               ij = ioff[i] + j;
               ijij = ioff[ij] + ij;
               value -= pert_param * tei[iijj];
               }
            }

         /* determine average K over spin-coupling set */
         num_unique = 0;
         for (a1=0; a1<na; a1++) unique_occs[num_unique++] = (int) alplist_local->occs[a1];
         /* for (j=0; j<num_unique; j++)
            outfile->Printf("unique_occs[%d] = %d\n",j,unique_occs[j]); */
            for (b1=0; b1<nb; b1++) {
               j = (int) betlist_local->occs[b1];
               for (a1=0; a1<na; a1++) {
                  if (j==unique_occs[a1]) break;
                  if (a1==(na-1)) unique_occs[num_unique++] = j;
                  }
               }
         /* outfile->Printf("num_unique = %d\n",num_unique);
         outfile->Printf("num_el = %d\n",num_el);
         */
         if (num_unique>num_el) outfile->Printf("WARNING: The number of explicit electrons" \
                             "!= num_el\n");

       /*
         for (j=0; j<na; j++)
            outfile->Printf("alp_occs[%d] = %d\n",j,(int)alplist_local->occs[j]);
         for (j=0; j<nb; j++)
            outfile->Printf("bet_occs[%d] = %d\n",j,(int)betlist_local->occs[j]);
         for (j=0; j<num_unique; j++)
            outfile->Printf("unique_occs[%d] = %d\n",j,unique_occs[j]);
       */

         Kave = 0.0;
         for (a1=0; a1<num_unique; a1++) {
            i = unique_occs[a1];
            for (b1=0; b1<a1; b1++) {
               j = unique_occs[b1];
               ij = ioff[MAX0(i,j)] + MIN0(i,j);
               ijij = ioff[ij] + ij;
               Kave += tei[ijij];
               /* outfile->Printf("tei[%d] = %lf\n",ijij,tei[ijij]); */
               }
            }

         /* outfile->Printf("num_unique = %d\n",num_unique);
         outfile->Printf("ioff[num_unique-1] = %d\n",ioff[num_unique]);
         outfile->Printf("k_total = %d\n",k_total);
         */
         if (num_unique > 1) Kave /= ioff[num_unique-1];
         value += 0.5 * Kave * k_total * pert_param;
         /* outfile->Printf("Kave = %lf\n",Kave); */

         if (print_lvl_ > 5) {
           outfile->Printf("acnt = %d\t bcnt = %d\n",acnt,bcnt);
           outfile->Printf("tval = %lf\n",tval);
           for(a1=0; a1<na; a1++)
             outfile->Printf(" %d",alplist_local->occs[a1]);
           outfile->Printf(" \n");
           for(b1=0; b1<nb; b1++)
             outfile->Printf(" %d",betlist_local->occs[b1]);
           outfile->Printf(" \n");
           }
         H0[acnt][bcnt] = value;
       /*
         outfile->Printf("H0[%d][%d] = %lf\n",acnt,bcnt,value);
       */
         betlist_local++;
         } /* end loop over bcnt */

      alplist_local++;
      }

}
double CIvect::ssq(struct stringwr *alplist, struct stringwr *betlist,
     double **CL, double **CR, int nas, int nbs,
     int Ja_list, int Jb_list)
{
   struct stringwr *Ia, *Ib ;
   unsigned int Ia_ex, Ib_ex;
   int Ia_idx, Ib_idx;
   int Ja_idx, Jb_idx;
   int Ja_sgn, Jb_sgn;
   int ij, ji, i1, j1, i2, j2;
   double tval, Ms, S2, smin_spls = 0.0;

   int Iacnt, Jbcnt, *Iaij, *Ibij;
   unsigned int *Iaridx, *Ibridx;
   signed char *Iasgn, *Ibsgn;

   /* <S^2> = <S_z> + <S_z>^2 + <S_S+> */
   /* First determine the expection value of <S_S+> */

   /* loop over Ia */
   if (print_lvl_ > 2) {
     outfile->Printf("number of alpha strings = %d\n",nas);
   }
   for (Ia=alplist,Ia_idx=0; Ia_idx < nas; Ia_idx++,Ia++) {

      /* loop over excitations E^a_{ji} from |A(I_a)> */
      Iacnt = Ia->cnt[Ja_list];
      Iaridx = Ia->ridx[Ja_list];
      Iasgn = Ia->sgn[Ja_list];
      Iaij = Ia->oij[Ja_list];
      for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
         ji = *Iaij++;
         Ja_idx = *Iaridx++;
         Ja_sgn = *Iasgn++;
         i1 = ji/CI_CalcInfo_->num_ci_orbs;
         j1 = ji%CI_CalcInfo_->num_ci_orbs;

         /* loop over Ib */
         if (print_lvl_ > 2) {
           outfile->Printf("number of beta strings = %d\n",nbs);
         }
         for (Ib=betlist, Ib_idx=0; Ib_idx < nbs; Ib_idx++, Ib++) {

            /* loop over excitations E^b_{ij} from |B(I_b)> */
            Jbcnt = Ib->cnt[Jb_list];
            Ibridx = Ib->ridx[Jb_list];
            Ibsgn = Ib->sgn[Jb_list];
            Ibij = Ib->oij[Jb_list];

            tval = 0.0;
            for (Ib_ex=0; Ib_ex < Jbcnt; Ib_ex++) {
               ij = *Ibij++;
               Jb_idx = *Ibridx++;
               Jb_sgn = *Ibsgn++;
               i2 = ij/CI_CalcInfo_->num_ci_orbs;
               j2 = ij%CI_CalcInfo_->num_ci_orbs;
               if (i1!=j2 || i2!=j1) continue;
               tval += CR[Ia_idx][Ib_idx] * CL[Ja_idx][Jb_idx] *
                   (double) Ja_sgn * (double) Jb_sgn;
               if (print_lvl_ > 3) {
                 outfile->Printf("\n\nIa_idx = %d\n",Ia_idx);
                 outfile->Printf("Ib_idx = %d\n",Ib_idx);
                 outfile->Printf("Ja_idx = %d\n",Ja_idx);
                 outfile->Printf("Jb_idx = %d\n",Jb_idx);
                 outfile->Printf("tval_ssq = %lf\n",-tval);
                 outfile->Printf("CR = %lf\n",CR[Ia_idx][Ib_idx]);
                 outfile->Printf("LR = %lf\n",CL[Ja_idx][Jb_idx]);
                 outfile->Printf("Ja_sgn = %lf\n",Ja_sgn);
                 outfile->Printf("Jb_sgn = %lf\n",Jb_sgn);
               }
            }
            smin_spls += tval;

          } /* end loop over Ib */
       } /* end loop over Ia excitations */
     } /* end loop over Ia */

   S2 = -smin_spls;

   return(S2);
}

}} // namespace psi::detci
