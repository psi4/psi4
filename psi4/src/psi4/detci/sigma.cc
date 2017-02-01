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
**  \brief Routines to compute sigma = H * c
**
** Here collect the stuff to calculate the sigma vector within the
** framework of the CI vector class.
** Rewrote lots of stuff to handle the three out-of-core cases better.
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia
** February 1996
**
*/


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/vector.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/civect.h"
#include "psi4/detci/ciwave.h"

namespace psi { namespace detci {

extern void transp_sigma(double **a, int rows, int cols, int phase);
//extern void H0block_gather(double **mat, int al, int bl, int cscode,
//   int mscode, int phase);
extern void b2brepl(unsigned char **occs, int *Jcnt, int **Jij, int **Joij,
   int **Jridx, signed char **Jsgn, struct olsen_graph *Graph,
   int Ilist, int Jlist, int len, struct calcinfo *Cinfo);
extern void b2brepl_test(unsigned char ***occs, int *Jcnt, int **Jij,
   int **Joij, int **Jridx, signed char **Jsgn, struct olsen_graph *Graph,
   struct calcinfo *Cinfo);
extern void s3_block_bz(int Ialist, int Iblist, int Jalist,
   int Jblist, int nas, int nbs, int cnas,
   double *tei, double **C, double **S, double **Cprime, double **Sprime,
   struct calcinfo *CalcInfo, int ***OV);
extern void set_row_ptrs(int rows, int cols, double **matrix);

extern void s1_block_vfci(struct stringwr **alplist,
   struct stringwr **betlist,
   double **C, double **S, double *oei, double *tei, double *F,
   int nlists, int nas, int nbs, int Ib_list, int Jb_list,
   int Jb_list_nbs);
extern void s1_block_vras(struct stringwr **alplist,
   struct stringwr **betlist,
   double **C, double **S, double *oei, double *tei, double *F,
   int nlists, int nas, int nbs, int sbc, int cbc, int cnbs);
extern void s1_block_vras_rotf(int *Cnt[2], int **Ij[2], int **Oij[2],
   int **Ridx[2], signed char **Sgn[2], unsigned char **Toccs,
   double **C, double **S,
   double *oei, double *tei, double *F, int nlists, int nas, int nbs,
   int Ib_list, int Jb_list, int Jb_list_nbs, struct olsen_graph *BetaG,
   struct calcinfo *CIinfo, unsigned char ***Occs);
extern void s2_block_vfci(struct stringwr **alplist,
   struct stringwr **betlist,
   double **C, double **S, double *oei, double *tei, double *F,
   int nlists, int nas, int nbs, int Ia_list, int Ja_list,
   int Ja_list_nas);
extern void s2_block_vras(struct stringwr **alplist,
   struct stringwr **betlist,
   double **C, double **S, double *oei, double *tei, double *F,
   int nlists, int nas, int nbs, int sac, int cac, int cnas);
extern void s2_block_vras_rotf(int *Cnt[2], int **Ij[2], int **Oij[2],
   int **Ridx[2], signed char **Sgn[2], unsigned char **Toccs,
   double **C, double **S,
   double *oei, double *tei, double *F, int nlists, int nas, int nbs,
   int Ia_list, int Ja_list, int Ja_list_nbs, struct olsen_graph *AlphaG,
   struct olsen_graph *BetaG, struct calcinfo *CIinfo, unsigned char ***Occs);
extern void s3_block_vdiag(struct stringwr *alplist,
   struct stringwr *betlist,
   double **C, double **S, double *tei, int nas, int nbs, int cnas,
   int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
   double **Cprime, double *F, double *V, double *Sgn, int *L, int *R,
   int norbs, int *orbsym);
extern void s3_block_v(struct stringwr *alplist,struct stringwr *betlist,
   double **C, double **S, double *tei, int nas, int nbs, int cnas,
   int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
   double **Cprime, double *F, double *V, double *Sgn, int *L, int *R,
   int norbs, int *orbsym);
extern void s3_block_vrotf(int *Cnt[2], int **Ij[2], int **Ridx[2],
   signed char **Sn[2], double **C, double **S,
   double *tei, int nas, int nbs, int cnas,
   int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
   double **Cprime, double *F, double *V, double *Sgn, int *L, int *R,
   int norbs, int *orbsym);
extern void s3_block_vdiag_rotf(int *Cnt[2], int **Ij[2], int **Ridx[2],
   signed char **Sn[2], double **C, double **S,
   double *tei, int nas, int nbs, int cnas,
   int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
   double **Cprime, double *F, double *V, double *Sgn, int *L, int *R,
   int norbs, int *orbsym);



/*
** sigma_init()
**
** This function initializes all the globals associated with calculating
** the sigma vector.
**
*/
void CIWavefunction::sigma_init(CIvect& C, CIvect &S)
{
   int i,j;
   int maxcols=0, maxrows=0;
   int nsingles, max_dim=0;
   unsigned long int bufsz=0;

   SigmaData_->transp_tmp = NULL;
   SigmaData_->cprime = NULL;
   SigmaData_->sprime = NULL;

   if (CalcInfo_->sigma_initialized) {
       //outfile->Printf("(sigma_init): sigma_initialized already set to 1\n");
       return;
   }

   for (i = 0; i < C.num_blocks_; i++) {
       if (C.Ib_size_[i] > max_dim) max_dim = C.Ib_size_[i];
       if (C.Ia_size_[i] > max_dim) max_dim = C.Ia_size_[i];
   }
   SigmaData_->max_dim = max_dim;
   SigmaData_->F = init_array(max_dim);

   SigmaData_->Sgn = init_array(max_dim);
   SigmaData_->V = init_array(max_dim);
   SigmaData_->L = init_int_array(max_dim);
   SigmaData_->R = init_int_array(max_dim);

   if (Parameters_->repl_otf) {
       max_dim += AlphaG_->num_el_expl;
       nsingles = AlphaG_->num_el_expl * AlphaG_->num_orb;
       for (i = 0; i < 2; i++) {
           SigmaData_->Jcnt[i] = init_int_array(max_dim);
           SigmaData_->Jij[i] = init_int_matrix(max_dim, nsingles);
           SigmaData_->Joij[i] = init_int_matrix(max_dim, nsingles);
           SigmaData_->Jridx[i] = init_int_matrix(max_dim, nsingles);
           SigmaData_->Jsgn[i] =
               (signed char **)malloc(max_dim * sizeof(signed char *));
           for (j = 0; j < max_dim; j++) {
               SigmaData_->Jsgn[i][j] =
                   (signed char *)malloc(nsingles * sizeof(signed char));
           }
       }

       SigmaData_->Toccs =
           (unsigned char **)malloc(sizeof(unsigned char *) * nsingles);

       /* test out the on-the-fly replacement routines */
       /*
       b2brepl_test(Occs_,SigmaData_->Jcnt[0],SigmaData_->Jij[0],
                    SigmaData_->Joij[0],SigmaData_->Jridx[0],SigmaData_->Jsgn[0],AlphaG);
       */
   }

   /* figure out which C blocks contribute to s */
   s1_contrib_ = init_int_matrix(S.num_blocks_, C.num_blocks_);
   s2_contrib_ = init_int_matrix(S.num_blocks_, C.num_blocks_);
   s3_contrib_ = init_int_matrix(S.num_blocks_, C.num_blocks_);
   if (Parameters_->repl_otf)
       sigma_get_contrib_rotf(C, S, s1_contrib_, s2_contrib_, s3_contrib_,
                              SigmaData_->Jcnt, SigmaData_->Jij,
                              SigmaData_->Joij, SigmaData_->Jridx,
                              SigmaData_->Jsgn, SigmaData_->Toccs);
   else
       sigma_get_contrib(alplist_, betlist_, C, S, s1_contrib_, s2_contrib_,
                         s3_contrib_);

   if ((C.icore_ == 2 && C.Ms0_ && CalcInfo_->ref_sym != 0) ||
       (C.icore_ == 0 && C.Ms0_)) {
       for (i = 0, maxrows = 0, maxcols = 0; i < C.num_blocks_; i++) {
           if (C.Ia_size_[i] > maxrows) maxrows = C.Ia_size_[i];
           if (C.Ib_size_[i] > maxcols) maxcols = C.Ib_size_[i];
       }
       if (maxcols > maxrows) maxrows = maxcols;
       SigmaData_->transp_tmp = (double **)malloc(maxrows * sizeof(double *));
       if (SigmaData_->transp_tmp == NULL) {
          outfile->Printf(
               "(sigma_init): Trouble with malloc'ing "
               "SigmaData_->transp_tmp\n");
       }
       bufsz = C.get_max_blk_size();
       SigmaData_->transp_tmp[0] = init_array(bufsz);
       if (SigmaData_->transp_tmp[0] == NULL) {
          outfile->Printf(
               "(sigma_init): Trouble with malloc'ing "
               "SigmaData_->transp_tmp[0]\n");
       }
   }

   /* make room for SigmaData_->cprime and SigmaData_->sprime if necessary */
   for (i = 0, maxrows = 0; i < C.num_blocks_; i++) {
       if (C.Ia_size_[i] > maxrows) maxrows = C.Ia_size_[i];
       if (C.Ib_size_[i] > maxcols) maxcols = C.Ib_size_[i];
   }
   if ((C.icore_ == 2 && C.Ms0_ && CalcInfo_->ref_sym != 0) ||
       (C.icore_ == 0 && C.Ms0_)) {
       if (maxcols > maxrows) maxrows = maxcols;
   }
   bufsz = C.get_max_blk_size();

   SigmaData_->cprime = (double **)malloc(maxrows * sizeof(double *));
   if (SigmaData_->cprime == NULL) {
      outfile->Printf("(sigma_init): Trouble with malloc'ing SigmaData_->cprime\n");
   }
   if (C.icore_ == 0 && C.Ms0_ && SigmaData_->transp_tmp != NULL &&
       SigmaData_->transp_tmp[0] != NULL)
       SigmaData_->cprime[0] = SigmaData_->transp_tmp[0];
   else
       SigmaData_->cprime[0] = init_array(bufsz);

   if (SigmaData_->cprime[0] == NULL) {
      outfile->Printf("(sigma_init): Trouble with malloc'ing SigmaData_->cprime[0]\n");
   }

   if (Parameters_->bendazzoli) {
       SigmaData_->sprime = (double **)malloc(maxrows * sizeof(double *));
       if (SigmaData_->sprime == NULL) {
          outfile->Printf("(sigma_init): Trouble with malloc'ing SigmaData_->sprime\n");
       }
       SigmaData_->sprime[0] = init_array(bufsz);
       if (SigmaData_->sprime[0] == NULL) {
          outfile->Printf(
               "(sigma_init): Trouble with malloc'ing SigmaData_->sprime[0]\n");
       }
   }

   CalcInfo_->sigma_initialized = 1;
}

void CIWavefunction::sigma_free()
{
   free(SigmaData_->F);
   free(SigmaData_->Sgn);
   free(SigmaData_->V);
   free(SigmaData_->L);
   free(SigmaData_->R);
   if (Parameters_->repl_otf) {
      for (int i=0; i<2; i++) {
         free(SigmaData_->Jcnt[i]);
         free_int_matrix(SigmaData_->Jij[i]);
         free_int_matrix(SigmaData_->Joij[i]);
         free_int_matrix(SigmaData_->Jridx[i]);
         for (int j=0; j<SigmaData_->max_dim; j++) {
            free(SigmaData_->Jsgn[i][j]);
         }
         free(SigmaData_->Jsgn[i]);
      }
   }
   CalcInfo_->sigma_initialized = false;
// DGAS: Not sure how to free these yet
//      SigmaData_->Toccs = (unsigned char **) malloc (sizeof(unsigned char *) * nsingles);
//unsigned char **SigmaData_->Toccs;
//double **SigmaData_->transp_tmp, **SigmaData_->cprime, **SigmaData_->sprime;
}

/*
** sigma()
**
** Routine to get the sigma vector using the CI vector class
**
** Changed into a master function which calls the appropriate subfunction
**
*/
void CIWavefunction::sigma(CIvect &C, CIvect &S, double *oei, double *tei, int ivec) {
    if (!CalcInfo_->sigma_initialized) sigma_init(C, S);
    int fci = Parameters_->fci;

    switch (C.icore_) {
        case 0:
            sigma_a(alplist_, betlist_, C, S, oei, tei, fci, ivec);
            break;
        case 1:
            sigma_b(alplist_, betlist_, C, S, oei, tei, fci, ivec);
            break;
        case 2:
            sigma_c(alplist_, betlist_, C, S, oei, tei, fci, ivec);
            break;
        default:
            outfile->Printf("(sigma): Error, invalid icore option\n");
            break;
    }
}
void CIWavefunction::sigma(SharedCIVector C, SharedCIVector S, int cvec, int svec) {
    C->cur_vect_ = cvec;
    double *oei;
    if (Parameters_->fci)
        oei = CalcInfo_->tf_onel_ints->pointer();
    else
        oei = CalcInfo_->gmat->pointer();
    sigma(*(C.get()), *(S.get()), oei, CalcInfo_->twoel_ints->pointer(), svec);
}
void CIWavefunction::sigma(SharedCIVector C, SharedCIVector S, int cvec,
                           int svec, SharedVector oei, SharedVector tei) {
    if ((oei->nirrep() != 1) || (tei->nirrep() != 1)){
      throw PSIEXCEPTION("CIWavefunction::sigma: Electron integrals cannot have irreps");
    }
    C->cur_vect_ = cvec;
    sigma(*(C.get()), *(S.get()), oei->pointer(), tei->pointer(), svec);
}

/*
** sigma_a(): This function computes the sigma vector for a given C vector
**    using the CIvector class.   Is somewhat intelligent about constructing
**    the sigma vector blockwise; will attempt to reduce I/O for out-of-core
**    cases.  This version is for icore==0.
**
** Parameters:
**    alplist  = list of alpha strings with replacements
**    betlist  = same for beta strings
**    C        = current ci vector
**    S        = sigma vector to be computed
**    oei      = one-electron integrals
**    tei      = two-electron integrals
**    fci      = full-ci flag (helps determine which sigma1 routine called)
**    ivec     = sigma vector number (for write call)
**
** Notes: assumes M_s = 0 for now
*/
void CIWavefunction::sigma_a(struct stringwr **alplist, struct stringwr **betlist,
      CIvect& C, CIvect& S, double *oei, double *tei, int fci, int ivec)
{

   int buf, cbuf;
   int sblock, cblock, cblock2;  /* id of sigma and C blocks */
   int sac, sbc, nas, nbs;
   int cac, cbc, cnas, cnbs;
   int do_cblock, do_cblock2;
   int cairr, cbirr, sbirr;
   int did_sblock = 0;
   int phase;

   if (!Parameters_->Ms0) phase = 1;
   else phase = ((int) Parameters_->S % 2) ? -1 : 1;

   /* this does a sigma subblock at a time: icore==0 */
   for (buf=0; buf<S.buf_per_vect_; buf++) {
      S.zero();
      did_sblock = 0;
      sblock = S.buf2blk_[buf];
      sac = S.Ia_code_[sblock];
      sbc = S.Ib_code_[sblock];
      nas = S.Ia_size_[sblock];
      nbs = S.Ib_size_[sblock];
      sbirr = sbc / BetaG_->subgr_per_irrep;
      if (SigmaData_->sprime != NULL) set_row_ptrs(nas, nbs, SigmaData_->sprime);

      for (cbuf=0; cbuf<C.buf_per_vect_; cbuf++) {
         do_cblock=0; do_cblock2=0;
         cblock=C.buf2blk_[cbuf];
         cblock2 = -1;
         cac = C.Ia_code_[cblock];
         cbc = C.Ib_code_[cblock];
         cbirr = cbc / BetaG_->subgr_per_irrep;
         cairr = cac / AlphaG_->subgr_per_irrep;
         if (C.Ms0_) cblock2 = C.decode_[cbc][cac];
         cnas = C.Ia_size_[cblock];
         cnbs = C.Ib_size_[cblock];
         if (s1_contrib_[sblock][cblock] || s2_contrib_[sblock][cblock] ||
             s3_contrib_[sblock][cblock]) do_cblock = 1;
         if (C.buf_offdiag_[cbuf] && (s1_contrib_[sblock][cblock2] ||
             s2_contrib_[sblock][cblock2] || s3_contrib_[sblock][cblock2]))
            do_cblock2 = 1;
         if (C.check_zero_block(cblock)) do_cblock = 0;
         if (cblock2 >= 0 && C.check_zero_block(cblock2)) do_cblock2 = 0;
         if (!do_cblock && !do_cblock2) continue;

         C.read(C.cur_vect_, cbuf);

         if (do_cblock) {
            if (SigmaData_->cprime != NULL) set_row_ptrs(cnas, cnbs, SigmaData_->cprime);
            sigma_block(alplist, betlist, C.blocks_[cblock], S.blocks_[sblock],
               oei, tei, fci, cblock, sblock, nas, nbs, sac, sbc, cac, cbc,
               cnas, cnbs, C.num_alpcodes_, C.num_betcodes_, sbirr, cbirr,
               S.Ms0_);
            did_sblock = 1;
            }

         /* what's with this bcopy stuff?  what's going on?  -DS 6/11/96 */
         /* I think I should copy to cblock2 not cblock */
         if (do_cblock2) {
            C.transp_block(cblock, SigmaData_->transp_tmp);
//          bcopy((char *) SigmaData_->transp_tmp[0], (char *) C.blocks_[cblock][0],
//            cnas * cnbs * sizeof(double));
//          bcopy is non-ANSI.  memcpy reverses the arguments.
            memcpy((void *) C.blocks_[cblock][0], (void *) SigmaData_->transp_tmp[0],
              cnas * cnbs * sizeof(double));
            /* set_row_ptrs(cnbs, cnas, C.blocks_[cblock]); */
            if (SigmaData_->cprime != NULL) set_row_ptrs(cnbs, cnas, SigmaData_->cprime);
            sigma_block(alplist, betlist, C.blocks_[cblock2], S.blocks_[sblock],
               oei, tei, fci, cblock2, sblock, nas, nbs, sac, sbc,
               cbc, cac, cnbs, cnas, C.num_alpcodes_, C.num_betcodes_, sbirr,
               cairr, S.Ms0_);
            did_sblock = 1;
            }

         } /* end loop over c buffers */

      if (did_sblock) {
         S.set_zero_block(sblock, 0);
         if (S.Ms0_) S.set_zero_block(S.decode_[sbc][sac], 0);
         }

      if (S.Ms0_ && (sac==sbc))
         transp_sigma(S.blocks_[sblock], nas, nbs, phase);

      H0block_gather(S.blocks_[sblock], sac, sbc, 1, Parameters_->Ms0, phase);

      if (S.Ms0_) {
         if ((int) Parameters_->S % 2) S.symmetrize(-1.0, sblock);
         else S.symmetrize(1.0, sblock);
         }
      S.write(ivec, buf);

      } /* end loop over sigma buffers */

}



/*
** sigma_b(): This function computes the sigma vector for a given C vector
**    using the CIvector class.   Is somewhat intelligent about constructing
**    the sigma vector blockwise; will attempt to reduce I/O for out-of-core
**    cases.  This version is for icore=1 (whole vector in-core)
**
** Parameters:
**    alplist  = list of alpha strings with replacements
**    betlist  = same for beta strings
**    C        = current ci vector
**    S        = sigma vector to be computed
**    oei      = one-electron integrals
**    tei      = two-electron integrals
**    fci      = full-ci flag (helps determine which sigma1 routine called)
**    ivec     = sigma vector number (for write call)
**
** Notes: I think I removed the M_s = 0 assumption from this one
*/
void CIWavefunction::sigma_b(struct stringwr **alplist, struct stringwr **betlist,
      CIvect& C, CIvect& S, double *oei, double *tei, int fci, int ivec)
{

   int sblock, cblock;  /* id of sigma and C blocks */
   int sac, sbc, nas, nbs;
   int cac, cbc, cnas, cnbs;
   int sbirr, cbirr;
   int did_sblock = 0;
   int phase;

   if (!Parameters_->Ms0) phase = 1;
   else phase = ((int) Parameters_->S % 2) ? -1 : 1;

   S.zero();
   C.read(C.cur_vect_, 0);

   /* loop over unique sigma subblocks */
   for (sblock=0; sblock<S.num_blocks_; sblock++) {
      //if (Parameters_->cc && !cc_reqd_sblocks[sblock]) continue;
      did_sblock = 0;
      sac = S.Ia_code_[sblock];
      sbc = S.Ib_code_[sblock];
      nas = S.Ia_size_[sblock];
      nbs = S.Ib_size_[sblock];
      if (nas==0 || nbs==0) continue;
      if (S.Ms0_ && sbc > sac) continue;
      sbirr = sbc / BetaG_->subgr_per_irrep;
      if (SigmaData_->sprime != NULL) set_row_ptrs(nas, nbs, SigmaData_->sprime);

      for (cblock=0; cblock<C.num_blocks_; cblock++) {
         if (C.check_zero_block(cblock)) continue;
         cac = C.Ia_code_[cblock];
         cbc = C.Ib_code_[cblock];
         cnas = C.Ia_size_[cblock];
         cnbs = C.Ib_size_[cblock];
         cbirr = cbc / BetaG_->subgr_per_irrep;
         if (s1_contrib_[sblock][cblock] || s2_contrib_[sblock][cblock] ||
             s3_contrib_[sblock][cblock]) {
            if (SigmaData_->cprime != NULL) set_row_ptrs(cnas, cnbs, SigmaData_->cprime);
            sigma_block(alplist, betlist, C.blocks_[cblock], S.blocks_[sblock],
               oei, tei, fci, cblock, sblock, nas, nbs, sac, sbc,
               cac, cbc, cnas, cnbs, C.num_alpcodes_, C.num_betcodes_, sbirr,
               cbirr, S.Ms0_);
            did_sblock = 1;
            }
         } /* end loop over c blocks */

      if (did_sblock) S.set_zero_block(sblock, 0);

      if (S.Ms0_ && (sac==sbc))
         transp_sigma(S.blocks_[sblock], nas, nbs, phase);
      H0block_gather(S.blocks_[sblock], sac, sbc, 1, Parameters_->Ms0,
         phase);
      } /* end loop over sigma blocks */

      if (S.Ms0_) {
         if ((int) Parameters_->S % 2) S.symmetrize(-1.0, 0);
         else S.symmetrize(1.0, 0);
         }

   S.write(ivec, 0);

}


/*
** sigma_c(): This function computes the sigma vector for a given C vector
**    using the CIvector class.   Is somewhat intelligent about constructing
**    the sigma vector blockwise; will attempt to reduce I/O for out-of-core
**    cases.  This version is for icore=2 (irrep at a time)
**
** Parameters:
**    alplist  = list of alpha strings with replacements
**    betlist  = same for beta strings
**    C        = current ci vector
**    S        = sigma vector to be computed
**    oei      = one-electron integrals
**    tei      = two-electron integrals
**    fci      = full-ci flag (helps determine which sigma1 routine called)
**    ivec     = sigma vector number (for write call)
**
** Notes: tried to remove Ms=0 assumption
*/
void CIWavefunction::sigma_c(struct stringwr **alplist, struct stringwr **betlist,
      CIvect& C, CIvect& S, double *oei, double *tei, int fci, int ivec)
{

   int buf, cbuf;
   int sblock, cblock, cblock2;  /* id of sigma and C blocks */
   int sairr;                    /* irrep of alpha string for sigma block */
   int cairr;                    /* irrep of alpha string for C block */
   int sbirr, cbirr;
   int sac, sbc, nas, nbs;
   int cac, cbc, cnas, cnbs;
   int did_sblock = 0;
   int phase;

   if (!Parameters_->Ms0) phase = 1;
   else phase = ((int) Parameters_->S % 2) ? -1 : 1;


   for (buf=0; buf<S.buf_per_vect_; buf++) {
      sairr = S.buf2blk_[buf];
      sbirr = sairr ^ CalcInfo_->ref_sym;
      S.zero();
      for (cbuf=0; cbuf<C.buf_per_vect_; cbuf++) {
         C.read(C.cur_vect_, cbuf); /* go ahead and assume it will contrib */
         cairr = C.buf2blk_[cbuf];
         cbirr = cairr ^ CalcInfo_->ref_sym;

         for (sblock=S.first_ablk_[sairr];sblock<=S.last_ablk_[sairr];sblock++){
            sac = S.Ia_code_[sblock];
            sbc = S.Ib_code_[sblock];
            nas = S.Ia_size_[sblock];
            nbs = S.Ib_size_[sblock];
            did_sblock = 0;

            if (S.Ms0_ && (sac < sbc)) continue;
            if (SigmaData_->sprime != NULL) set_row_ptrs(nas, nbs, SigmaData_->sprime);

            for (cblock=C.first_ablk_[cairr]; cblock <= C.last_ablk_[cairr];
                  cblock++) {

               cac = C.Ia_code_[cblock];
               cbc = C.Ib_code_[cblock];
               cnas = C.Ia_size_[cblock];
               cnbs = C.Ib_size_[cblock];

               if ((s1_contrib_[sblock][cblock] || s2_contrib_[sblock][cblock] ||
                    s3_contrib_[sblock][cblock]) &&
                    !C.check_zero_block(cblock)) {
      if (SigmaData_->cprime != NULL) set_row_ptrs(cnas, cnbs, SigmaData_->cprime);
                  sigma_block(alplist, betlist, C.blocks_[cblock],
                     S.blocks_[sblock], oei, tei, fci, cblock,
                     sblock, nas, nbs, sac, sbc, cac, cbc, cnas, cnbs,
                     C.num_alpcodes_, C.num_betcodes_, sbirr, cbirr, S.Ms0_);
                  did_sblock = 1;
                  }

               if (C.buf_offdiag_[cbuf]) {
                  cblock2 = C.decode_[cbc][cac];
                  if ((s1_contrib_[sblock][cblock2] ||
                       s2_contrib_[sblock][cblock2] ||
                       s3_contrib_[sblock][cblock2]) &&
                      !C.check_zero_block(cblock2)) {
                     C.transp_block(cblock, SigmaData_->transp_tmp);
         if (SigmaData_->cprime != NULL) set_row_ptrs(cnbs, cnas, SigmaData_->cprime);
                     sigma_block(alplist, betlist, SigmaData_->transp_tmp,S.blocks_[sblock],
                        oei, tei, fci, cblock2, sblock, nas, nbs, sac, sbc,
                        cbc, cac, cnbs, cnas, C.num_alpcodes_, C.num_betcodes_,
                        sbirr, cairr, S.Ms0_);
                     did_sblock = 1;
                     }
                  }
               } /* end loop over C blocks in this irrep */

            if (did_sblock) S.set_zero_block(sblock, 0);
            } /* end loop over sblock */

         } /* end loop over cbuf */

      /* transpose the diagonal sigma subblocks in this irrep */
      for (sblock=S.first_ablk_[sairr];sblock<=S.last_ablk_[sairr];sblock++){
         sac = S.Ia_code_[sblock];
         sbc = S.Ib_code_[sblock];
         nas = S.Ia_size_[sblock];
         nbs = S.Ib_size_[sblock];
         if (S.Ms0_ && (sac==sbc)) transp_sigma(S.blocks_[sblock], nas, nbs,
            phase);

         /* also gather the contributions from sigma to the H0block */
         if (!S.Ms0_ || sac >= sbc) {
            H0block_gather(S.blocks_[sblock], sac, sbc, 1, Parameters_->Ms0,
               phase);
            }
         }

      if (S.Ms0_) {
         if ((int) Parameters_->S % 2) S.symmetrize(-1.0, sairr);
         else S.symmetrize(1.0, sairr);
         }
     S.write(ivec, buf);

     } /* end loop over sigma irrep */
}


/*
** sigma_block()
**
** Calculate the contribution to sigma block sblock from C block cblock
**
*/
void CIWavefunction::sigma_block(struct stringwr **alplist, struct stringwr **betlist,
      double **cmat, double **smat, double *oei, double *tei, int fci,
      int cblock, int sblock, int nas, int nbs, int sac, int sbc,
      int cac, int cbc, int cnas, int cnbs, int cnac, int cnbc,
      int sbirr, int cbirr, int Ms0)
{

   /* SIGMA2 CONTRIBUTION */
  if (s2_contrib_[sblock][cblock]) {

    timer_on("CIWave: s2");

      if (fci) {
          s2_block_vfci(alplist, betlist, cmat, smat, oei, tei, SigmaData_->F, cnac,
                            nas, nbs, sac, cac, cnas);
        }
      else {
          if (Parameters_->repl_otf) {
              s2_block_vras_rotf(SigmaData_->Jcnt, SigmaData_->Jij, SigmaData_->Joij,
                                 SigmaData_->Jridx, SigmaData_->Jsgn,
                                 SigmaData_->Toccs, cmat, smat, oei, tei, SigmaData_->F, cnac,
                                 nas, nbs, sac, cac, cnas, AlphaG_, BetaG_, CalcInfo_, Occs_);
            }
          else {
              s2_block_vras(alplist, betlist, cmat, smat,
                            oei, tei, SigmaData_->F, cnac, nas, nbs, sac, cac, cnas);
            }
        }
    timer_off("CIWave: s2");

    } /* end sigma2 */


   if (print_ > 3) {
     outfile->Printf( "s2: Contribution to sblock=%d from cblock=%d\n",
        sblock, cblock);
     print_mat(smat, nas, nbs, "outfile");
   }

   /* SIGMA1 CONTRIBUTION */
   if (!Ms0 || (sac != sbc)) {
    timer_on("CIWave: s1");

      if (s1_contrib_[sblock][cblock]) {
          if (fci) {
             s1_block_vfci(alplist, betlist, cmat, smat, oei, tei, SigmaData_->F, cnbc,
                                nas, nbs, sbc, cbc, cnbs);
            }
         else {
            if (Parameters_->repl_otf) {
               s1_block_vras_rotf(SigmaData_->Jcnt, SigmaData_->Jij, SigmaData_->Joij,
                  SigmaData_->Jridx, SigmaData_->Jsgn,
                  SigmaData_->Toccs, cmat, smat, oei, tei, SigmaData_->F, cnbc, nas, nbs,
                  sbc, cbc, cnbs, BetaG_, CalcInfo_, Occs_);
               }
            else {
               s1_block_vras(alplist, betlist, cmat, smat, oei, tei, SigmaData_->F, cnbc,
                  nas, nbs, sbc, cbc, cnbs);
               }
            }
         }

      timer_off("CIWave: s1");
   } /* end sigma1 */

   if (print_ > 3) {
     outfile->Printf( "s1: Contribution to sblock=%d from cblock=%d\n",
        sblock, cblock);
     print_mat(smat, nas, nbs, "outfile");
   }

   /* SIGMA3 CONTRIBUTION */
   if (s3_contrib_[sblock][cblock]) {
      timer_on("CIWave: s3");

      /* zero_mat(smat, nas, nbs); */

      if (!Ms0 || (sac != sbc)) {
         if (Parameters_->repl_otf) {
            b2brepl(Occs_[sac], SigmaData_->Jcnt[0], SigmaData_->Jij[0],
               SigmaData_->Joij[0], SigmaData_->Jridx[0],
               SigmaData_->Jsgn[0], AlphaG_, sac, cac, nas, CalcInfo_);
            b2brepl(Occs_[sbc], SigmaData_->Jcnt[1], SigmaData_->Jij[1],
                    SigmaData_->Joij[1], SigmaData_->Jridx[1],
                    SigmaData_->Jsgn[1], BetaG_, sbc, cbc, nbs, CalcInfo_);
            s3_block_vrotf(SigmaData_->Jcnt, SigmaData_->Jij, SigmaData_->Jridx,
                           SigmaData_->Jsgn, cmat, smat, tei, nas, nbs,
                           cnas, sbc, cac, cbc, sbirr, cbirr, SigmaData_->cprime,
                           SigmaData_->F, SigmaData_->V, SigmaData_->Sgn, SigmaData_->L,
                           SigmaData_->R, CalcInfo_->num_ci_orbs,
                           CalcInfo_->orbsym + CalcInfo_->num_drc_orbs);
            }
         else {
            s3_block_v(alplist[sac], betlist[sbc], cmat, smat, tei,
               nas, nbs, cnas, sbc, cac, cbc, sbirr, cbirr,
               SigmaData_->cprime, SigmaData_->F, SigmaData_->V,
               SigmaData_->Sgn, SigmaData_->L, SigmaData_->R,
               CalcInfo_->num_ci_orbs, CalcInfo_->orbsym + CalcInfo_->num_drc_orbs);
            }
         }

      else if (Parameters_->bendazzoli) {
         s3_block_bz(sac, sbc, cac, cbc, nas, nbs, cnas, tei, cmat, smat,
            SigmaData_->cprime, SigmaData_->sprime, CalcInfo_, OV_);
         }

      else {
         if (Parameters_->repl_otf) {
            b2brepl(Occs_[sac], SigmaData_->Jcnt[0], SigmaData_->Jij[0],
                    SigmaData_->Joij[0], SigmaData_->Jridx[0],
                    SigmaData_->Jsgn[0], AlphaG_, sac, cac, nas, CalcInfo_);
            b2brepl(Occs_[sbc], SigmaData_->Jcnt[1], SigmaData_->Jij[1],
                    SigmaData_->Joij[1], SigmaData_->Jridx[1],
                    SigmaData_->Jsgn[1], BetaG_, sbc, cbc, nbs, CalcInfo_);
            s3_block_vdiag_rotf(SigmaData_->Jcnt, SigmaData_->Jij, SigmaData_->Jridx,
                                SigmaData_->Jsgn, cmat, smat, tei, nas, nbs, cnas, sbc,
                                cac, cbc, sbirr, cbirr, SigmaData_->cprime, SigmaData_->F,
                                SigmaData_->V, SigmaData_->Sgn, SigmaData_->L,
                                SigmaData_->R, CalcInfo_->num_ci_orbs,
                                CalcInfo_->orbsym + CalcInfo_->num_drc_orbs);
            }
         else {
            s3_block_vdiag(alplist[sac], betlist[sbc], cmat, smat, tei, nas, nbs,
                           cnas, sbc, cac, cbc, sbirr, cbirr, SigmaData_->cprime,
                           SigmaData_->F, SigmaData_->V, SigmaData_->Sgn,
                           SigmaData_->L, SigmaData_->R, CalcInfo_->num_ci_orbs,
                           CalcInfo_->orbsym + CalcInfo_->num_drc_orbs);
            }
         }

      if (print_ > 3) {
        outfile->Printf( "s3: Contribution to sblock=%d from cblock=%d\n",
           sblock, cblock);
        print_mat(smat, nas, nbs, "outfile");
      }

      timer_off("CIWave: s3");

      } /* end sigma3 */
}

void CIWavefunction::sigma_get_contrib(struct stringwr **alplist, struct stringwr **betlist,
      CIvect &C, CIvect &S, int **s1_contrib, int **s2_contrib,
      int **s3_contrib)
{

   int sblock,cblock;
   int sac, sbc, cac, cbc;
   int nas, nbs;
   struct stringwr *Ib, *Ia, *Kb, *Ka;
   unsigned int Ibidx, Iaidx, Kbidx, Kaidx, Ib_ex, Ia_ex;
   unsigned int Ibcnt, Iacnt, *Ibridx, *Iaridx;
   int Kb_list, Ka_list;
   int found,i,j;

   for (sblock=0; sblock<S.num_blocks_; sblock++) {
      sac = S.Ia_code_[sblock];
      sbc = S.Ib_code_[sblock];
      nas = S.Ia_size_[sblock];
      nbs = S.Ib_size_[sblock];
      for (cblock=0; cblock<C.num_blocks_; cblock++) {
         cac = C.Ia_code_[cblock];
         cbc = C.Ib_code_[cblock];


         /* does this c block contribute to sigma1? */
         if (sac == cac) {
            for (Ib=betlist[sbc], Ibidx=0, found=0; Ibidx < nbs && !found;
               Ibidx++, Ib++) {
               /* loop over excitations E^b_{kl} from |B(I_b)> */
               for (Kb_list=0; Kb_list < S.num_betcodes_ && !found; Kb_list++) {
                  Ibcnt = Ib->cnt[Kb_list];
                  Ibridx = Ib->ridx[Kb_list];
                  for (Ib_ex=0; Ib_ex < Ibcnt; Ib_ex++) {
                     Kbidx = *Ibridx++;
                     Kb = betlist[Kb_list] + Kbidx;
                     if (Kb->cnt[cbc]) { found=1;  break; }
                     }
                  }
               }
            if (found) s1_contrib[sblock][cblock] = 1;
            }

         /* does this c block contribute to sigma2? */
         if (sbc == cbc) {
         for (Ia=alplist[sac], Iaidx=0, found=0; Iaidx < nas && !found;
               Iaidx++, Ia++) {
               /* loop over excitations E^a_{kl} from |A(I_a)> */
               for (Ka_list=0; Ka_list < S.num_alpcodes_ && !found; Ka_list++) {
                  Iacnt = Ia->cnt[Ka_list];
                  Iaridx = Ia->ridx[Ka_list];
                  for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
                     Kaidx = *Iaridx++;
                     Ka = alplist[Ka_list] + Kaidx;
                     if (Ka->cnt[cac]) { found=1;  break; }
                     }
                  }
               }
            if (found) s2_contrib[sblock][cblock] = 1;
            }

         /* does this c block contribute to sigma3? */
         for (Iaidx=0,found=0; Iaidx<S.Ia_size_[sblock]; Iaidx++) {
            if (alplist[sac][Iaidx].cnt[cac]) found=1;
            }
         if (found) { /* see if beta is ok */
            found=0;
            for (Ibidx=0; Ibidx<S.Ib_size_[sblock]; Ibidx++) {
               if (betlist[sbc][Ibidx].cnt[cbc]) found=1;
               }
            if (found)
               s3_contrib[sblock][cblock] = 1;
            }

         } /* end loop over c blocks */
      } /* end loop over sigma blocks */

   if (print_ > 4) {
    outfile->Printf("\nSigma 1:\n");
     for (i=0; i<S.num_blocks_; i++) {
       outfile->Printf( "Contributions to sigma block %d\n", i);
       for (j=0; j<C.num_blocks_; j++) {
         if (s1_contrib[i][j]) outfile->Printf( "%3d ", j);
       }
       outfile->Printf( "\n");
     }

    outfile->Printf("\n\nSigma 2:\n");
     for (i=0; i<S.num_blocks_; i++) {
       outfile->Printf( "Contributions to sigma block %d\n", i);
       for (j=0; j<C.num_blocks_; j++) {
         if (s2_contrib[i][j]) outfile->Printf( "%3d ", j);
       }
       outfile->Printf( "\n");
     }

    outfile->Printf("\n\nSigma 3:\n");
     for (i=0; i<S.num_blocks_; i++) {
       outfile->Printf( "Contributions to sigma block %d\n", i);
       for (j=0; j<C.num_blocks_; j++) {
         if (s3_contrib[i][j]) outfile->Printf( "%3d ", j);
       }
       outfile->Printf( "\n");
     }
   }

}

void CIWavefunction::sigma_get_contrib_rotf(CIvect &C, CIvect &S,
      int **s1_contrib, int **s2_contrib, int **s3_contrib,
      int *Cnt[2], int **Ij[2], int **Oij[2], int **Ridx[2],
      signed char **Sgn[2], unsigned char **Toccs)
{

   int sblock,cblock;
   int sac, sbc, cac, cbc;
   int nas, nbs;
   int Ibidx, Iaidx, Ib_ex, Ia_ex;
   int Ibcnt, Iacnt;
   int Kb_list, Ka_list;
   int found,i,j;

   for (sblock=0; sblock<S.num_blocks_; sblock++) {
      sac = S.Ia_code_[sblock];
      sbc = S.Ib_code_[sblock];
      nas = S.Ia_size_[sblock];
      nbs = S.Ib_size_[sblock];
      for (cblock=0; cblock<C.num_blocks_; cblock++) {
         cac = C.Ia_code_[cblock];
         cbc = C.Ib_code_[cblock];


         /* does this c block contribute to sigma1? */
         if (sac == cac) {
            found = 0;
            for (Kb_list=0; Kb_list < S.num_betcodes_ && !found; Kb_list++) {
               b2brepl(Occs_[sbc], Cnt[0], Ij[0], Oij[0], Ridx[0],
                  Sgn[0], BetaG_, sbc, Kb_list, nbs, CalcInfo_);
               for (Ibidx=0; Ibidx < nbs && !found; Ibidx++) {
                  Ibcnt = Cnt[0][Ibidx];
                  if (Ibcnt) {
                     for (i=0; i<Ibcnt; i++) {
                        j = Ridx[0][Ibidx][i];
                        Toccs[i] = Occs_[Kb_list][j];
                        }
                     b2brepl(Toccs, Cnt[1], Ij[1], Oij[1], Ridx[1], Sgn[1],
                        BetaG_, Kb_list, cbc, Ibcnt, CalcInfo_);
                     for (Ib_ex=0; Ib_ex < Ibcnt; Ib_ex++) {
                        if (Cnt[1][Ib_ex]) { found=1; break; }
                        }
                     }
                  }
               }
            if (found) s1_contrib[sblock][cblock] = 1;
            }

         /* does this c block contribute to sigma2? */
         if (sbc == cbc) {
            found = 0;
            for (Ka_list=0; Ka_list < S.num_alpcodes_ && !found; Ka_list++) {
               b2brepl(Occs_[sac], Cnt[0], Ij[0], Oij[0], Ridx[0],
                  Sgn[0], AlphaG_, sac, Ka_list, nas, CalcInfo_);
               for (Iaidx=0; Iaidx < nas && !found; Iaidx++) {
                  Iacnt = Cnt[0][Iaidx];
                  if (Iacnt) {
                     for (i=0; i<Iacnt; i++) {
                        j = Ridx[0][Iaidx][i];
                        Toccs[i] = Occs_[Ka_list][j];
                        }
                     b2brepl(Toccs, Cnt[1], Ij[1], Oij[1], Ridx[1], Sgn[1],
                        AlphaG_, Ka_list, cac, Iacnt, CalcInfo_);
                     for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
                        if (Cnt[1][Ia_ex]) { found=1; break; }
                        }
                     }
                  }
               }
            if (found) s2_contrib[sblock][cblock] = 1;
            }

         /* does this c block contribute to sigma3? */
         b2brepl(Occs_[sac], Cnt[0], Ij[0], Oij[0], Ridx[0],
            Sgn[0], AlphaG_, sac, cac, nas, CalcInfo_);
         for (Iaidx=0,found=0; Iaidx<S.Ia_size_[sblock]; Iaidx++) {
            if (Cnt[0][Iaidx]) found=1;
            }
         if (found) { /* see if beta is ok */
            found=0;
            b2brepl(Occs_[sbc], Cnt[0], Ij[0], Oij[0], Ridx[0],
               Sgn[0], BetaG_, sbc, cbc, nbs, CalcInfo_);
            for (Ibidx=0; Ibidx<S.Ib_size_[sblock]; Ibidx++) {
               if (Cnt[0][Ibidx]) found=1;
               }
            if (found) s3_contrib[sblock][cblock] = 1;
            }
         } /* end loop over c blocks */
      } /* end loop over sigma blocks */

   if (print_ > 3) {
    outfile->Printf("\nSigma 1:\n");
     for (i=0; i<S.num_blocks_; i++) {
       outfile->Printf( "Contributions to sigma block %d\n", i);
       for (j=0; j<C.num_blocks_; j++) {
         if (s1_contrib[i][j]) outfile->Printf( "%3d ", j);
       }
       outfile->Printf( "\n");
     }

    outfile->Printf("\n\nSigma 2:\n");
     for (i=0; i<S.num_blocks_; i++) {
       outfile->Printf( "Contributions to sigma block %d\n", i);
       for (j=0; j<C.num_blocks_; j++) {
         if (s2_contrib[i][j]) outfile->Printf( "%3d ", j);
       }
       outfile->Printf( "\n");
     }

    outfile->Printf("\n\nSigma 3:\n");
     for (i=0; i<S.num_blocks_; i++) {
       outfile->Printf( "Contributions to sigma block %d\n", i);
       for (j=0; j<C.num_blocks_; j++) {
         if (s3_contrib[i][j]) outfile->Printf( "%3d ", j);
       }
       outfile->Printf( "\n");
     }
   }

}


}} // namespace psi::detci
