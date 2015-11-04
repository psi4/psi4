/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */


#include <cstdio>
#include <cmath>
#include <cstring>
#include <psifiles.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libqt/slaterdset.h>

#include "structs.h"
#include "slaterd.h"
#include "civect.h"
#include "ciwave.h"


namespace psi { namespace detci {


/*
** diag_h(): Function diagonalizes the hamiltonian
**
** Parameters:
**    alplist = list of alpha strings
**    betlist = list of beta strings
**
** Returns: none
*/
void CIWavefunction::diag_h()
{
   BIGINT size;
   int nroots, i, j;
   double conv_rms, conv_e, *evals, **evecs, nucrep, edrc, tval;
   int *tptr;
   double *cbuf;
   char e_label[PSIO_KEYLEN]; /* 80... */

   nroots = Parameters_->num_roots;

   conv_rms = Parameters_->convergence;
   conv_e = Parameters_->energy_convergence;

   if (Parameters_->have_special_conv) {
     tval = sqrt(conv_rms) * 10.0;
     if (tval > conv_e) conv_e = tval;
     if (Parameters_->special_conv > conv_rms)
       conv_rms = Parameters_->special_conv;
   }

   size = CIblks_->vectlen;
   if ((BIGINT) Parameters_->nprint > size) Parameters_->nprint = (int) size;
   nucrep = CalcInfo_->enuc;
   edrc = CalcInfo_->edrc;
   // tmp_ras_array = init_array(1024);

   if (Parameters_->bendazzoli)
      outfile->Printf( "\nBendazzoli algorithm selected for sigma3\n");

   /* Direct Method --- use RSP diagonalization routine */
   if (Parameters_->diag_method == METHOD_RSP) {

      CIvect Cvec(1, 1, 0, 0, CIblks_);
      //CIvect Cvec(CIblks_->vectlen, CIblks_->num_blocks, 1, Parameters_->Ms0,
      //   CIblks_->Ia_code, CIblks_->Ib_code, CIblks_->Ia_size, CIblks_->Ib_size,
      //   CIblks_->offset, CIblks_->num_alp_codes, CIblks_->num_bet_codes,
      //   CalcInfo_->nirreps, AlphaG_->subgr_per_irrep, 1, 0, 0,
      //   CIblks_->first_iablk, CIblks_->last_iablk, CIblks_->decode);
      // shouldn't need to open I/O files for this fake CIvec, unit=0

      double **H, **rsp_evecs;
      int Iarel, Ialist, Ibrel, Iblist;
      BIGINT ii, jj;
      SlaterDeterminant I, J;
      int *mi_iac, *mi_ibc, *mi_iaidx, *mi_ibidx;
      double *mi_coeff;

      if (Parameters_->print_lvl) {
         outfile->Printf( "\nFind all roots with RSP\n") ;
         outfile->Printf( "\n") ;
         }

      /* construct and print one block at a time for debugging */
      /*
      int ii2, jj2, blk, blk2, det1, det2;
      double **Hpart;

      for (blk = 0; blk < CIblks_->num_blocks; blk++) {
        for (blk2 = 0; blk2 < CIblks_->num_blocks; blk2++) {
          Hpart = init_matrix(CIblks_->Ia_size[blk]*CIblks_->Ib_size[blk],
                              CIblks_->Ia_size[blk2]*CIblks_->Ib_size[blk2]);
          for (ii=0,det1=0; ii<CIblks_->Ia_size[blk]; ii++) {
            for (jj=0; jj<CIblks_->Ib_size[blk]; jj++, det1++) {
              I.set(CalcInfo_->num_alp_expl,alplist[CIblks_->Ia_code[blk]][ii].occs,
                   CalcInfo_->num_bet_expl,betlist[CIblks_->Ib_code[blk]][jj].occs);
              for (ii2=0,det2=0; ii2<CIblks_->Ia_size[blk2]; ii2++) {
                for (jj2=0; jj2<CIblks_->Ib_size[blk2]; jj2++,det2++) {
                  J.set(CalcInfo_->num_alp_expl,
                        alplist[CIblks_->Ia_code[blk2]][ii2].occs,
                        CalcInfo_->num_bet_expl,
                        betlist[CIblks_->Ib_code[blk2]][jj2].occs);
                  Hpart[det1][det2] = matrix_element(&I,&J);
                }
              }
            }
          }
          if (Parameters_->print_lvl > 4 && size < 200) {
            outfile->Printf( "\nBlock %d %d of ", blk, blk2);
            outfile->Printf( "Hamiltonian matrix:\n");
            print_mat(Hpart, CIblks_->Ia_size[blk]*CIblks_->Ib_size[blk],
                             CIblks_->Ia_size[blk2]*CIblks_->Ib_size[blk2],
                      outfile);
          }
          free_matrix(Hpart, CIblks_->Ia_size[blk]*CIblks_->Ib_size[blk]);
        }
      }
      */
      /* end block-at-a-time stuff */

      H = init_matrix(size, size);
      rsp_evecs = init_matrix(size, size) ;
      evals = init_array(size) ;
      for (ii=0; ii<size; ii++) {
         Cvec.det2strings(ii, &Ialist, &Iarel, &Iblist, &Ibrel);
         I.set(CalcInfo_->num_alp_expl,
               alplist_[Ialist][Iarel].occs, CalcInfo_->num_bet_expl,
               betlist_[Iblist][Ibrel].occs);
         /* introduce symmetry or other restrictions here */
         for (jj=0; jj<ii; jj++) {
            Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
            J.set(CalcInfo_->num_alp_expl,
               alplist_[Ialist][Iarel].occs, CalcInfo_->num_bet_expl,
               betlist_[Iblist][Ibrel].occs);
            H[ii][jj] = H[jj][ii] = matrix_element(&I, &J);
            }
         Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
         J.set(CalcInfo_->num_alp_expl,
            alplist_[Ialist][Iarel].occs, CalcInfo_->num_bet_expl,
            betlist_[Iblist][Ibrel].occs);
         H[jj][jj] = matrix_element(&J, &J) + CalcInfo_->edrc;
         }

      if (Parameters_->print_lvl > 4 && size < 200) {
         outfile->Printf( "\nHamiltonian matrix:\n");
         print_mat(H, size, size, "outfile");
         }

      sq_rsp(size, size, H, evals, 1, rsp_evecs, 1.0E-14);
      if (Parameters_->print_lvl > 4) {
         eivout(rsp_evecs, evals, size, nroots, "outfile") ;
         }
      evecs = init_matrix(nroots, size);
      for (ii=0; ii<nroots; ii++) {
         for (jj=0; jj<size; jj++) {
            evecs[ii][jj] = rsp_evecs[jj][ii];
            }
         }
      free_matrix(rsp_evecs, size);
      free_matrix(H, size);

      if (Parameters_->print_lvl) {
         mi_iac = init_int_array(Parameters_->nprint);
         mi_ibc = init_int_array(Parameters_->nprint);
         mi_iaidx = init_int_array(Parameters_->nprint);
         mi_ibidx = init_int_array(Parameters_->nprint);
         mi_coeff = init_array(Parameters_->nprint);

         for (i=0; i<nroots; i++) {
            outfile->Printf(
               "\n\n* ROOT %2d CI total energy = %17.13lf\n",
               i+1,evals[i]+nucrep);
            Cvec.setarray(evecs[i], size);
            zero_arr(mi_coeff, Parameters_->nprint);
            Cvec.max_abs_vals(Parameters_->nprint, mi_iac, mi_ibc, mi_iaidx,
               mi_ibidx, mi_coeff, Parameters_->neg_only);
            print_vec(Parameters_->nprint, mi_iac, mi_ibc, mi_iaidx, mi_ibidx,
               mi_coeff);
            }

         free(mi_iac);  free(mi_ibc);
         free(mi_iaidx);  free(mi_ibidx);
         free(mi_coeff);
         }

      //if (Parameters_->write_energy) write_energy(nroots, evals, nucrep);

      /* Dump the vector to a PSIO file
         Added by Edward valeev (June 2002) */
      if (Parameters_->export_ci_vector) {
        StringSet alphastrings, betastrings;
        SlaterDetSet dets;
        SlaterDetVector vec;
        short int *drc_occ;
        unsigned char *newocc;

        if (CalcInfo_->num_drc_orbs > 0) {
          drc_occ = (short int *)
            malloc(CalcInfo_->num_drc_orbs*sizeof(short int));
          for (int l=0; l<CalcInfo_->num_drc_orbs; l++) {
            drc_occ[l] = CalcInfo_->order[l]; /* put it in Pitzer order */
          }
        }

        newocc = (unsigned char *)
          malloc(((AlphaG_->num_el > BetaG_->num_el) ?
            AlphaG_->num_el : BetaG_->num_el)*sizeof(unsigned char));

        stringset_init(&alphastrings,AlphaG_->num_str,AlphaG_->num_el,
                       CalcInfo_->num_drc_orbs, drc_occ);
        int list_gr = 0;
        int offset = 0;
        for(int irrep=0; irrep<AlphaG_->nirreps; irrep++) {
          for(int gr=0; gr<AlphaG_->subgr_per_irrep; gr++,list_gr++) {
            int nlists_per_gr = AlphaG_->sg[irrep][gr].num_strings;
            for(int l=0; l<nlists_per_gr; l++) {
              /* convert occs to Pitzer order */
              for (int n=0; n<AlphaG_->num_el; n++) {
                newocc[n] = (unsigned char)
                  CalcInfo_->order[alplist_[list_gr][l].occs[n] +
                                CalcInfo_->num_drc_orbs];
              }
              stringset_add(&alphastrings,l+offset,newocc);
            }
            offset += nlists_per_gr;
          }
        }

        stringset_init(&betastrings,BetaG_->num_str,BetaG_->num_el,
                       CalcInfo_->num_drc_orbs, drc_occ);
        list_gr = 0;
        offset = 0;
        for(int irrep=0; irrep<BetaG_->nirreps; irrep++) {
          for(int gr=0; gr<BetaG_->subgr_per_irrep; gr++,list_gr++) {
            int nlists_per_gr = BetaG_->sg[irrep][gr].num_strings;
            for(int l=0; l<nlists_per_gr; l++) {
              /* convert occs to Pitzer order */
              for (int n=0; n<BetaG_->num_el; n++) {
                newocc[n] = (unsigned char)
                  CalcInfo_->order[betlist_[list_gr][l].occs[n] +
                                CalcInfo_->num_drc_orbs];
              }
              stringset_add(&betastrings,l+offset,newocc);
            }
            offset += nlists_per_gr;
          }
        }
        free(newocc);
        if (CalcInfo_->num_drc_orbs > 0)
          free(drc_occ);

        int Iarel, Ialist, Ibrel, Iblist;
        // the slaterdetset code below will fail if size > int
        // but that should be ok b/c we won't be running RSP in that case...
        slaterdetset_init(&dets,size,&alphastrings,&betastrings);
        for (int ii=0; ii<size; ii++) {
          Cvec.det2strings(ii, &Ialist, &Iarel, &Iblist, &Ibrel);
          int irrep = Ialist/AlphaG_->subgr_per_irrep;
          int gr = Ialist%AlphaG_->subgr_per_irrep;
          int Ia = Iarel + AlphaG_->list_offset[Ialist];
          irrep = Iblist/BetaG_->subgr_per_irrep;
          gr = Iblist%BetaG_->subgr_per_irrep;
          int Ib = Ibrel + BetaG_->list_offset[Iblist];
          slaterdetset_add(&dets, ii, Ia, Ib);
        }

        slaterdetvector_init(&vec, &dets);
        slaterdetvector_set(&vec, evecs[0]);
        slaterdetvector_write(PSIF_CIVECT,"CI vector",&vec);

        slaterdetvector_delete(&vec);
        slaterdetset_delete(&dets);
        stringset_delete(&alphastrings);
        stringset_delete(&betastrings);
      }
    } /* end RSP section */

   /* RSP test of Davidson/Liu (SEM) diagonalization routine */
   else if (Parameters_->diag_method == METHOD_RSPTEST_OF_SEM) {

      // in-core CIvectors, shouldn't need to open files
      CIvect Cvec(1, 1, 0, 0, CIblks_);
      CIvect Hd(1, 1, 0, 0, CIblks_);
      //CIvect Cvec(CIblks_->vectlen, CIblks_->num_blocks, 1, Parameters_->Ms0,
      //   CIblks_->Ia_code, CIblks_->Ib_code, CIblks_->Ia_size, CIblks_->Ib_size,
      //   CIblks_->offset, CIblks_->num_alp_codes, CIblks_->num_bet_codes,
      //   CalcInfo_->nirreps, AlphaG_->subgr_per_irrep, 1, 0, 0,
      //   CIblks_->first_iablk, CIblks_->last_iablk, CIblks_->decode);
      //CIvect Hd(CIblks_->vectlen, CIblks_->num_blocks, 1, Parameters_->Ms0,
      //   CIblks_->Ia_code, CIblks_->Ib_code, CIblks_->Ia_size, CIblks_->Ib_size,
      //   CIblks_->offset, CIblks_->num_alp_codes, CIblks_->num_bet_codes,
      //   CalcInfo_->nirreps, AlphaG_->subgr_per_irrep, 1, 0, 0,
      //   CIblks_->first_iablk, CIblks_->last_iablk, CIblks_->decode);

      double **H, **b;
      int Ia, Ib, Iarel, Ialist, Ibrel, Iblist, ij, k, l, tmpi, L;
      unsigned long int ii, jj;
      SlaterDeterminant I, J;
      int *mi_iac, *mi_ibc, *mi_iaidx, *mi_ibidx;
      double *mi_coeff;
      int sm_tridim;
      double *sm_evals, *sm_mat, **sm_evecs, tval;

      if (Parameters_->print_lvl) {
         outfile->Printf( "\nFind the roots by the SEM Test Method\n");
         outfile->Printf( "(n.b. this is for debugging purposes only!)\n");
         //outfile->Printf( "Energy convergence = %3g\n", conv_e);
         //outfile->Printf( "RMS CI vector convergence = %3g\n\n", conv_rms);
         }

      H0block_init(size);

      /* get the diagonal elements of H into an array Hd */

      Hd.diag_mat_els(alplist_, betlist_, CalcInfo_->onel_ints,
         CalcInfo_->twoel_ints, edrc, CalcInfo_->num_alp_expl,
         CalcInfo_->num_bet_expl, CalcInfo_->num_ci_orbs, Parameters_->hd_ave);

      /* get the biggest elements and put in H0block */
      if (H0block_->size) {
         Hd.max_abs_vals(H0block_->size, H0block_->alplist, H0block_->betlist,
            H0block_->alpidx, H0block_->betidx, H0block_->H00, Parameters_->neg_only);
         }

    /* MLL added this line 5-21-98 */
    /*
      if (Parameters_->hd_otf) rclose(Parameters_->first_hd_tmp_unit,4);
    */

      H0block_setup(CIblks_->num_blocks, CIblks_->Ia_code, CIblks_->Ib_code);
      if (Parameters_->hd_ave) {
        H0block_spin_cpl_chk();
         if (H0block_->osize - H0block_->size) {
            outfile->Printf("H0block size reduced by %d to %d to ensure"
             "completion of spin-coupling sets\n",
             (H0block_->osize - H0block_->size), H0block_->size);

            }
        }
      if (Parameters_->Ms0) {
         H0block_pairup(0);
         if (H0block_->osize - H0block_->size) {
            outfile->Printf("H0block size reduced by %d to ensure pairing.\n",
               (H0block_->osize - H0block_->size));

            }
         }

      if (Parameters_->print_lvl > 4 && Parameters_->hd_otf == FALSE) {
         outfile->Printf( "\nDiagonal elements of the Hamiltonian\n");
         Hd.print("outfile");
         }


      if (H0block_->size) {
         H0block_fill();
         }

      if (Parameters_->print_lvl > 2 && H0block_->size) {
         H0block_print();
         }

      if (Parameters_->print_lvl > 3 && H0block_->size) {
         outfile->Printf( "\n\nH0 Block:\n");
         print_mat(H0block_->H0b, H0block_->size, H0block_->size, "outfile");
         }

      H = init_matrix(size, size);
      evals = init_array(size) ;
      for (ii=0; ii<size; ii++) {
         Cvec.det2strings(ii, &Ialist, &Iarel, &Iblist, &Ibrel);
         I.set(CalcInfo_->num_alp_expl,
               alplist_[Ialist][Iarel].occs, CalcInfo_->num_bet_expl,
               betlist_[Iblist][Ibrel].occs);
         /* introduce symmetry or other restrictions here */
         for (jj=0; jj<ii; jj++) {
            Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
            J.set(CalcInfo_->num_alp_expl,
               alplist_[Ialist][Iarel].occs, CalcInfo_->num_bet_expl,
               betlist_[Iblist][Ibrel].occs);
            H[ii][jj] = H[jj][ii] = matrix_element(&I, &J);
            }
         Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
         J.set(CalcInfo_->num_alp_expl,
            alplist_[Ialist][Iarel].occs, CalcInfo_->num_bet_expl,
            betlist_[Iblist][Ibrel].occs);
         H[jj][jj] = matrix_element(&J, &J) + CalcInfo_->edrc;
         }

      /* obtain a set of L orthonormal trial vectors, L > nroots */
      b = (double **) malloc (Parameters_->maxnvect * sizeof(double *)) ;
      for (i=0; i<Parameters_->maxnvect; i++) {
         if (i<Parameters_->num_init_vecs) b[i] = init_array(size) ;
         else b[i] = NULL ;
         }

      evecs = init_matrix(Parameters_->num_roots, size);

      L = H0block_->size;
      sm_tridim = L * (L + 1) / 2 ;
      sm_mat = init_array(sm_tridim) ;
      sm_evals = init_array(L) ;
      sm_evecs = init_matrix(L, L) ;
      for (i=0, ij=0; i<L; i++)
         for (j=0; j<=i; j++, ij++)
            sm_mat[ij] = H0block_->H0b[i][j] ;
      rsp(L, L, sm_tridim, sm_mat, sm_evals, 1, sm_evecs, 1E-14) ;

      /*
      if (Parameters_->precon == PRECON_GEN_DAVIDSON) {
        for (i=0; i<H0block_->size; i++) {
           H0block_->H0b_eigvals[i] = sm_evals[i];
           for (j=0; j<H0block_->size; i++)
              H0block_->H0b_diag[i][j] = sm_evecs[i][j];
           }
        }
      */

      /* need to fill out sm_evecs into b (pad w/ 0's) */
      cbuf = *(Cvec.blockptr(0));
      Cvec.buf_unlock();
      for (i=0,k=0; i<L && k < Parameters_->num_init_vecs; i++) {

         /* check sm_evecs[i] to see if it has the correct spin symmetry */
         for (j=0,tmpi=0; Parameters_->Ms0 && j<L && !tmpi; j++) {
            l = H0block_->pair[j];
            if (l == -1) {
               throw PsiException("(diag_h sem_test): unpaired H0block member!",__FILE__,__LINE__);
               }
            tval = sm_evecs[l][i];
            if ((int) Parameters_->S%2) tval = -tval;
            if (sm_evecs[j][i] - tval > 1.0E-12) tmpi=1;
            }
         if (tmpi) continue;

         for (j=0; j<L; j++) sm_evals[j] = sm_evecs[j][i];

         Cvec.buf_lock(b[k]);
         Cvec.init_vals(0, L, H0block_->alplist, H0block_->alpidx,
            H0block_->betlist, H0block_->betidx, H0block_->blknum, sm_evals);
         Cvec.buf_unlock();
         k++;
         }
      Cvec.buf_lock(cbuf);
      free(sm_mat);
      free(sm_evals);
      free_matrix(sm_evecs, L);
      for (i=k; i<Parameters_->num_init_vecs; i++) free(b[i]);
      L = k;
      if (L < Parameters_->num_roots) {
         throw PsiException("(diag_h sem_test): Ooops! L < num_roots!",__FILE__,__LINE__);
         }
      sem_test(H, size, Parameters_->num_roots, L, evecs, evals, b, conv_e,
         conv_rms, Parameters_->maxiter, (nucrep+CalcInfo_->edrc), &i,
         Parameters_->maxnvect);

      outfile->Printf( "SEM used %d expansion vectors\n", i);

      if (Parameters_->print_lvl > 4) {
         eivout(evecs, evals, size, nroots, "outfile") ;
         }
      free_matrix(H, size);

      if (Parameters_->print_lvl) {
         mi_iac = init_int_array(Parameters_->nprint);
         mi_ibc = init_int_array(Parameters_->nprint);
         mi_iaidx = init_int_array(Parameters_->nprint);
         mi_ibidx = init_int_array(Parameters_->nprint);
         mi_coeff = init_array(Parameters_->nprint);

         for (i=0; i<nroots; i++) {
            outfile->Printf(
               "\n\n* ROOT %2d CI total energy = %17.13lf\n",
               i+1,evals[i]+nucrep);
            Cvec.setarray(evecs[i], size);
            zero_arr(mi_coeff, Parameters_->nprint);
            Cvec.max_abs_vals(Parameters_->nprint, mi_iac, mi_ibc, mi_iaidx,
               mi_ibidx, mi_coeff, Parameters_->neg_only);
            print_vec(Parameters_->nprint, mi_iac, mi_ibc, mi_iaidx, mi_ibidx,
               mi_coeff);
            }
         free(mi_iac);  free(mi_ibc);
         free(mi_iaidx);  free(mi_ibidx);
         free(mi_coeff);
         free_matrix(evecs, Parameters_->num_roots);
         }

      //if (Parameters_->write_energy) write_energy(nroots, evals, nucrep);

      } /* end Test of Davidson/Liu section */


   /*
    * Davidson/Liu Simultaneous Expansion Method OR
    * Mitrushenkov's Olsen-modified Davidson Algorithm
    */

   else {

      /* prepare the H0 block */
      H0block_init(size);

      CIvect Hd(Parameters_->icore, 1, Parameters_->num_hd_tmp_units,
                Parameters_->first_hd_tmp_unit, CIblks_);
      //CIvect Hd(CIblks_->vectlen, CIblks_->num_blocks, Parameters_->icore,
      //   Parameters_->Ms0, CIblks_->Ia_code, CIblks_->Ib_code, CIblks_->Ia_size,
      //   CIblks_->Ib_size, CIblks_->offset, CIblks_->num_alp_codes,
      //   CIblks_->num_bet_codes, CalcInfo_->nirreps, AlphaG_->subgr_per_irrep, 1,
      //   Parameters_->num_hd_tmp_units, Parameters_->first_hd_tmp_unit,
      //   CIblks_->first_iablk, CIblks_->last_iablk, CIblks_->decode);

      bool open_old = false;
      if (Parameters_->restart) open_old = true;
      Hd.init_io_files(open_old);

      /* get the diagonal elements of H into an array Hd */
      if (!Parameters_->restart || (Parameters_->restart && Parameters_->hd_otf)) {
         if (Parameters_->print_lvl > 1) {
            outfile->Printf( "\nForming diagonal elements of H\n");

           }
         Hd.diag_mat_els(alplist_, betlist_, CalcInfo_->onel_ints,
            CalcInfo_->twoel_ints, edrc, CalcInfo_->num_alp_expl,
            CalcInfo_->num_bet_expl, CalcInfo_->num_ci_orbs, Parameters_->hd_ave);
         }
      else {
         Hd.read(0,0);
         }

      /* get the biggest elements and put in H0block */
      if (H0block_->size) {

         if (Parameters_->print_lvl > 1) {
            outfile->Printf( "\nForming H0 block\n");

           }

         if (!Parameters_->hd_otf)
           Hd.max_abs_vals(H0block_->size+H0block_->coupling_size,
              H0block_->alplist, H0block_->betlist, H0block_->alpidx,
              H0block_->betidx, H0block_->H00, Parameters_->neg_only);
         }

      //if (Parameters_->hd_otf) rclose(Parameters_->first_hd_tmp_unit,4);
      if (Parameters_->hd_otf) psio_close(Parameters_->first_hd_tmp_unit,1);

      H0block_setup(CIblks_->num_blocks, CIblks_->Ia_code, CIblks_->Ib_code);
      if (Parameters_->filter_guess) H0block_filter_setup();
      if (Parameters_->hd_ave) {
        H0block_spin_cpl_chk();
         if ((H0block_->osize - H0block_->size) && Parameters_->print_lvl > 1) {
            outfile->Printf("H0block size reduced by %d to ensure "
             "completion of spin-coupling sets\n",
             (H0block_->osize - H0block_->size));
            H0block_->osize = H0block_->size;
            }
         if ((H0block_->oguess_size - H0block_->guess_size) &&
             Parameters_->print_lvl > 1) {
           outfile->Printf("H0block guess size reduced by %d to ensure "
             "completion of spin-coupling sets\n",
             (H0block_->oguess_size - H0block_->guess_size));
            H0block_->oguess_size = H0block_->guess_size;
           }
         if ((H0block_->ocoupling_size - H0block_->coupling_size) &&
             Parameters_->print_lvl > 1) {
           outfile->Printf("H0block coupling size reduced by %d to ensure "
             "completion of spin-coupling sets\n",
             (H0block_->ocoupling_size - H0block_->coupling_size));
            H0block_->ocoupling_size = H0block_->coupling_size;
           }

        }
      if (Parameters_->Ms0) {
         /* if (H0block_->guess_size < H0block_->size) */
         H0block_pairup(0); /* pairup h0block size */
         H0block_pairup(1); /* pairup guess_size */
         H0block_pairup(2); /* pairup coupling size */
         if ((H0block_->osize - H0block_->size) && Parameters_->print_lvl > 1) {
           outfile->Printf("H0block size reduced by %d to ensure pairing"
               "and spin-coupling.\n", (H0block_->osize - H0block_->size));
            }
         if ((H0block_->oguess_size - H0block_->guess_size) &&
             Parameters_->print_lvl > 1) {
           outfile->Printf("H0block guess size reduced by %d to "
               "ensure pairing and spin-coupling.\n",
               (H0block_->oguess_size - H0block_->guess_size));
            }
         if ((H0block_->ocoupling_size - H0block_->coupling_size) &&
             Parameters_->print_lvl > 1) {
           outfile->Printf("H0block coupling size reduced by %d to "
               "ensure pairing and spin-coupling.\n",
               (H0block_->ocoupling_size - H0block_->coupling_size));
            }

         }

      Parameters_->neg_only = 0; /* MLL 7-2-97 */
      if (Parameters_->print_lvl > 4) {
         outfile->Printf( "\nDiagonal elements of the Hamiltonian\n");
         Hd.print("outfile");
         }

      if (H0block_->size) {
         H0block_fill();
         }

      if (Parameters_->print_lvl > 2 && H0block_->size) {
         H0block_print();
         }

      if (Parameters_->print_lvl > 3 && H0block_->size) {
         outfile->Printf( "\n\nH0 Block:\n");
         print_mat(H0block_->H0b, H0block_->size, H0block_->size, "outfile");
         }

      /* Davidson/Liu Simultaneous Expansion Method */
      if (Parameters_->diag_method == METHOD_DAVIDSON_LIU_SEM) {

         if (Parameters_->print_lvl) {
            outfile->Printf(
               "\nFind the roots by the Simultaneous Expansion Method ");
            outfile->Printf( "(Block Davidson Method)\n");
            //outfile->Printf( "Energy convergence = %3g\n", conv_e);
            //outfile->Printf( "RMS CI vector convergence = %3g\n\n", conv_rms);
            }

         evals = init_array(nroots);

         sem_iter(Hd, alplist_, betlist_, evals, conv_e, conv_rms,
            nucrep, edrc, nroots, Parameters_->maxiter,
            Parameters_->maxnvect, "outfile", Parameters_->print_lvl);
         }

      /* Mitrushenkov's Olsen Method */
      else {

         if (Parameters_->print_lvl) {
            if (Parameters_->diag_method == METHOD_MITRUSHENKOV)
              outfile->Printf(
                "\nFind the roots with Mitrushenkov's two vector algorithm\n");
            else if (Parameters_->diag_method == METHOD_OLSEN)
              outfile->Printf(
                "\nFind the roots with Olsen's single vector algorithm\n");
            //outfile->Printf( "Energy convergence = %3g\n", conv_e);
            //outfile->Printf( "RMS CI vector convergence = %3g\n", conv_rms);
            }

         evals = init_array(nroots);

         mitrush_iter(Hd, alplist_, betlist_, nroots, evals, conv_rms, conv_e,
            nucrep, edrc, Parameters_->maxiter, Parameters_->maxnvect, "outfile",
            Parameters_->print_lvl);

         H0block_free();
         }

      //if (Parameters_->write_energy) write_energy(nroots, evals, nucrep+edrc);

      } /* end the Davidson-Liu/Mitrushenkov-Olsen-Davidson section */

   /* write the CI energy to PSIF_CHKPT: later fix this to loop over roots */
   chkpt_init(PSIO_OPEN_OLD);
   tval = evals[Parameters_->root]+edrc+nucrep;
   chkpt_wt_etot(tval);

   Process::environment.globals["CURRENT ENERGY"] = tval;
   Process::environment.globals["CURRENT CORRELATION ENERGY"] = tval - CalcInfo_->escf;
   Process::environment.globals["CURRENT REFERENCE ENERGY"] = CalcInfo_->escf;
   Process::environment.globals["CI TOTAL ENERGY"] = tval;
   // eref seems wrong for open shells so replace it with escf below
   // until I fix it ---CDS 11/5/11
   Process::environment.globals["CI CORRELATION ENERGY"] = tval - CalcInfo_->escf;

   if (Parameters_->fci) {
     Process::environment.globals["FCI TOTAL ENERGY"] = tval;
     Process::environment.globals["FCI CORRELATION ENERGY"] = tval - CalcInfo_->escf;
   }
   else {
     if (Parameters_->ex_lvl == 2) {
       Process::environment.globals["CISD TOTAL ENERGY"] = tval;
       Process::environment.globals["CISD CORRELATION ENERGY"] = tval - CalcInfo_->escf;
     }
     else if (Parameters_->ex_lvl == 3) {
       Process::environment.globals["CISDT TOTAL ENERGY"] = tval;
       Process::environment.globals["CISDT CORRELATION ENERGY"] = tval - CalcInfo_->escf;
     }
     else if (Parameters_->ex_lvl == 4) {
       Process::environment.globals["CISDTQ TOTAL ENERGY"] = tval;
       Process::environment.globals["CISDTQ CORRELATION ENERGY"] = tval - CalcInfo_->escf;
     }
     else {
       /*- Process::environment.globals["CIn TOTAL ENERGY"] -*/
       /*- Process::environment.globals["CIn CORRELATION ENERGY"] -*/
       std::stringstream s;
       s << "CI" << Parameters_->ex_lvl << " TOTAL ENERGY";
       Process::environment.globals[s.str()] = tval;
       s.str(std::string());
       s << "CI" << Parameters_->ex_lvl << " CORRELATION ENERGY";
       Process::environment.globals[s.str()] = tval - CalcInfo_->escf;
     }
   }

   for (i=0; i<nroots; i++) {
     sprintf(e_label,"Root %2d energy",i);
     tval = evals[i]+edrc+nucrep;
     chkpt_wt_e_labeled(e_label, tval);

     /*- Process::environment.globals["CI ROOT n TOTAL ENERGY"] -*/
     /*- Process::environment.globals["CI ROOT n CORRELATION ENERGY"] -*/
     std::stringstream s;
     s << "CI ROOT " << (i+1) << " TOTAL ENERGY";
     Process::environment.globals[s.str()] = tval;
     s.str(std::string());
     s << "CI ROOT " << (i+1) << " CORRELATION ENERGY";
     // eref seems wrong for open shells so replace it with escf below
     // until I fix it ---CDS 11/5/11
     Process::environment.globals[s.str()] = tval - CalcInfo_->escf;
   }

   if (Parameters_->average_num > 1) {
     tval = 0.0;
     for (i=0; i<Parameters_->average_num; i++)
       tval += Parameters_->average_weights[i] *
               (edrc+nucrep+evals[Parameters_->average_states[i]]);
     chkpt_wt_e_labeled("State averaged energy",tval);
     Process::environment.globals["CI STATE-AVERAGED TOTAL ENERGY"] = tval;
     // eref seems wrong for open shells so replace it with escf below
     // until I fix it ---CDS 11/5/11
     Process::environment.globals["CI STATE-AVERAGED CORRELATION ENERGY"] =
       tval - CalcInfo_->escf;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] =
       Process::environment.globals["CI STATE-AVERAGED CORRELATION ENERGY"];
   }

   // Set the energy as MCSCF would find it
   if (Parameters_->average_num > 1) { // state average
     Process::environment.globals["MCSCF TOTAL ENERGY"] =
       Process::environment.globals["CI STATE-AVERAGED TOTAL ENERGY"];
   }
   else if (Parameters_->root != 0) { // follow some specific root != lowest
     std::stringstream s;
     s << "CI ROOT " << (Parameters_->root+1) << " TOTAL ENERGY";
     Process::environment.globals["MCSCF TOTAL ENERGY"] =
      Process::environment.globals[s.str()];
   }
   else {
     Process::environment.globals["MCSCF TOTAL ENERGY"] =
      Process::environment.globals["CI TOTAL ENERGY"];
   }

   chkpt_close();

}



}} // namespace psi::detci


