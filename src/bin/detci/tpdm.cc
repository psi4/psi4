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

/*! \file
**  \ingroup DETCI
**  \brief Compute the two-particle density matrix (TPDM)
**
**  C. David Sherrill
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
/* may no longer need #include <libc.h> */
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libmints/mints.h>
#include "structs.h"
#include "civect.h"
#include "ciwave.h"

namespace psi { namespace detci {

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

// DGAS this is still awkward, I think the TPDM code can be less general than the OPDM one for now.
void CIWavefunction::form_tpdm(void)
{
  std::vector<SharedVector> tpdm_list;
  tpdm_list = tpdm(Parameters_->num_roots, Parameters_->first_d_tmp_unit, Parameters_->first_d_tmp_unit);
  tpdm_called_ = true;
  tpdm_aa_ = tpdm_list[0];
  tpdm_ab_ = tpdm_list[1];
  tpdm_bb_ = tpdm_list[2];
  tpdm_    = tpdm_list[3];
}

// We always return the state-averaged tpdm here, not sure what else we really need.
std::vector<SharedVector> CIWavefunction::tpdm(int nroots, int Ifirstunit, int Jfirstunit) 
{

   //CIvect Ivec, Jvec;
   struct iwlbuf TBuff;
   struct iwlbuf TBuff_aa;
   struct iwlbuf TBuff_bb;
   struct iwlbuf TBuff_ab;
   int i, j, k, l, lmax, ij, kl, ijkl, ijksym;
   int i2, j2, k2, l2, ndrc, populated_orbs;
   int *orbsym;
   int maxrows, maxcols, ntri, ntri2;
   unsigned long bufsz;
   double **transp_tmp = NULL;
   double **transp_tmp2 = NULL;
   double *buffer1, *buffer2, value;
   double **onepdm_a, **onepdm_b;
   double *twopdm_aa, *twopdm_bb, *twopdm_ab;
   int Iroot, Jroot;
   int Iblock, Iblock2, Ibuf, Iac, Ibc, Inas, Inbs, Iairr;
   int Jblock, Jblock2, Jbuf, Jac, Jbc, Jnas, Jnbs, Jairr;
   int do_Jblock, do_Jblock2;
   char opdm_key[80];
   int root_idx;   /* what root we're on */
   double weight;  /* the weight of that root */

   ndrc = CalcInfo_->num_drc_orbs;
   populated_orbs = CalcInfo_->nmo - CalcInfo_->num_drv_orbs;
   int writeflag = 1;
   int printflag = 0;
   int targetfile = Parameters_->tpdm_file;


    
   CIvect Ivec(Parameters_->icore, nroots, 1, Ifirstunit, CIblks_, CalcInfo_, Parameters_,
               H0block_, false);
   Ivec.init_io_files(true);

   CIvect Jvec(Parameters_->icore, nroots, 1, Jfirstunit, CIblks_, CalcInfo_, Parameters_,
               H0block_, false);
   Jvec.init_io_files(true);

   buffer1 = Ivec.buf_malloc();
   buffer2 = Jvec.buf_malloc();
   Ivec.buf_lock(buffer1);
   Jvec.buf_lock(buffer2);

   ntri = CalcInfo_->num_ci_orbs * CalcInfo_->num_ci_orbs;
   ntri2 = (ntri * (ntri + 1)) / 2;
   SharedVector tpdm_aa(new Vector("MO-basis TPDM AA", ntri2));
   SharedVector tpdm_ab(new Vector("MO-basis TPDM AB", ntri*ntri));
   SharedVector tpdm_bb(new Vector("MO-basis TPDM BB", ntri2));
   double* tpdm_aap = tpdm_aa->pointer();
   double* tpdm_abp = tpdm_ab->pointer();
   double* tpdm_bbp = tpdm_bb->pointer();

   if ((Ivec.icore_==2 && Ivec.Ms0_ && CalcInfo_->ref_sym != 0) || 
       (Ivec.icore_==0 && Ivec.Ms0_)) {
     for (i=0, maxrows=0, maxcols=0; i<Ivec.num_blocks_; i++) {
       if (Ivec.Ia_size_[i] > maxrows) maxrows = Ivec.Ia_size_[i];
       if (Ivec.Ib_size_[i] > maxcols) maxcols = Ivec.Ib_size_[i];
     }
     if (maxcols > maxrows) maxrows = maxcols;
     transp_tmp = (double **) malloc (maxrows * sizeof(double *));
     transp_tmp2 = (double **) malloc (maxrows * sizeof(double *));
     if (transp_tmp == NULL || transp_tmp2 == NULL) {
       printf("(tpdm): Trouble with malloc'ing transp_tmp\n");
     }
     bufsz = Ivec.get_max_blk_size();
     transp_tmp[0] = init_array(bufsz);
     transp_tmp2[0] = init_array(bufsz);
     if (transp_tmp[0] == NULL || transp_tmp2[0] == NULL) {
       printf("(tpdm): Trouble with malloc'ing transp_tmp[0]\n");
     }
   }


   if (Parameters_->icore == 0) {

     /* loop over all the roots requested */
     for (root_idx=0; root_idx<Parameters_->average_num; root_idx++)  {
       Iroot = Parameters_->average_states[root_idx];
       weight = Parameters_->average_weights[root_idx];
       Jroot = Iroot;  /* change later if need transition matrix elements */

       for (Ibuf=0; Ibuf<Ivec.buf_per_vect_; Ibuf++) {
         Ivec.read(Iroot, Ibuf);
         Iblock = Ivec.buf2blk_[Ibuf];
         Iac = Ivec.Ia_code_[Iblock];
         Ibc = Ivec.Ib_code_[Iblock];
         Inas = Ivec.Ia_size_[Iblock];
         Inbs = Ivec.Ib_size_[Iblock];
       
         for (Jbuf=0; Jbuf<Jvec.buf_per_vect_; Jbuf++) {
           do_Jblock=0; do_Jblock2=0;
           Jblock = Jvec.buf2blk_[Jbuf];
           Jblock2 = -1;
           Jac = Jvec.Ia_code_[Jblock];
           Jbc = Jvec.Ib_code_[Jblock];
           if (Jvec.Ms0_) Jblock2 = Jvec.decode_[Jbc][Jac];
             Jnas = Jvec.Ia_size_[Jblock];
             Jnbs = Jvec.Ib_size_[Jblock];
             if (s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock]
               || s3_contrib_[Iblock][Jblock]) 
             do_Jblock = 1;
           if (Jvec.buf_offdiag_[Jbuf] && (s1_contrib_[Iblock][Jblock2] ||
             s2_contrib_[Iblock][Jblock2] ||
             s3_contrib_[Iblock][Jblock2]))
             do_Jblock2 = 1;
           if (!do_Jblock && !do_Jblock2) continue;

           Jvec.read(Jroot, Jbuf);
	 
           if (do_Jblock) {
             tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs, 
               Ivec.num_alpcodes_, Ivec.num_betcodes_, tpdm_aap, tpdm_bbp, tpdm_abp, 
               Jvec.blocks_[Jblock], Ivec.blocks_[Iblock], 
               Jac, Jbc, Jnas, Jnbs, Iac, Ibc, Inas, Inbs, weight);
           }
	 
           if (do_Jblock2) {
             Jvec.transp_block(Jblock, transp_tmp);
             tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
               Ivec.num_alpcodes_, Ivec.num_betcodes_, 
               tpdm_aap, tpdm_bbp, tpdm_abp, transp_tmp, Ivec.blocks_[Iblock], 
               Jbc, Jac, Jnbs, Jnas, Iac, Ibc, Inas, Inbs, weight);
           }
	 
         } /* end loop over Jbuf */
       
         if (Ivec.buf_offdiag_[Ibuf]) { /* need to get contrib of transpose */
           Iblock2 = Ivec.decode_[Ibc][Iac];
           Iac = Ivec.Ia_code_[Iblock2];
           Ibc = Ivec.Ib_code_[Iblock2];
           Inas = Ivec.Ia_size_[Iblock2];
           Inbs = Ivec.Ib_size_[Iblock2];
       
           Ivec.transp_block(Iblock, transp_tmp2);

           for (Jbuf=0; Jbuf<Jvec.buf_per_vect_; Jbuf++) {
	     do_Jblock=0; do_Jblock2=0;
	     Jblock = Jvec.buf2blk_[Jbuf];
	     Jblock2 = -1;
	     Jac = Jvec.Ia_code_[Jblock];
	     Jbc = Jvec.Ib_code_[Jblock];
	     if (Jvec.Ms0_) Jblock2 = Jvec.decode_[Jbc][Jac];
	     Jnas = Jvec.Ia_size_[Jblock];
	     Jnbs = Jvec.Ib_size_[Jblock];
	     if (s1_contrib_[Iblock2][Jblock] || s2_contrib_[Iblock2][Jblock] ||
	         s3_contrib_[Iblock2][Jblock]) 
	       do_Jblock = 1;
	     if (Jvec.buf_offdiag_[Jbuf] && (s1_contrib_[Iblock2][Jblock2] ||
               s2_contrib_[Iblock2][Jblock2] ||
               s3_contrib_[Iblock2][Jblock2]))
               do_Jblock2 = 1;
             if (!do_Jblock && !do_Jblock2) continue;
	   
             Jvec.read(Jroot, Jbuf);
	 
             if (do_Jblock) {
               tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs, 
                 Ivec.num_alpcodes_, Ivec.num_betcodes_, 
                 tpdm_aap, tpdm_bbp, tpdm_abp, Jvec.blocks_[Jblock], 
                 transp_tmp2, Jac, Jbc, Jnas,
                 Jnbs, Iac, Ibc, Inas, Inbs, weight);
             }
	   
             if (do_Jblock2) {
               Jvec.transp_block(Jblock, transp_tmp);
               tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                 Ivec.num_alpcodes_, Ivec.num_betcodes_,
                 tpdm_aap, tpdm_bbp, tpdm_abp, transp_tmp, transp_tmp2, 
                 Jbc, Jac, Jnbs, Jnas, Iac, Ibc, Inas, Inbs, weight);
             }
           } /* end loop over Jbuf */

         } /* end loop over Ibuf transpose */
       } /* end loop over Ibuf */
     } /* end loop over roots */
   } /* end icore==0 */

   else if (Parameters_->icore==1) { /* whole vectors in-core */

     for (root_idx=0; root_idx<Parameters_->average_num; root_idx++)  {
       Iroot = Parameters_->average_states[root_idx];
       weight = Parameters_->average_weights[root_idx];
       Jroot = Iroot;  /* change later if need transition matrix elements */

       Ivec.read(Iroot, 0);
       Jvec.read(Jroot, 0);
       for (Iblock=0; Iblock<Ivec.num_blocks_; Iblock++) {
         Iac = Ivec.Ia_code_[Iblock];
         Ibc = Ivec.Ib_code_[Iblock];
         Inas = Ivec.Ia_size_[Iblock];
         Inbs = Ivec.Ib_size_[Iblock];
         if (Inas==0 || Inbs==0) continue;
         for (Jblock=0; Jblock<Jvec.num_blocks_; Jblock++) {
           Jac = Jvec.Ia_code_[Jblock];
           Jbc = Jvec.Ib_code_[Jblock];
           Jnas = Jvec.Ia_size_[Jblock];
           Jnbs = Jvec.Ib_size_[Jblock];
           if (s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock] ||
             s3_contrib_[Iblock][Jblock])
             tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
               Ivec.num_alpcodes_, Ivec.num_betcodes_, 
               tpdm_aap, tpdm_bbp, tpdm_abp, Jvec.blocks_[Jblock], Ivec.blocks_[Iblock], 
               Jac, Jbc, Jnas, Jnbs, Iac, Ibc, Inas, Inbs, weight);
         }
       } /* end loop over Iblock */
     } /* end loop over roots */
   } /* end icore==1 */

   else if (Parameters_->icore==2) { /* icore==2 */
     for (root_idx=0; root_idx<Parameters_->average_num; root_idx++)  {
       Iroot = Parameters_->average_states[root_idx];
       weight = Parameters_->average_weights[root_idx];
       Jroot = Iroot;  /* change later if need transition matrix elements */

       for (Ibuf=0; Ibuf<Ivec.buf_per_vect_; Ibuf++) {
         Ivec.read(Iroot, Ibuf);
         Iairr = Ivec.buf2blk_[Ibuf];

         for (Jbuf=0; Jbuf<Jvec.buf_per_vect_; Jbuf++) {
           Jvec.read(Jroot, Jbuf);
           Jairr = Jvec.buf2blk_[Jbuf];

           for (Iblock=Ivec.first_ablk_[Iairr]; Iblock<=Ivec.last_ablk_[Iairr];
             Iblock++) {
             Iac = Ivec.Ia_code_[Iblock];
             Ibc = Ivec.Ib_code_[Iblock];
             Inas = Ivec.Ia_size_[Iblock];
             Inbs = Ivec.Ib_size_[Iblock];
   
             for (Jblock=Jvec.first_ablk_[Jairr]; Jblock<=Jvec.last_ablk_[Jairr];
               Jblock++) {
               Jac = Jvec.Ia_code_[Jblock];
               Jbc = Jvec.Ib_code_[Jblock];
               Jnas = Jvec.Ia_size_[Jblock];
               Jnbs = Jvec.Ib_size_[Jblock];
   
               if (s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock] ||
                 s3_contrib_[Iblock][Jblock])
               tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs, 
                 Ivec.num_alpcodes_, Ivec.num_betcodes_,
                 tpdm_aap, tpdm_bbp, tpdm_abp, Jvec.blocks_[Jblock], Ivec.blocks_[Iblock], 
                 Jac, Jbc, Jnas, Jnbs, Iac, Ibc, Inas, Inbs, weight);

               if (Jvec.buf_offdiag_[Jbuf]) {
                 Jblock2 = Jvec.decode_[Jbc][Jac];
                 if (s1_contrib_[Iblock][Jblock2] ||
                   s2_contrib_[Iblock][Jblock2] ||
                   s3_contrib_[Iblock][Jblock2]) {
                   Jvec.transp_block(Jblock, transp_tmp);
                   tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                     Ivec.num_alpcodes_, Ivec.num_betcodes_,
                     tpdm_aap, tpdm_bbp, tpdm_abp, transp_tmp, Ivec.blocks_[Iblock], 
                     Jbc, Jac, Jnbs, Jnas, Iac, Ibc, Inas, Inbs, weight);
                 }
               }

             } /* end loop over Jblock */

             if (Ivec.buf_offdiag_[Ibuf]) {
               Iblock2 = Ivec.decode_[Ibc][Iac];
               Ivec.transp_block(Iblock, transp_tmp2);
               Iac = Ivec.Ia_code_[Iblock2];
               Ibc = Ivec.Ib_code_[Iblock2];
               Inas = Ivec.Ia_size_[Iblock2];
               Inbs = Ivec.Ib_size_[Iblock2];
	   
               for (Jblock=Jvec.first_ablk_[Jairr]; 
                 Jblock<=Jvec.last_ablk_[Jairr]; Jblock++) {
                 Jac = Jvec.Ia_code_[Jblock];
                 Jbc = Jvec.Ib_code_[Jblock];
                 Jnas = Jvec.Ia_size_[Jblock];
                 Jnbs = Jvec.Ib_size_[Jblock];
	   
                 if (s1_contrib_[Iblock2][Jblock] || s2_contrib_[Iblock2][Jblock]
                   || s3_contrib_[Iblock2][Jblock])
                   tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                     Ivec.num_alpcodes_, Ivec.num_betcodes_, 
                     tpdm_aap, tpdm_bbp, tpdm_abp, Jvec.blocks_[Jblock], transp_tmp2, 
                     Jac, Jbc, Jnas, Jnbs, Iac, Ibc, Inas, Inbs, weight);

                 if (Jvec.buf_offdiag_[Jbuf]) {
                   Jblock2 = Jvec.decode_[Jbc][Jac];
                   if (s1_contrib_[Iblock][Jblock2] ||
                     s2_contrib_[Iblock][Jblock2] ||
                     s3_contrib_[Iblock][Jblock2]) {
                     Jvec.transp_block(Jblock, transp_tmp);
                     tpdm_block(alplist_, betlist_, CalcInfo_->num_ci_orbs,
                     Ivec.num_alpcodes_, Ivec.num_betcodes_,
                     tpdm_aap, tpdm_bbp, tpdm_abp, transp_tmp, transp_tmp2, Jbc, Jac,
                     Jnbs, Jnas, Iac, Ibc, Inas, Inbs, weight);
                   }
                 }

               } /* end loop over Jblock */
             } /* end Ivec offdiag */

           } /* end loop over Iblock */
         } /* end loop over Jbuf */
       } /* end loop over Ibuf */
     } /* end loop over roots */
   } /* end icore==2 */

   else {
     throw PSIEXCEPTION("CIWavefunction::tpdm: unrecognized core option!\n");
   }

   /* write and/or print the tpdm
      total and spincases are written out. however, total 2-pdm must be scaled by 1/2
      to conform the stupid psi convention. thus the sum of spincases does not equal
      the total 2-pdm.
    */



   // Build summed TPDM
   SharedVector tpdm_ns(new Vector("MO-basis TPDM (unsymmetrized)", ntri2));
   double* tpdm_nsp = tpdm_ns->pointer();
   orbsym = CalcInfo_->orbsym + ndrc;
   for (i=0; i<CalcInfo_->num_ci_orbs; i++) {
   for (j=0; j<CalcInfo_->num_ci_orbs; j++) {
   for (k=0; k<=i; k++) {

     if (k==i) lmax = j+1;
     else lmax = CalcInfo_->num_ci_orbs;
     ijksym = orbsym[i] ^ orbsym[j] ^ orbsym[k];

     for (l=0; l<lmax; l++) {
       if ((orbsym[l] ^ ijksym) != 0) continue;
       ij = i * CalcInfo_->num_ci_orbs + j;
       kl = k * CalcInfo_->num_ci_orbs + l;
       long int ijkl_tri = INDEX(ij,kl);
       long int ijkl = ij * ntri + kl;
       long int klij = kl * ntri + ij;
       
       double value_aa = tpdm_aap[ijkl_tri];
       double value_bb = tpdm_bbp[ijkl_tri];
       double value_ab = tpdm_abp[ijkl];
       double value_ba = tpdm_abp[klij];
       // For PSI the factor of 1/2 is not pulled outside...so put
       // it back inside now and then write out to the IWL buffer.
       value = 0.5 * (value_aa + value_bb + value_ab + value_ba);
       tpdm_nsp[INDEX(ij, kl)] = value;
     }
   }}}

   //outfile->Printf("Wrote MO-basis TPDM\n");
   //outfile->Printf("Build MO-basis TPDM\n");
   SharedVector tpdm(new Vector("MO-basis TPDM", (ntri * (ntri + 1))/2 ));
   double* tpdmp = tpdm->pointer();

   // Symmetrize and reorder
   for (int p=0, target=0; p<CalcInfo_->num_ci_orbs; p++) {
     for (int q=0; q<=p; q++) {
       for (int r=0; r<=p; r++) {
         int smax = (p == r) ? q+1 : r+1;
         for (int s=0; s<smax; s++) {

          int r_p = CalcInfo_->act_order[p];
          int r_q = CalcInfo_->act_order[q];
          int r_r = CalcInfo_->act_order[r];
          int r_s = CalcInfo_->act_order[s];
          int r_pq = INDEX(r_p, r_q);
          int r_rs = INDEX(r_r, r_s);
          int tpdm_idx = INDEX(r_pq, r_rs);

          int pq = p * CalcInfo_->num_ci_orbs + q;
          int qp = q * CalcInfo_->num_ci_orbs + p;
          int rs = r * CalcInfo_->num_ci_orbs + s;
          int sr = s * CalcInfo_->num_ci_orbs + r;
          int pqrs = INDEX(pq,rs);
          int qprs = INDEX(qp,rs);
          int pqsr = INDEX(pq,sr);
          int qpsr = INDEX(qp,sr); 
          /* would be 0.25 but the formulae I used for the diag hessian
           * seem to define the TPDM with the 1/2 back outside */
          tpdmp[tpdm_idx] = 0.5 * (tpdm_nsp[pqrs] + tpdm_nsp[qprs] +
                         tpdm_nsp[pqsr] + tpdm_nsp[qpsr]);

         }
       }
     }
   }
   tpdm_ns.reset();

   Ivec.buf_unlock();
   Jvec.buf_unlock();
   if (transp_tmp != NULL) free_block(transp_tmp);
   if (transp_tmp2 != NULL) free_block(transp_tmp2);
   free(buffer1);
   free(buffer2);

  std::vector<SharedVector> ret_list;
  ret_list.push_back(tpdm_aa); 
  ret_list.push_back(tpdm_ab); 
  ret_list.push_back(tpdm_bb); 
  ret_list.push_back(tpdm); 
  return ret_list;
}



void CIWavefunction::tpdm_block(struct stringwr **alplist, struct stringwr **betlist,
	        int nbf, int nalplists, int nbetlists,
		double *twopdm_aa, double *twopdm_bb, double *twopdm_ab, double **CJ, double **CI, int Ja_list, 
		int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list, 
		int Inas, int Inbs, double weight)
{
   const int nbf2 = nbf * nbf;
   int Ia_idx, Ib_idx, Ja_idx, Jb_idx, Ja_ex, Jb_ex, Jbcnt, Jacnt; 
   int Kbcnt, Kacnt, Kb_ex, Ka_ex, Kb_list, Ka_list, Kb_idx, Ka_idx;
   struct stringwr *Jb, *Ja, *Kb, *Ka;
   signed char *Jbsgn, *Jasgn, *Kbsgn, *Kasgn;
   unsigned int *Jbridx, *Jaridx, *Kbridx, *Karidx;
   double C1, C2, Ib_sgn, Ia_sgn, Kb_sgn, Ka_sgn, tval;
   int i, j, k, l, ij, kl, ijkl, oij, okl, *Jboij, *Jaoij, *Kboij, *Kaoij;
 
  /* loop over Ia in Ia_list */
  if (Ia_list == Ja_list) {
    for (Ia_idx=0; Ia_idx<Inas; Ia_idx++) {
      for (Jb=betlist[Jb_list], Jb_idx=0; Jb_idx<Jnbs; Jb_idx++, Jb++) {
	C1 = CJ[Ia_idx][Jb_idx] * weight;

	/* loop over excitations E^b_{kl} from |B(J_b)> */
	for (Kb_list=0; Kb_list < nbetlists; Kb_list++) {
	  Jbcnt = Jb->cnt[Kb_list];
	  Jbridx = Jb->ridx[Kb_list];
	  Jbsgn = Jb->sgn[Kb_list];
	  Jboij = Jb->oij[Kb_list];
	  for (Jb_ex=0; Jb_ex < Jbcnt; Jb_ex++) {
	    okl = *Jboij++;
	    Kb_idx = *Jbridx++;
	    Kb_sgn = (double) *Jbsgn++;

            Kb = betlist[Kb_list] + Kb_idx;
	    if (Kb_list == Ib_list) {
	      C2 = CI[Ia_idx][Kb_idx];
	      i = okl / nbf;
	      l = okl % nbf;
	      for (j=0; j<nbf && j<=i; j++) {
		ij = i * nbf + j;
		kl = j * nbf + l;
		if (ij >= kl) {
		  ijkl = INDEX(ij,kl);
		  twopdm_bb[ijkl] -= Kb_sgn * C1 * C2;
		}
	      }
	    }

	    /* loop over excitations E^b_{ij} from |B(K_b)> */
	    /* Ib_list pre-determined because of C blocking */
	    Kbcnt = Kb->cnt[Ib_list];
	    Kbridx = Kb->ridx[Ib_list];
	    Kbsgn = Kb->sgn[Ib_list];
	    Kboij = Kb->oij[Ib_list];
	    for (Kb_ex=0; Kb_ex<Kbcnt; Kb_ex++) {
	      Ib_idx = *Kbridx++;
	      Ib_sgn = (double) *Kbsgn++;
	      oij = *Kboij++;
	      if (oij >= okl) {
		C2 = CI[Ia_idx][Ib_idx];
		ijkl = INDEX(oij,okl);
		twopdm_bb[ijkl] += Ib_sgn * Kb_sgn * C1 * C2;
	      }
	    }

	  } /* end loop over Jb_ex */
	} /* end loop over Kb_list */
      } /* end loop over Jb_idx */
    } /* end loop over Ia_idx */
  } /* end case Ia_list == Ja_list */

  /* loop over Ib in Ib_list */
  if (Ib_list == Jb_list) {
    for (Ib_idx=0; Ib_idx<Inbs; Ib_idx++) {
      for (Ja=alplist[Ja_list], Ja_idx=0; Ja_idx<Jnas; Ja_idx++, Ja++) {
	C1 = CJ[Ja_idx][Ib_idx] * weight;

	/* loop over excitations E^a_{kl} from |A(J_a)> */
	for (Ka_list=0; Ka_list < nalplists; Ka_list++) {
	  Jacnt = Ja->cnt[Ka_list];
	  Jaridx = Ja->ridx[Ka_list];
	  Jasgn = Ja->sgn[Ka_list];
	  Jaoij = Ja->oij[Ka_list];
	  for (Ja_ex=0; Ja_ex < Jacnt; Ja_ex++) {
	    okl = *Jaoij++;
	    Ka_idx = *Jaridx++;
	    Ka_sgn = (double) *Jasgn++;

            Ka = alplist[Ka_list] + Ka_idx;
	    if (Ka_list == Ia_list) {
	      C2 = CI[Ka_idx][Ib_idx];
	      i = okl / nbf;
	      l = okl % nbf;
	      for (j=0; j<nbf && j<=i; j++) {
		ij = i * nbf + j;
		kl = j * nbf + l;
		if (ij >= kl) {
		  ijkl = INDEX(ij,kl);
		  twopdm_aa[ijkl] -= Ka_sgn * C1 * C2;
		}
	      }
	    }

	    /* loop over excitations E^a_{ij} from |A(K_a)> */
	    /* Ia_list pre-determined because of C blocking */
	    Kacnt = Ka->cnt[Ia_list];
	    Karidx = Ka->ridx[Ia_list];
	    Kasgn = Ka->sgn[Ia_list];
	    Kaoij = Ka->oij[Ia_list];
	    for (Ka_ex=0; Ka_ex<Kacnt; Ka_ex++) {
	      Ia_idx = *Karidx++;
	      Ia_sgn = (double) *Kasgn++;
	      oij = *Kaoij++;
	      if (oij >= okl) {
		C2 = CI[Ia_idx][Ib_idx];
		ijkl = INDEX(oij,okl);
		twopdm_aa[ijkl] += Ia_sgn * Ka_sgn * C1 * C2;
	      }
	    }

	  } /* end loop over Ja_ex */
	} /* end loop over Ka_list */
      } /* end loop over Ja_idx */
    } /* end loop over Ib_idx */
  } /* end case Ib_list == Jb_list */


  /* now do the sigma3 looking (alpha-beta) part */
  /* loop over Ja                                */
  for (Ja=alplist[Ja_list], Ja_idx=0; Ja_idx<Jnas; Ja_idx++, Ja++) {

    /* loop over excitations E^a_{kl} from |A(I_a)> */
    Jacnt = Ja->cnt[Ia_list];
    Jaridx = Ja->ridx[Ia_list];
    Jasgn = Ja->sgn[Ia_list];
    Jaoij = Ja->oij[Ia_list];
    for (Ja_ex=0; Ja_ex < Jacnt; Ja_ex++) {
      okl = *Jaoij++;
      Ia_idx = *Jaridx++;
      Ia_sgn = (double) *Jasgn++;

      /* loop over Jb */
      for (Jb=betlist[Jb_list], Jb_idx=0; Jb_idx<Jnbs; Jb_idx++, Jb++) {

	C1 = CJ[Ja_idx][Jb_idx] * weight;

	/* loop over excitations E^b_{ij} from |B(J_b)> */
	Jbcnt = Jb->cnt[Ib_list];
	Jbridx = Jb->ridx[Ib_list];
	Jbsgn = Jb->sgn[Ib_list];
	Jboij = Jb->oij[Ib_list];

	for (Jb_ex=0; Jb_ex < Jbcnt; Jb_ex++) {
	  oij = *Jboij++;
	  Ib_idx = *Jbridx++;
	  Ib_sgn = (double) *Jbsgn++;
	  C2 = CI[Ia_idx][Ib_idx];
	  // alpha-beta matrix is stored without packing bra and ket together
      ijkl = oij * nbf2 + okl;
	  tval = Ib_sgn * Ia_sgn * C1 * C2;
	  // in orbital (i.e. non-spi-orbital) code had to scale the diagonal by 2
	  // because d(ij,kl) += d_ab(ij,kl) + d_ab(kl,ij), hence
	  // d(ij,ij) += 2 d_ab(ij,ij)
	  //if (oij == okl) tval *= 2.0;
	  twopdm_ab[ijkl] += tval;
	}
      } /* end loop over Jb */
    } /* end loop over Ja_ex */
  } /* end loop over Ja */

}

}} // namespace psi::detci
