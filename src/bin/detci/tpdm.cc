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
#include "structs.h"
#define EXTERN
#include "globals.h"
#include "civect.h"

namespace psi { namespace detci {

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define TPDMINDEX2(i,j,n) ((i)*(n) + (j))


void tpdm_block(struct stringwr **alplist, struct stringwr **betlist,
		int nbf, int nalpcodes, int nbetcodes,
		double *twopdm_aa, double *twopdm_bb, double *twopdm_ab, double **CJ, double **CI, int Ja_list, 
		int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list, 
		int Inas, int Inbs, double weight);


void tpdm(struct stringwr **alplist, struct stringwr **betlist, 
          int Inroots, int Inunits, int Ifirstunit, 
	  int Jnroots, int Jnunits, int Jfirstunit, 
	  int targetfile, int writeflag, int printflag)
{

   CIvect Ivec, Jvec;
   struct iwlbuf TBuff;
   struct iwlbuf TBuff_aa;
   struct iwlbuf TBuff_bb;
   struct iwlbuf TBuff_ab;
   int i, j, k, l, lmax, ij, kl, ijkl, ijksym;
   int i2, j2, k2, l2, nfzc, populated_orbs;
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

   nfzc = CalcInfo.num_fzc_orbs;
   populated_orbs = CalcInfo.nmo - CalcInfo.num_fzv_orbs;

   if (nfzc) {
     psio_open(Parameters.opdm_file, PSIO_OPEN_OLD);
     onepdm_a = block_matrix(populated_orbs, populated_orbs);
     onepdm_b = block_matrix(populated_orbs, populated_orbs);
   }

   Ivec.set(CIblks.vectlen, CIblks.num_blocks, Parameters.icore, 1,
	     CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
	     CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
	     CalcInfo.nirreps, AlphaG->subgr_per_irrep, Inroots, Inunits,
	     Ifirstunit, CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);

   Jvec.set(CIblks.vectlen, CIblks.num_blocks, Parameters.icore, 1,
	     CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
	     CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
	     CalcInfo.nirreps, AlphaG->subgr_per_irrep, Jnroots, Jnunits,
	     Jfirstunit, CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);

   buffer1 = Ivec.buf_malloc();
   buffer2 = Jvec.buf_malloc();
   Ivec.buf_lock(buffer1);
   Jvec.buf_lock(buffer2);

   ntri = CalcInfo.num_ci_orbs * CalcInfo.num_ci_orbs;
   ntri2 = (ntri * (ntri + 1)) / 2;
   twopdm_aa = init_array(ntri2);
   twopdm_bb = init_array(ntri2);
   twopdm_ab = init_array(ntri * ntri);

   if ((Ivec.icore==2 && Ivec.Ms0 && CalcInfo.ref_sym != 0) || 
       (Ivec.icore==0 && Ivec.Ms0)) {
     for (i=0, maxrows=0, maxcols=0; i<Ivec.num_blocks; i++) {
       if (Ivec.Ia_size[i] > maxrows) maxrows = Ivec.Ia_size[i];
       if (Ivec.Ib_size[i] > maxcols) maxcols = Ivec.Ib_size[i];
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


   if (Parameters.icore == 0) {

     /* loop over all the roots requested */
     for (root_idx=0; root_idx<Parameters.average_num; root_idx++)  {
       Iroot = Parameters.average_states[root_idx];
       weight = Parameters.average_weights[root_idx];
       Jroot = Iroot;  /* change later if need transition matrix elements */

       for (Ibuf=0; Ibuf<Ivec.buf_per_vect; Ibuf++) {
         Ivec.read(Iroot, Ibuf);
         Iblock = Ivec.buf2blk[Ibuf];
         Iac = Ivec.Ia_code[Iblock];
         Ibc = Ivec.Ib_code[Iblock];
         Inas = Ivec.Ia_size[Iblock];
         Inbs = Ivec.Ib_size[Iblock];
       
         for (Jbuf=0; Jbuf<Jvec.buf_per_vect; Jbuf++) {
           do_Jblock=0; do_Jblock2=0;
           Jblock = Jvec.buf2blk[Jbuf];
           Jblock2 = -1;
           Jac = Jvec.Ia_code[Jblock];
           Jbc = Jvec.Ib_code[Jblock];
           if (Jvec.Ms0) Jblock2 = Jvec.decode[Jbc][Jac];
             Jnas = Jvec.Ia_size[Jblock];
             Jnbs = Jvec.Ib_size[Jblock];
             if (s1_contrib[Iblock][Jblock] || s2_contrib[Iblock][Jblock]
               || s3_contrib[Iblock][Jblock]) 
             do_Jblock = 1;
           if (Jvec.buf_offdiag[Jbuf] && (s1_contrib[Iblock][Jblock2] ||
             s2_contrib[Iblock][Jblock2] ||
             s3_contrib[Iblock][Jblock2]))
             do_Jblock2 = 1;
           if (!do_Jblock && !do_Jblock2) continue;

           Jvec.read(Jroot, Jbuf);
	 
           if (do_Jblock) {
             tpdm_block(alplist, betlist, CalcInfo.num_ci_orbs, 
               Ivec.num_alpcodes, Ivec.num_betcodes, twopdm_aa, twopdm_bb, twopdm_ab, 
               Jvec.blocks[Jblock], Ivec.blocks[Iblock], 
               Jac, Jbc, Jnas, Jnbs, Iac, Ibc, Inas, Inbs, weight);
           }
	 
           if (do_Jblock2) {
             Jvec.transp_block(Jblock, transp_tmp);
             tpdm_block(alplist, betlist, CalcInfo.num_ci_orbs,
               Ivec.num_alpcodes, Ivec.num_betcodes, 
               twopdm_aa, twopdm_bb, twopdm_ab, transp_tmp, Ivec.blocks[Iblock], 
               Jbc, Jac, Jnbs, Jnas, Iac, Ibc, Inas, Inbs, weight);
           }
	 
         } /* end loop over Jbuf */
       
         if (Ivec.buf_offdiag[Ibuf]) { /* need to get contrib of transpose */
           Iblock2 = Ivec.decode[Ibc][Iac];
           Iac = Ivec.Ia_code[Iblock2];
           Ibc = Ivec.Ib_code[Iblock2];
           Inas = Ivec.Ia_size[Iblock2];
           Inbs = Ivec.Ib_size[Iblock2];
       
           Ivec.transp_block(Iblock, transp_tmp2);

           for (Jbuf=0; Jbuf<Jvec.buf_per_vect; Jbuf++) {
	     do_Jblock=0; do_Jblock2=0;
	     Jblock = Jvec.buf2blk[Jbuf];
	     Jblock2 = -1;
	     Jac = Jvec.Ia_code[Jblock];
	     Jbc = Jvec.Ib_code[Jblock];
	     if (Jvec.Ms0) Jblock2 = Jvec.decode[Jbc][Jac];
	     Jnas = Jvec.Ia_size[Jblock];
	     Jnbs = Jvec.Ib_size[Jblock];
	     if (s1_contrib[Iblock2][Jblock] || s2_contrib[Iblock2][Jblock] ||
	         s3_contrib[Iblock2][Jblock]) 
	       do_Jblock = 1;
	     if (Jvec.buf_offdiag[Jbuf] && (s1_contrib[Iblock2][Jblock2] ||
               s2_contrib[Iblock2][Jblock2] ||
               s3_contrib[Iblock2][Jblock2]))
               do_Jblock2 = 1;
             if (!do_Jblock && !do_Jblock2) continue;
	   
             Jvec.read(Jroot, Jbuf);
	 
             if (do_Jblock) {
               tpdm_block(alplist, betlist, CalcInfo.num_ci_orbs, 
                 Ivec.num_alpcodes, Ivec.num_betcodes, 
                 twopdm_aa, twopdm_bb, twopdm_ab, Jvec.blocks[Jblock], 
                 transp_tmp2, Jac, Jbc, Jnas,
                 Jnbs, Iac, Ibc, Inas, Inbs, weight);
             }
	   
             if (do_Jblock2) {
               Jvec.transp_block(Jblock, transp_tmp);
               tpdm_block(alplist, betlist, CalcInfo.num_ci_orbs,
                 Ivec.num_alpcodes, Ivec.num_betcodes,
                 twopdm_aa, twopdm_bb, twopdm_ab, transp_tmp, transp_tmp2, 
                 Jbc, Jac, Jnbs, Jnas, Iac, Ibc, Inas, Inbs, weight);
             }
           } /* end loop over Jbuf */

         } /* end loop over Ibuf transpose */
       } /* end loop over Ibuf */
     } /* end loop over roots */
   } /* end icore==0 */

   else if (Parameters.icore==1) { /* whole vectors in-core */

     for (root_idx=0; root_idx<Parameters.average_num; root_idx++)  {
       Iroot = Parameters.average_states[root_idx];
       weight = Parameters.average_weights[root_idx];
       Jroot = Iroot;  /* change later if need transition matrix elements */

       Ivec.read(Iroot, 0);
       Jvec.read(Jroot, 0);
       for (Iblock=0; Iblock<Ivec.num_blocks; Iblock++) {
         Iac = Ivec.Ia_code[Iblock];
         Ibc = Ivec.Ib_code[Iblock];
         Inas = Ivec.Ia_size[Iblock];
         Inbs = Ivec.Ib_size[Iblock];
         if (Inas==0 || Inbs==0) continue;
         for (Jblock=0; Jblock<Jvec.num_blocks; Jblock++) {
           Jac = Jvec.Ia_code[Jblock];
           Jbc = Jvec.Ib_code[Jblock];
           Jnas = Jvec.Ia_size[Jblock];
           Jnbs = Jvec.Ib_size[Jblock];
           if (s1_contrib[Iblock][Jblock] || s2_contrib[Iblock][Jblock] ||
             s3_contrib[Iblock][Jblock])
             tpdm_block(alplist, betlist, CalcInfo.num_ci_orbs,
               Ivec.num_alpcodes, Ivec.num_betcodes, 
               twopdm_aa, twopdm_bb, twopdm_ab, Jvec.blocks[Jblock], Ivec.blocks[Iblock], 
               Jac, Jbc, Jnas, Jnbs, Iac, Ibc, Inas, Inbs, weight);
         }
       } /* end loop over Iblock */
     } /* end loop over roots */
   } /* end icore==1 */

   else if (Parameters.icore==2) { /* icore==2 */
     for (root_idx=0; root_idx<Parameters.average_num; root_idx++)  {
       Iroot = Parameters.average_states[root_idx];
       weight = Parameters.average_weights[root_idx];
       Jroot = Iroot;  /* change later if need transition matrix elements */

       for (Ibuf=0; Ibuf<Ivec.buf_per_vect; Ibuf++) {
         Ivec.read(Iroot, Ibuf);
         Iairr = Ivec.buf2blk[Ibuf];

         for (Jbuf=0; Jbuf<Jvec.buf_per_vect; Jbuf++) {
           Jvec.read(Jroot, Jbuf);
           Jairr = Jvec.buf2blk[Jbuf];

           for (Iblock=Ivec.first_ablk[Iairr]; Iblock<=Ivec.last_ablk[Iairr];
             Iblock++) {
             Iac = Ivec.Ia_code[Iblock];
             Ibc = Ivec.Ib_code[Iblock];
             Inas = Ivec.Ia_size[Iblock];
             Inbs = Ivec.Ib_size[Iblock];
   
             for (Jblock=Jvec.first_ablk[Jairr]; Jblock<=Jvec.last_ablk[Jairr];
               Jblock++) {
               Jac = Jvec.Ia_code[Jblock];
               Jbc = Jvec.Ib_code[Jblock];
               Jnas = Jvec.Ia_size[Jblock];
               Jnbs = Jvec.Ib_size[Jblock];
   
               if (s1_contrib[Iblock][Jblock] || s2_contrib[Iblock][Jblock] ||
                 s3_contrib[Iblock][Jblock])
               tpdm_block(alplist, betlist, CalcInfo.num_ci_orbs, 
                 Ivec.num_alpcodes, Ivec.num_betcodes,
                 twopdm_aa, twopdm_bb, twopdm_ab, Jvec.blocks[Jblock], Ivec.blocks[Iblock], 
                 Jac, Jbc, Jnas, Jnbs, Iac, Ibc, Inas, Inbs, weight);

               if (Jvec.buf_offdiag[Jbuf]) {
                 Jblock2 = Jvec.decode[Jbc][Jac];
                 if (s1_contrib[Iblock][Jblock2] ||
                   s2_contrib[Iblock][Jblock2] ||
                   s3_contrib[Iblock][Jblock2]) {
                   Jvec.transp_block(Jblock, transp_tmp);
                   tpdm_block(alplist, betlist, CalcInfo.num_ci_orbs,
                     Ivec.num_alpcodes, Ivec.num_betcodes,
                     twopdm_aa, twopdm_bb, twopdm_ab, transp_tmp, Ivec.blocks[Iblock], 
                     Jbc, Jac, Jnbs, Jnas, Iac, Ibc, Inas, Inbs, weight);
                 }
               }

             } /* end loop over Jblock */

             if (Ivec.buf_offdiag[Ibuf]) {
               Iblock2 = Ivec.decode[Ibc][Iac];
               Ivec.transp_block(Iblock, transp_tmp2);
               Iac = Ivec.Ia_code[Iblock2];
               Ibc = Ivec.Ib_code[Iblock2];
               Inas = Ivec.Ia_size[Iblock2];
               Inbs = Ivec.Ib_size[Iblock2];
	   
               for (Jblock=Jvec.first_ablk[Jairr]; 
                 Jblock<=Jvec.last_ablk[Jairr]; Jblock++) {
                 Jac = Jvec.Ia_code[Jblock];
                 Jbc = Jvec.Ib_code[Jblock];
                 Jnas = Jvec.Ia_size[Jblock];
                 Jnbs = Jvec.Ib_size[Jblock];
	   
                 if (s1_contrib[Iblock2][Jblock] || s2_contrib[Iblock2][Jblock]
                   || s3_contrib[Iblock2][Jblock])
                   tpdm_block(alplist, betlist, CalcInfo.num_ci_orbs,
                     Ivec.num_alpcodes, Ivec.num_betcodes, 
                     twopdm_aa, twopdm_bb, twopdm_ab, Jvec.blocks[Jblock], transp_tmp2, 
                     Jac, Jbc, Jnas, Jnbs, Iac, Ibc, Inas, Inbs, weight);

                 if (Jvec.buf_offdiag[Jbuf]) {
                   Jblock2 = Jvec.decode[Jbc][Jac];
                   if (s1_contrib[Iblock][Jblock2] ||
                     s2_contrib[Iblock][Jblock2] ||
                     s3_contrib[Iblock][Jblock2]) {
                     Jvec.transp_block(Jblock, transp_tmp);
                     tpdm_block(alplist, betlist, CalcInfo.num_ci_orbs,
                     Ivec.num_alpcodes, Ivec.num_betcodes,
                     twopdm_aa, twopdm_bb, twopdm_ab, transp_tmp, transp_tmp2, Jbc, Jac,
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
     printf("tpdm: unrecognized core option!\n");
     return;
   }

   /* write and/or print the tpdm
      total and spincases are written out. however, total 2-pdm must be scaled by 1/2
      to conform the stupid psi convention. thus the sum of spincases does not equal
      the total 2-pdm.
    */

   if (writeflag) {
     
     if (printflag) fprintf(outfile, "\nTwo-particle density matrix\n\n");
     iwl_buf_init(&TBuff, targetfile, 0.0, 0, 0);
     iwl_buf_init(&TBuff_aa, PSIF_MO_AA_TPDM, 0.0, 0, 0);
     iwl_buf_init(&TBuff_bb, PSIF_MO_BB_TPDM, 0.0, 0, 0);
     iwl_buf_init(&TBuff_ab, PSIF_MO_AB_TPDM, 0.0, 0, 0);
     orbsym = CalcInfo.orbsym + nfzc;

     /* do the core-core and core-active part here */
     if (nfzc) {
         if (Parameters.average_num == 1) 
           sprintf(opdm_key, "MO-basis Alpha OPDM Root %d",
             Parameters.average_states[0]);
         else
           sprintf(opdm_key, "MO-basis Alpha OPDM");

         psio_read_entry(Parameters.opdm_file, opdm_key, (char *) onepdm_a[0],
           populated_orbs * populated_orbs * sizeof(double));

         if (Parameters.average_num == 1) 
           sprintf(opdm_key, "MO-basis Beta OPDM Root %d",
             Parameters.average_states[0]);
         else
           sprintf(opdm_key, "MO-basis Beta OPDM");

         psio_read_entry(Parameters.opdm_file, opdm_key, (char *) onepdm_b[0],
           populated_orbs * populated_orbs * sizeof(double));

       /* core-core part */
       for (i=0; i<nfzc; i++) {
         for (j=0; j<i; j++) {
           iwl_buf_wrt_val(&TBuff,i,i,j,j, 2.00,printflag,outfile,0);
           iwl_buf_wrt_val(&TBuff_aa,i,i,j,j, 1.00,0,outfile,0);
           iwl_buf_wrt_val(&TBuff_bb,i,i,j,j, 1.00,0,outfile,0);
           iwl_buf_wrt_val(&TBuff_ab,i,i,j,j, 1.00,0,outfile,0);
           iwl_buf_wrt_val(&TBuff_ab,j,j,i,i, 1.00,0,outfile,0);

           iwl_buf_wrt_val(&TBuff,i,j,j,i,-1.00,printflag,outfile,0);
           iwl_buf_wrt_val(&TBuff_aa,i,j,j,i,-1.00,0,outfile,0);
           iwl_buf_wrt_val(&TBuff_bb,i,j,j,i,-1.00,0,outfile,0);
         }
         iwl_buf_wrt_val(&TBuff,i,i,i,i,1.0,printflag,outfile,0);
         iwl_buf_wrt_val(&TBuff_ab,i,i,i,i,1.0,0,outfile,0);
       }

       /* core-active part */
       for (i=nfzc; i<populated_orbs; i++) {
         for (j=nfzc; j<populated_orbs; j++) {
           const double value_a = onepdm_a[i][j];
           const double value_b = onepdm_b[i][j];
           for (k=0; k<nfzc; k++) {
             iwl_buf_wrt_val(&TBuff,i,j,k,k,value_a + value_b,printflag,outfile,0);
             iwl_buf_wrt_val(&TBuff_aa,i,j,k,k,value_a,0,outfile,0);
             iwl_buf_wrt_val(&TBuff_bb,i,j,k,k,value_b,0,outfile,0);
             iwl_buf_wrt_val(&TBuff_ab,i,j,k,k,value_a,0,outfile,0);
             iwl_buf_wrt_val(&TBuff_ab,k,k,i,j,value_b,0,outfile,0);

             iwl_buf_wrt_val(&TBuff,i,k,k,j,-0.5*(value_a+value_b),printflag,outfile,0);
             iwl_buf_wrt_val(&TBuff_aa,i,k,k,j,-value_a,0,outfile,0);
             iwl_buf_wrt_val(&TBuff_bb,i,k,k,j,-value_b,0,outfile,0);
           }
         }
       }
     }
 
     for (i=0; i<CalcInfo.num_ci_orbs; i++) {
       i2 = i+ nfzc;
       for (j=0; j<CalcInfo.num_ci_orbs; j++) {
         j2 = j + nfzc;
	 for (k=0; k<=i; k++) {
           k2 = k + nfzc;
	   if (k==i) lmax = j+1;
	   else lmax = CalcInfo.num_ci_orbs;
	   ijksym = orbsym[i] ^ orbsym[j] ^ orbsym[k];
	   for (l=0; l<lmax; l++) {
             l2 = l + nfzc;
	     if ((orbsym[l] ^ ijksym) != 0) continue;
	     ij = i * CalcInfo.num_ci_orbs + j;
	     kl = k * CalcInfo.num_ci_orbs + l;
	     const long int ijkl_tri = INDEX(ij,kl);
         const long int ijkl = TPDMINDEX2(ij,kl,ntri);
         const long int klij = TPDMINDEX2(kl,ij,ntri);
	     
         const double value_aa = twopdm_aa[ijkl_tri];
         const double value_bb = twopdm_bb[ijkl_tri];
         const double value_ab = twopdm_ab[ijkl];
         const double value_ba = twopdm_ab[klij];
             // For PSI the factor of 1/2 is not pulled outside...so put
             // it back inside now and then write out to the IWL buffer.
             value = 0.5 * (value_aa + value_bb + value_ab + value_ba);
	     iwl_buf_wrt_val(&TBuff,i2,j2,k2,l2,value,printflag,outfile,0);
         iwl_buf_wrt_val(&TBuff_aa,i2,j2,k2,l2,value_aa,0,outfile,0);
         iwl_buf_wrt_val(&TBuff_bb,i2,j2,k2,l2,value_bb,0,outfile,0);
         iwl_buf_wrt_val(&TBuff_ab,i2,j2,k2,l2,value_ab,0,outfile,0);
         if (ij != kl)
           iwl_buf_wrt_val(&TBuff_ab,k2,l2,i2,j2,value_ba,0,outfile,0);
	   }
	 }
       }
     }
     iwl_buf_flush(&TBuff, 1);
     iwl_buf_close(&TBuff, 1);
     iwl_buf_flush(&TBuff_aa, 1);
     iwl_buf_close(&TBuff_aa, 1);
     iwl_buf_flush(&TBuff_bb, 1);
     iwl_buf_close(&TBuff_bb, 1);
     iwl_buf_flush(&TBuff_ab, 1);
     iwl_buf_close(&TBuff_ab, 1);
     fprintf(outfile, "\n");
   }


   if (nfzc) {
     psio_close(Parameters.opdm_file, 1);
     free_block(onepdm_a);
     free_block(onepdm_b);
   }

   Ivec.buf_unlock();
   Jvec.buf_unlock();
   free(twopdm_aa);
   free(twopdm_bb);
   free(twopdm_ab);
   if (transp_tmp != NULL) free_block(transp_tmp);
   if (transp_tmp2 != NULL) free_block(transp_tmp2);
   free(buffer1);
   free(buffer2);
}



void tpdm_block(struct stringwr **alplist, struct stringwr **betlist,
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
	  ijkl = TPDMINDEX2(oij, okl, nbf2);
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

