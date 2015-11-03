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

/*! \defgroup DETCI detci: The Determinant CI code */

/*! \file
    \ingroup DETCI
    \brief Determinant-based CI program

   DETCI

   DETERMINANT CI Program, incorporating Abelian point-group symmetry

   C. David Sherrill
   Center for Computational Quantum Chemistry
   University of Georgia
   August 1994

   Updated 3/95 to do frozen core and virtuals correctly
   Updated 5/95 to do RAS CI's again
   Updated 2/96 to clean up code and rename DETCI

*/

//#include <sstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <cstring>
#include <psifiles.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libqt/slaterdset.h>

#include <libpsi4util/libpsi4util.h>
#include <pthread.h>
#include <iostream>
#include "structs.h"
#include "globals.h"
#include "globaldefs.h"
#include "ci_tol.h"
#include "tpool.h"
#include "odometer.h"
#include "slaterd.h"
#include "civect.h"
#include "ciwave.h"

namespace psi { namespace detci {

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

unsigned char ***Occs;
struct stringwr **alplist;
struct stringwr **betlist;

extern void eivout_t(double **evecs, double *evals, int rows, int cols,
   std::string OutFileRMR);
//extern void read_integrals(void);
//extern void zapt_shift(double *TEI, int nirreps, int nmo, int *doccpi,
//   int *soccpi, int *orbspi, int *frzdoccpi, int *reorder);
extern void print_vec(unsigned int nprint, int *Iacode, int *Ibcode,
   int *Iaidx, int *Ibidx, double *coeff,
   struct olsen_graph *AlphaG, struct olsen_graph *BetaG,
   struct stringwr **alplist, struct stringwr **betlist,
   std::string OutFileRMR);
//extern void print_config(int nbf, int num_alp_el, int num_bet_el,
//   struct stringwr *stralp, struct stringwr *strbet,
//   int num_drc_orbs, char *outstring);
//extern void init_stringwr_temps(int nel, int num_ci_orbs, int nsym);
//extern void free_stringwr_temps(int nsym);
//extern void str_abs2rel(int absidx, int *relidx, int *listnum,
//   struct olsen_graph *Graph);
//extern int str_rel2abs(int relidx, int listnum, struct olsen_graph *Graph);
//extern void H0block_init(unsigned int size);
//extern void H0block_fill(struct stringwr **alplist,
//   struct stringwr **betlist);
//extern void H0block_free(void);
//extern void H0block_print(void);
//extern void H0block_setup(int num_blocks, int *Ia_code, int *Ib_code);
//extern void H0block_pairup(int guess);
//extern void H0block_spin_cpl_chk(void);
//extern void H0block_filter_setup(void);
extern void sem_test(double **A, int N, int M, int L, double **evecs,
   double *evals, double **b, double conv_e, double conv_rms,
   int maxiter, double offst, int *vu, int maxnvect, std::string OutFileRMR);
extern void form_ov(struct stringwr **alplist);
extern void write_energy(int nroots, double *evals, double offset);

//void cleanup(void);
void quote(void);
//void mpn(struct stringwr **strlista, struct stringwr **strlistb);
//void form_opdm(void);
//void form_tpdm(void);
//extern void mitrush_iter(CIvect &Hd,
//   struct stringwr **alplist, struct stringwr **betlist,
//   int nroots, double *evals, double conv_rms, double conv_e, double enuc,
//   double edrc, int maxiter, int maxnvect, std::string OutFileRMR,
//   int print_lvl);
//extern void sem_iter(CIvect &Hd, struct stringwr **alplist, struct stringwr
//   **betlist, double *evals, double conv_e,
//   double conv_rms, double enuc, double edrc,
//   int nroots, int maxiter, int maxnvect, std::string OutFileRMR, int print_lvl);
//extern void mpn_generator(CIvect &Hd, struct stringwr **alplist,
//   struct stringwr **betlist);
//extern void opdm(struct stringwr **alplist, struct stringwr **betlist,
//   int transdens, int dipmom,
//   int Inroots, int Iroot, int Inunits, int Ifirstunit,
//   int Jnroots, int Jroot, int Jnunits, int Jfirstunit,
//   int targetfile, int writeflag, int printflag);
//extern void tpdm(struct stringwr **alplist, struct stringwr **betlist,
//   int Inroots, int Inunits, int Ifirstunit,
//   int Jnroots, int Jnunits, int Jfirstunit,
//   int targetfile, int writeflag, int printflag);
//extern void compute_cc(void);
//extern void calc_mrpt(void);

PsiReturnType detci(Options &options);

}} // namespace psi::detci


namespace psi { namespace detci {


PsiReturnType detci(Options &options)
{

   boost::shared_ptr<Wavefunction> refwfn = Process::environment.wavefunction();
   boost::shared_ptr<CIWavefunction> ciwfn(new CIWavefunction(refwfn, options));

   if (Parameters.nthreads > 1)
     tpool_init(&thread_pool, Parameters.nthreads, CalcInfo.num_alp_str, 0);
                                /* initialize thread pool */
   init_time_new(detci_time);             /* initialize timing routines */

   if (Parameters.istop) {      /* Print size of space, other stuff, only   */
     ciwfn->cleanup();
     Process::environment.globals["CURRENT ENERGY"] = 0.0;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = 0.0;
     Process::environment.globals["CI TOTAL ENERGY"] = 0.0;
     Process::environment.globals["CI CORRELATION ENERGY"] = 0.0;

     return Success;
   }

   // MCSCF is special, we let it handle a lot of its own issues
   if (Parameters.mcscf){
     ciwfn->compute_mcscf();
   }
   else{
     // Transform and set ci integrals
     ciwfn->transform_ci_integrals();

     if (Parameters.mpn){
       ciwfn->compute_mpn();
       }
     else if (Parameters.cc)
       ciwfn->compute_cc();
     else
       ciwfn->diag_h();
   }

   // Finished CI, setting wavefunction parameters
   if(!Parameters.zaptn & Parameters.opdm){
     ciwfn->form_opdm();
     ciwfn->set_opdm();
   }
   else{
     ciwfn->set_opdm(true);
   }

   if (Parameters.tpdm) ciwfn->form_tpdm();
   if (Parameters.print_lvl) print_time_new(detci_time);
   if (Parameters.nthreads > 1) tpool_destroy(thread_pool, 1);
   if (Parameters.print_lvl > 0) quote();


   ciwfn->cleanup();
   Process::environment.set_wavefunction((static_cast<boost::shared_ptr<Wavefunction> > (ciwfn)));
   return Success;
}

/*
** init_ioff(): Set up the ioff array for quick indexing
**
*/
void CIWavefunction::init_ioff(void)
{
   int i;

   /* set offsets for ij-type canonical ordering */
   ioff = (int *) malloc (IOFF_MAX * sizeof(int)) ;
   ioff[0] = 0;
   for (i = 1; i < IOFF_MAX ; i++) {
      ioff[i] = ioff[i-1] + i;
      }
}



/*
** cleanup(): Free any allocated memory that wasn't already freed elsewhere
*/
void CIWavefunction::cleanup(void)
{

}

/*
** title(): Function prints a program identification
*/
void CIWavefunction::title(void)
{
  // if (Parameters.print_lvl) {
   outfile->Printf("\n");
   outfile->Printf("         ---------------------------------------------------------\n");
   outfile->Printf("                                 D E T C I  \n");
   outfile->Printf("\n");
   outfile->Printf("                             C. David Sherrill\n") ;
   outfile->Printf("                             Matt L. Leininger\n") ;
   outfile->Printf("                               18 June 1999\n") ;
   outfile->Printf("         ---------------------------------------------------------\n");
   outfile->Printf("\n");
}







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
// void diag_h(struct stringwr **alplist, struct stringwr **betlist)
{
   BIGINT size;
   int nroots, i, j;
   double conv_rms, conv_e, *evals, **evecs, nucrep, edrc, tval;
   int *tptr;
   double *cbuf;
   char e_label[PSIO_KEYLEN]; /* 80... */

   nroots = Parameters.num_roots;

   conv_rms = Parameters.convergence;
   conv_e = Parameters.energy_convergence;

   if (Parameters.have_special_conv) {
     tval = sqrt(conv_rms) * 10.0;
     if (tval > conv_e) conv_e = tval;
     if (Parameters.special_conv > conv_rms)
       conv_rms = Parameters.special_conv;
   }
   outfile->Printf("About to print vectlen!\n");
   outfile->Printf("Vectlen: %ld", CIblks_->vectlen);
   size = CIblks_->vectlen;
   outfile->Printf("here!\n");
   if ((BIGINT) Parameters.nprint > size) Parameters.nprint = (int) size;
   nucrep = CalcInfo.enuc;
   edrc = CalcInfo.edrc;
   // tmp_ras_array = init_array(1024);

   H0block.size = 0;
   H0block.osize = 0;

   if (Parameters.bendazzoli)
      outfile->Printf( "\nBendazzoli algorithm selected for sigma3\n");

   /* Direct Method --- use RSP diagonalization routine */
   if (Parameters.diag_method == METHOD_RSP) {

      CIvect Cvec(CIblks_->vectlen, CIblks_->num_blocks, 1, Parameters.Ms0,
         CIblks_->Ia_code, CIblks_->Ib_code, CIblks_->Ia_size, CIblks_->Ib_size,
         CIblks_->offset, CIblks_->num_alp_codes, CIblks_->num_bet_codes,
         CalcInfo.nirreps, AlphaG_->subgr_per_irrep, 1, 0, 0,
         CIblks_->first_iablk, CIblks_->last_iablk, CIblks_->decode);
      // shouldn't need to open I/O files for this fake CIvec, unit=0

      double **H, **rsp_evecs;
      int Iarel, Ialist, Ibrel, Iblist;
      BIGINT ii, jj;
      SlaterDeterminant I, J;
      int *mi_iac, *mi_ibc, *mi_iaidx, *mi_ibidx;
      double *mi_coeff;

      if (Parameters.print_lvl) {
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
              I.set(CalcInfo.num_alp_expl,alplist[CIblks_->Ia_code[blk]][ii].occs,
                   CalcInfo.num_bet_expl,betlist[CIblks_->Ib_code[blk]][jj].occs);
              for (ii2=0,det2=0; ii2<CIblks_->Ia_size[blk2]; ii2++) {
                for (jj2=0; jj2<CIblks_->Ib_size[blk2]; jj2++,det2++) {
                  J.set(CalcInfo.num_alp_expl,
                        alplist[CIblks_->Ia_code[blk2]][ii2].occs,
                        CalcInfo.num_bet_expl,
                        betlist[CIblks_->Ib_code[blk2]][jj2].occs);
                  Hpart[det1][det2] = matrix_element(&I,&J);
                }
              }
            }
          }
          if (Parameters.print_lvl > 4 && size < 200) {
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
         I.set(CalcInfo.num_alp_expl,
               alplist_[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
               betlist_[Iblist][Ibrel].occs);
         /* introduce symmetry or other restrictions here */
         for (jj=0; jj<ii; jj++) {
            Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
            J.set(CalcInfo.num_alp_expl,
               alplist_[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
               betlist_[Iblist][Ibrel].occs);
            H[ii][jj] = H[jj][ii] = matrix_element(&I, &J);
            }
         Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
         J.set(CalcInfo.num_alp_expl,
            alplist_[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
            betlist_[Iblist][Ibrel].occs);
         H[jj][jj] = matrix_element(&J, &J) + CalcInfo.edrc;
         }

      if (Parameters.print_lvl > 4 && size < 200) {
         outfile->Printf( "\nHamiltonian matrix:\n");
         print_mat(H, size, size, "outfile");
         }

      sq_rsp(size, size, H, evals, 1, rsp_evecs, 1.0E-14);
      if (Parameters.print_lvl > 4) {
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

      if (Parameters.print_lvl) {
         mi_iac = init_int_array(Parameters.nprint);
         mi_ibc = init_int_array(Parameters.nprint);
         mi_iaidx = init_int_array(Parameters.nprint);
         mi_ibidx = init_int_array(Parameters.nprint);
         mi_coeff = init_array(Parameters.nprint);

         for (i=0; i<nroots; i++) {
            outfile->Printf(
               "\n\n* ROOT %2d CI total energy = %17.13lf\n",
               i+1,evals[i]+nucrep);
            Cvec.setarray(evecs[i], size);
            zero_arr(mi_coeff, Parameters.nprint);
            Cvec.max_abs_vals(Parameters.nprint, mi_iac, mi_ibc, mi_iaidx,
               mi_ibidx, mi_coeff, Parameters.neg_only);
            print_vec(Parameters.nprint, mi_iac, mi_ibc, mi_iaidx, mi_ibidx,
               mi_coeff, AlphaG_, BetaG_, alplist_, betlist_, "outfile");
            }

         free(mi_iac);  free(mi_ibc);
         free(mi_iaidx);  free(mi_ibidx);
         free(mi_coeff);
         }

      if (Parameters.write_energy) write_energy(nroots, evals, nucrep);

      /* Dump the vector to a PSIO file
         Added by Edward valeev (June 2002) */
      if (Parameters.export_ci_vector) {
        StringSet alphastrings, betastrings;
        SlaterDetSet dets;
        SlaterDetVector vec;
        short int *drc_occ;
        unsigned char *newocc;

        if (CalcInfo.num_drc_orbs > 0) {
          drc_occ = (short int *)
            malloc(CalcInfo.num_drc_orbs*sizeof(short int));
          for (int l=0; l<CalcInfo.num_drc_orbs; l++) {
            drc_occ[l] = CalcInfo.order[l]; /* put it in Pitzer order */
          }
        }

        newocc = (unsigned char *)
          malloc(((AlphaG_->num_el > BetaG_->num_el) ?
            AlphaG_->num_el : BetaG_->num_el)*sizeof(unsigned char));

        stringset_init(&alphastrings,AlphaG_->num_str,AlphaG_->num_el,
                       CalcInfo.num_drc_orbs, drc_occ);
        int list_gr = 0;
        int offset = 0;
        for(int irrep=0; irrep<AlphaG_->nirreps; irrep++) {
          for(int gr=0; gr<AlphaG_->subgr_per_irrep; gr++,list_gr++) {
            int nlists_per_gr = AlphaG_->sg[irrep][gr].num_strings;
            for(int l=0; l<nlists_per_gr; l++) {
              /* convert occs to Pitzer order */
              for (int n=0; n<AlphaG_->num_el; n++) {
                newocc[n] = (unsigned char)
                  CalcInfo.order[alplist_[list_gr][l].occs[n] +
                                CalcInfo.num_drc_orbs];
              }
              stringset_add(&alphastrings,l+offset,newocc);
            }
            offset += nlists_per_gr;
          }
        }

        stringset_init(&betastrings,BetaG_->num_str,BetaG_->num_el,
                       CalcInfo.num_drc_orbs, drc_occ);
        list_gr = 0;
        offset = 0;
        for(int irrep=0; irrep<BetaG_->nirreps; irrep++) {
          for(int gr=0; gr<BetaG_->subgr_per_irrep; gr++,list_gr++) {
            int nlists_per_gr = BetaG_->sg[irrep][gr].num_strings;
            for(int l=0; l<nlists_per_gr; l++) {
              /* convert occs to Pitzer order */
              for (int n=0; n<BetaG_->num_el; n++) {
                newocc[n] = (unsigned char)
                  CalcInfo.order[betlist_[list_gr][l].occs[n] +
                                CalcInfo.num_drc_orbs];
              }
              stringset_add(&betastrings,l+offset,newocc);
            }
            offset += nlists_per_gr;
          }
        }
        free(newocc);
        if (CalcInfo.num_drc_orbs > 0)
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
   else if (Parameters.diag_method == METHOD_RSPTEST_OF_SEM) {

      // in-core CIvectors, shouldn't need to open files
      CIvect Cvec(CIblks_->vectlen, CIblks_->num_blocks, 1, Parameters.Ms0,
         CIblks_->Ia_code, CIblks_->Ib_code, CIblks_->Ia_size, CIblks_->Ib_size,
         CIblks_->offset, CIblks_->num_alp_codes, CIblks_->num_bet_codes,
         CalcInfo.nirreps, AlphaG_->subgr_per_irrep, 1, 0, 0,
         CIblks_->first_iablk, CIblks_->last_iablk, CIblks_->decode);
      CIvect Hd(CIblks_->vectlen, CIblks_->num_blocks, 1, Parameters.Ms0,
         CIblks_->Ia_code, CIblks_->Ib_code, CIblks_->Ia_size, CIblks_->Ib_size,
         CIblks_->offset, CIblks_->num_alp_codes, CIblks_->num_bet_codes,
         CalcInfo.nirreps, AlphaG_->subgr_per_irrep, 1, 0, 0,
         CIblks_->first_iablk, CIblks_->last_iablk, CIblks_->decode);

      double **H, **b;
      int Ia, Ib, Iarel, Ialist, Ibrel, Iblist, ij, k, l, tmpi, L;
      unsigned long int ii, jj;
      SlaterDeterminant I, J;
      int *mi_iac, *mi_ibc, *mi_iaidx, *mi_ibidx;
      double *mi_coeff;
      int sm_tridim;
      double *sm_evals, *sm_mat, **sm_evecs, tval;

      if (Parameters.print_lvl) {
         outfile->Printf( "\nFind the roots by the SEM Test Method\n");
         outfile->Printf( "(n.b. this is for debugging purposes only!)\n");
         outfile->Printf( "Energy convergence = %3g\n", conv_e);
         outfile->Printf( "RMS CI vector convergence = %3g\n\n", conv_rms);
         }

      H0block_init(size);

      /* get the diagonal elements of H into an array Hd */

      Hd.diag_mat_els(alplist_, betlist_, CalcInfo.onel_ints,
         CalcInfo.twoel_ints, edrc, CalcInfo.num_alp_expl,
         CalcInfo.num_bet_expl, CalcInfo.num_ci_orbs, Parameters.hd_ave);

      /* get the biggest elements and put in H0block */
      if (H0block.size) {
         Hd.max_abs_vals(H0block.size, H0block.alplist, H0block.betlist,
            H0block.alpidx, H0block.betidx, H0block.H00, Parameters.neg_only);
         }

    /* MLL added this line 5-21-98 */
    /*
      if (Parameters.hd_otf) rclose(Parameters.first_hd_tmp_unit,4);
    */

      H0block_setup(CIblks_->num_blocks, CIblks_->Ia_code, CIblks_->Ib_code);
      if (Parameters.hd_ave) {
        H0block_spin_cpl_chk();
         if (H0block.osize - H0block.size) {
            outfile->Printf("H0block size reduced by %d to %d to ensure"
             "completion of spin-coupling sets\n",
             (H0block.osize - H0block.size), H0block.size);

            }
        }
      if (Parameters.Ms0) {
         H0block_pairup(0);
         if (H0block.osize - H0block.size) {
            outfile->Printf("H0block size reduced by %d to ensure pairing.\n",
               (H0block.osize - H0block.size));

            }
         }

      if (Parameters.print_lvl > 4 && Parameters.hd_otf == FALSE) {
         outfile->Printf( "\nDiagonal elements of the Hamiltonian\n");
         Hd.print("outfile");
         }


      if (H0block.size) {
         H0block_fill();
         }

      if (Parameters.print_lvl > 2 && H0block.size) {
         H0block_print();
         }

      if (Parameters.print_lvl > 3 && H0block.size) {
         outfile->Printf( "\n\nH0 Block:\n");
         print_mat(H0block.H0b, H0block.size, H0block.size, "outfile");
         }

      H = init_matrix(size, size);
      evals = init_array(size) ;
      for (ii=0; ii<size; ii++) {
         Cvec.det2strings(ii, &Ialist, &Iarel, &Iblist, &Ibrel);
         I.set(CalcInfo.num_alp_expl,
               alplist_[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
               betlist_[Iblist][Ibrel].occs);
         /* introduce symmetry or other restrictions here */
         for (jj=0; jj<ii; jj++) {
            Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
            J.set(CalcInfo.num_alp_expl,
               alplist_[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
               betlist_[Iblist][Ibrel].occs);
            H[ii][jj] = H[jj][ii] = matrix_element(&I, &J);
            }
         Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
         J.set(CalcInfo.num_alp_expl,
            alplist_[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
            betlist_[Iblist][Ibrel].occs);
         H[jj][jj] = matrix_element(&J, &J) + CalcInfo.edrc;
         }

      /* obtain a set of L orthonormal trial vectors, L > nroots */
      b = (double **) malloc (Parameters.maxnvect * sizeof(double *)) ;
      for (i=0; i<Parameters.maxnvect; i++) {
         if (i<Parameters.num_init_vecs) b[i] = init_array(size) ;
         else b[i] = NULL ;
         }

      evecs = init_matrix(Parameters.num_roots, size);

      L = H0block.size;
      sm_tridim = L * (L + 1) / 2 ;
      sm_mat = init_array(sm_tridim) ;
      sm_evals = init_array(L) ;
      sm_evecs = init_matrix(L, L) ;
      for (i=0, ij=0; i<L; i++)
         for (j=0; j<=i; j++, ij++)
            sm_mat[ij] = H0block.H0b[i][j] ;
      rsp(L, L, sm_tridim, sm_mat, sm_evals, 1, sm_evecs, 1E-14) ;

      /*
      if (Parameters.precon == PRECON_GEN_DAVIDSON) {
        for (i=0; i<H0block.size; i++) {
           H0block.H0b_eigvals[i] = sm_evals[i];
           for (j=0; j<H0block.size; i++)
              H0block.H0b_diag[i][j] = sm_evecs[i][j];
           }
        }
      */

      /* need to fill out sm_evecs into b (pad w/ 0's) */
      cbuf = *(Cvec.blockptr(0));
      Cvec.buf_unlock();
      for (i=0,k=0; i<L && k < Parameters.num_init_vecs; i++) {

         /* check sm_evecs[i] to see if it has the correct spin symmetry */
         for (j=0,tmpi=0; Parameters.Ms0 && j<L && !tmpi; j++) {
            l = H0block.pair[j];
            if (l == -1) {
               throw PsiException("(diag_h sem_test): unpaired H0block member!",__FILE__,__LINE__);
               }
            tval = sm_evecs[l][i];
            if ((int) Parameters.S%2) tval = -tval;
            if (sm_evecs[j][i] - tval > 1.0E-12) tmpi=1;
            }
         if (tmpi) continue;

         for (j=0; j<L; j++) sm_evals[j] = sm_evecs[j][i];

         Cvec.buf_lock(b[k]);
         Cvec.init_vals(0, L, H0block.alplist, H0block.alpidx,
            H0block.betlist, H0block.betidx, H0block.blknum, sm_evals);
         Cvec.buf_unlock();
         k++;
         }
      Cvec.buf_lock(cbuf);
      free(sm_mat);
      free(sm_evals);
      free_matrix(sm_evecs, L);
      for (i=k; i<Parameters.num_init_vecs; i++) free(b[i]);
      L = k;
      if (L < Parameters.num_roots) {
         throw PsiException("(diag_h sem_test): Ooops! L < num_roots!",__FILE__,__LINE__);
         }
      sem_test(H, size, Parameters.num_roots, L, evecs, evals, b, conv_e,
         conv_rms, Parameters.maxiter, (nucrep+CalcInfo.edrc), &i,
         Parameters.maxnvect, "outfile");

      outfile->Printf( "SEM used %d expansion vectors\n", i);

      if (Parameters.print_lvl > 4) {
         eivout(evecs, evals, size, nroots, "outfile") ;
         }
      free_matrix(H, size);

      if (Parameters.print_lvl) {
         mi_iac = init_int_array(Parameters.nprint);
         mi_ibc = init_int_array(Parameters.nprint);
         mi_iaidx = init_int_array(Parameters.nprint);
         mi_ibidx = init_int_array(Parameters.nprint);
         mi_coeff = init_array(Parameters.nprint);

         for (i=0; i<nroots; i++) {
            outfile->Printf(
               "\n\n* ROOT %2d CI total energy = %17.13lf\n",
               i+1,evals[i]+nucrep);
            Cvec.setarray(evecs[i], size);
            zero_arr(mi_coeff, Parameters.nprint);
            Cvec.max_abs_vals(Parameters.nprint, mi_iac, mi_ibc, mi_iaidx,
               mi_ibidx, mi_coeff, Parameters.neg_only);
            print_vec(Parameters.nprint, mi_iac, mi_ibc, mi_iaidx, mi_ibidx,
               mi_coeff, AlphaG_, BetaG_, alplist_, betlist_, "outfile");
            }
         free(mi_iac);  free(mi_ibc);
         free(mi_iaidx);  free(mi_ibidx);
         free(mi_coeff);
         free_matrix(evecs, Parameters.num_roots);
         }

      if (Parameters.write_energy) write_energy(nroots, evals, nucrep);

      } /* end Test of Davidson/Liu section */


   /*
    * Davidson/Liu Simultaneous Expansion Method OR
    * Mitrushenkov's Olsen-modified Davidson Algorithm
    */

   else {

      /* prepare the H0 block */
      H0block_init(size);

      CIvect Hd(CIblks_->vectlen, CIblks_->num_blocks, Parameters.icore,
         Parameters.Ms0, CIblks_->Ia_code, CIblks_->Ib_code, CIblks_->Ia_size,
         CIblks_->Ib_size, CIblks_->offset, CIblks_->num_alp_codes,
         CIblks_->num_bet_codes, CalcInfo.nirreps, AlphaG_->subgr_per_irrep, 1,
         Parameters.num_hd_tmp_units, Parameters.first_hd_tmp_unit,
         CIblks_->first_iablk, CIblks_->last_iablk, CIblks_->decode);

      bool open_old = false;
      if (Parameters.restart) open_old = true;
      Hd.init_io_files(open_old);

      /* get the diagonal elements of H into an array Hd */
      if (!Parameters.restart || (Parameters.restart && Parameters.hd_otf)) {
         if (Parameters.print_lvl > 1) {
            outfile->Printf( "\nForming diagonal elements of H\n");

           }
         Hd.diag_mat_els(alplist_, betlist_, CalcInfo.onel_ints,
            CalcInfo.twoel_ints, edrc, CalcInfo.num_alp_expl,
            CalcInfo.num_bet_expl, CalcInfo.num_ci_orbs, Parameters.hd_ave);
         }
      else {
         Hd.read(0,0);
         }

      /* get the biggest elements and put in H0block */
      if (H0block.size) {

         if (Parameters.print_lvl > 1) {
            outfile->Printf( "\nForming H0 block\n");

           }

         if (!Parameters.hd_otf)
           Hd.max_abs_vals(H0block.size+H0block.coupling_size,
              H0block.alplist, H0block.betlist, H0block.alpidx,
              H0block.betidx, H0block.H00, Parameters.neg_only);
         }

      //if (Parameters.hd_otf) rclose(Parameters.first_hd_tmp_unit,4);
      if (Parameters.hd_otf) psio_close(Parameters.first_hd_tmp_unit,1);

      H0block_setup(CIblks_->num_blocks, CIblks_->Ia_code, CIblks_->Ib_code);
      if (Parameters.filter_guess) H0block_filter_setup();
      if (Parameters.hd_ave) {
        H0block_spin_cpl_chk();
         if ((H0block.osize - H0block.size) && Parameters.print_lvl > 1) {
            outfile->Printf("H0block size reduced by %d to ensure "
             "completion of spin-coupling sets\n",
             (H0block.osize - H0block.size));
            H0block.osize = H0block.size;
            }
         if ((H0block.oguess_size - H0block.guess_size) &&
             Parameters.print_lvl > 1) {
           outfile->Printf("H0block guess size reduced by %d to ensure "
             "completion of spin-coupling sets\n",
             (H0block.oguess_size - H0block.guess_size));
            H0block.oguess_size = H0block.guess_size;
           }
         if ((H0block.ocoupling_size - H0block.coupling_size) &&
             Parameters.print_lvl > 1) {
           outfile->Printf("H0block coupling size reduced by %d to ensure "
             "completion of spin-coupling sets\n",
             (H0block.ocoupling_size - H0block.coupling_size));
            H0block.ocoupling_size = H0block.coupling_size;
           }

        }
      if (Parameters.Ms0) {
         /* if (H0block.guess_size < H0block.size) */
         H0block_pairup(0); /* pairup h0block size */
         H0block_pairup(1); /* pairup guess_size */
         H0block_pairup(2); /* pairup coupling size */
         if ((H0block.osize - H0block.size) && Parameters.print_lvl > 1) {
           outfile->Printf("H0block size reduced by %d to ensure pairing"
               "and spin-coupling.\n", (H0block.osize - H0block.size));
            }
         if ((H0block.oguess_size - H0block.guess_size) &&
             Parameters.print_lvl > 1) {
           outfile->Printf("H0block guess size reduced by %d to "
               "ensure pairing and spin-coupling.\n",
               (H0block.oguess_size - H0block.guess_size));
            }
         if ((H0block.ocoupling_size - H0block.coupling_size) &&
             Parameters.print_lvl > 1) {
           outfile->Printf("H0block coupling size reduced by %d to "
               "ensure pairing and spin-coupling.\n",
               (H0block.ocoupling_size - H0block.coupling_size));
            }

         }

      Parameters.neg_only = 0; /* MLL 7-2-97 */
      if (Parameters.print_lvl > 4) {
         outfile->Printf( "\nDiagonal elements of the Hamiltonian\n");
         Hd.print("outfile");
         }

      if (H0block.size) {
         H0block_fill();
         }

      if (Parameters.print_lvl > 2 && H0block.size) {
         H0block_print();
         }

      if (Parameters.print_lvl > 3 && H0block.size) {
         outfile->Printf( "\n\nH0 Block:\n");
         print_mat(H0block.H0b, H0block.size, H0block.size, "outfile");
         }

      /* Davidson/Liu Simultaneous Expansion Method */
      if (Parameters.diag_method == METHOD_DAVIDSON_LIU_SEM) {

         if (Parameters.print_lvl) {
            outfile->Printf(
               "\nFind the roots by the Simultaneous Expansion Method ");
            outfile->Printf( "(Block Davidson Method)\n");
            outfile->Printf( "Energy convergence = %3g\n", conv_e);
            outfile->Printf( "RMS CI vector convergence = %3g\n\n", conv_rms);

            }

         evals = init_array(nroots);

         sem_iter(Hd, alplist_, betlist_, evals, conv_e, conv_rms,
            nucrep, edrc, nroots, Parameters.maxiter,
            Parameters.maxnvect, "outfile", Parameters.print_lvl);
         }

      /* Mitrushenkov's Olsen Method */
      else {

         if (Parameters.print_lvl) {
            if (Parameters.diag_method == METHOD_MITRUSHENKOV)
              outfile->Printf(
                "\nFind the roots with Mitrushenkov's two vector algorithm\n");
            else if (Parameters.diag_method == METHOD_OLSEN)
              outfile->Printf(
                "\nFind the roots with Olsen's single vector algorithm\n");
            outfile->Printf( "Energy convergence = %3g\n", conv_e);
            outfile->Printf( "RMS CI vector convergence = %3g\n", conv_rms);

            }

         evals = init_array(nroots);

         mitrush_iter(Hd, alplist_, betlist_, nroots, evals, conv_rms, conv_e,
            nucrep, edrc, Parameters.maxiter, Parameters.maxnvect, "outfile",
            Parameters.print_lvl);

         H0block_free();
         }

      if (Parameters.write_energy) write_energy(nroots, evals, nucrep+edrc);

      } /* end the Davidson-Liu/Mitrushenkov-Olsen-Davidson section */

   /* write the CI energy to PSIF_CHKPT: later fix this to loop over roots */
   chkpt_init(PSIO_OPEN_OLD);
   tval = evals[Parameters.root]+edrc+nucrep;
   chkpt_wt_etot(tval);

   Process::environment.globals["CURRENT ENERGY"] = tval;
   Process::environment.globals["CURRENT CORRELATION ENERGY"] = tval - CalcInfo.escf;
   Process::environment.globals["CURRENT REFERENCE ENERGY"] = CalcInfo.escf;
   Process::environment.globals["CI TOTAL ENERGY"] = tval;
   // eref seems wrong for open shells so replace it with escf below
   // until I fix it ---CDS 11/5/11
   Process::environment.globals["CI CORRELATION ENERGY"] = tval - CalcInfo.escf;

   if (Parameters.fci) {
     Process::environment.globals["FCI TOTAL ENERGY"] = tval;
     Process::environment.globals["FCI CORRELATION ENERGY"] = tval - CalcInfo.escf;
   }
   else {
     if (Parameters.ex_lvl == 2) {
       Process::environment.globals["CISD TOTAL ENERGY"] = tval;
       Process::environment.globals["CISD CORRELATION ENERGY"] = tval - CalcInfo.escf;
     }
     else if (Parameters.ex_lvl == 3) {
       Process::environment.globals["CISDT TOTAL ENERGY"] = tval;
       Process::environment.globals["CISDT CORRELATION ENERGY"] = tval - CalcInfo.escf;
     }
     else if (Parameters.ex_lvl == 4) {
       Process::environment.globals["CISDTQ TOTAL ENERGY"] = tval;
       Process::environment.globals["CISDTQ CORRELATION ENERGY"] = tval - CalcInfo.escf;
     }
     else {
       /*- Process::environment.globals["CIn TOTAL ENERGY"] -*/
       /*- Process::environment.globals["CIn CORRELATION ENERGY"] -*/
       std::stringstream s;
       s << "CI" << Parameters.ex_lvl << " TOTAL ENERGY";
       Process::environment.globals[s.str()] = tval;
       s.str(std::string());
       s << "CI" << Parameters.ex_lvl << " CORRELATION ENERGY";
       Process::environment.globals[s.str()] = tval - CalcInfo.escf;
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
     Process::environment.globals[s.str()] = tval - CalcInfo.escf;
   }

   if (Parameters.average_num > 1) {
     tval = 0.0;
     for (i=0; i<Parameters.average_num; i++)
       tval += Parameters.average_weights[i] *
               (edrc+nucrep+evals[Parameters.average_states[i]]);
     chkpt_wt_e_labeled("State averaged energy",tval);
     Process::environment.globals["CI STATE-AVERAGED TOTAL ENERGY"] = tval;
     // eref seems wrong for open shells so replace it with escf below
     // until I fix it ---CDS 11/5/11
     Process::environment.globals["CI STATE-AVERAGED CORRELATION ENERGY"] =
       tval - CalcInfo.escf;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] =
       Process::environment.globals["CI STATE-AVERAGED CORRELATION ENERGY"];
   }

   // Set the energy as MCSCF would find it
   if (Parameters.average_num > 1) { // state average
     Process::environment.globals["MCSCF TOTAL ENERGY"] =
       Process::environment.globals["CI STATE-AVERAGED TOTAL ENERGY"];
   }
   else if (Parameters.root != 0) { // follow some specific root != lowest
     std::stringstream s;
     s << "CI ROOT " << (Parameters.root+1) << " TOTAL ENERGY";
     Process::environment.globals["MCSCF TOTAL ENERGY"] =
      Process::environment.globals[s.str()];
   }
   else {
     Process::environment.globals["MCSCF TOTAL ENERGY"] =
      Process::environment.globals["CI TOTAL ENERGY"];
   }

   chkpt_close();

}


void quote(void)
{
   outfile->Printf("\t\t \"A good bug is a dead bug\" \n\n");
   outfile->Printf("\t\t\t - Starship Troopers\n\n");
   outfile->Printf("\t\t \"I didn't write FORTRAN.  That's the problem.\"\n\n");
   outfile->Printf("\t\t\t - Edward Valeev\n\n");

}
//BIGINT strings2det(int alp_code, int alp_idx, int bet_code, int bet_idx) {
//
//   int blknum;
//   BIGINT addr;
//
//   blknum = CIblks_->decode[alp_code][bet_code];
//   addr = CIblks_->offset[blknum];
//   addr += alp_idx * CIblks_->Ib_size[blknum] + bet_idx;
//
//   return(addr);
//
//}




}} // namespace psi::detci


