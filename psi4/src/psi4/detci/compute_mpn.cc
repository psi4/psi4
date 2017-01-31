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

/*
** Moeller Plesset Perturbation Series Generator
**
** Matt L. Leininger
** August 11, 1998
**
*/

/* #define DEBUG */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/vector.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/ci_tol.h"
#include "psi4/detci/civect.h"
#include "psi4/detci/ciwave.h"

namespace psi { namespace detci {

extern void print_vec(unsigned int nprint, int *Iacode, int *Ibcode,
   int *Iaidx, int *Ibidx, double *coeff,
   struct olsen_graph *AlphaG, struct olsen_graph *BetaG,
   struct stringwr **alplist, struct stringwr **betlist,
   std::string );


/*
** mpn(): Function which sets up and generates the mpn series
**
** Returns: none
*/
void CIWavefunction::compute_mpn()
{
  int i, j, irrep, cnt;
  struct stringwr *stralp, *strbet;
  int **drc_orbs;
  double tval;

  if (print_){
      outfile->Printf("\n   ==> Starting MPn CI Computation <==\n\n");
  }

  if(Parameters_->zaptn){         /* Shift SCF eigenvalues for ZAPTn          */
    int h1, h2;
    int x, y, i, j;
    int offset, offset2;
    int ij, ijij;
    int docc, socc;
    int docc2, socc2;
    int totfzc;

    for(h1=0,totfzc=0; h1<CalcInfo_->nirreps; h1++)
        totfzc += CalcInfo_->dropped_docc[h1];

    for(h1 = 0,offset = 0; h1 < CalcInfo_->nirreps; h1++) {
        if(h1>0) offset += nmopi_[h1-1];
        docc = CalcInfo_->docc[h1];
        socc = CalcInfo_->socc[h1];
        for(x = offset+docc; x<offset+docc+socc; x++){
            for(h2 = 0,offset2 = 0; h2 < CalcInfo_->nirreps; h2++) {
                if(h2>0) offset2 += nmopi_[h2-1];
                docc2 = CalcInfo_->docc[h2];
                socc2 = CalcInfo_->socc[h2];
                for(y=offset2+docc2;y<offset2+docc2+socc2;y++) {
                    i = CalcInfo_->reorder[x] - totfzc;
                    j = CalcInfo_->reorder[y] - totfzc;
                    ij = INDEX2(i,j);
                    ijij = ioff[ij] + ij;
                    CalcInfo_->scfeigvala[i+totfzc] -= 0.5*CalcInfo_->twoel_ints->pointer()[ijij];
                    CalcInfo_->scfeigvalb[i+totfzc] += 0.5*CalcInfo_->twoel_ints->pointer()[ijij];
                }
            }
        }
    }

  }

  CIvect Hd(Parameters_->icore, 1, 1,
         Parameters_->hd_filenum, CIblks_, CalcInfo_, Parameters_, H0block_);

  Hd.init_io_files(false);

  /* Compute E0 from orbital energies */
  stralp = alplist_[CalcInfo_->ref_alp_list] + CalcInfo_->ref_alp_rel;
  strbet = betlist_[CalcInfo_->ref_bet_list] + CalcInfo_->ref_bet_rel;

  drc_orbs = init_int_matrix(CalcInfo_->nirreps, CalcInfo_->num_drc_orbs);
  cnt = 0;
  for (irrep=0; irrep<CalcInfo_->nirreps; irrep++)
     for (i=0; i<CalcInfo_->dropped_docc[irrep]; i++)
        drc_orbs[irrep][i] = cnt++;

  /* Loop over alp occs */
  //CalcInfo_->e0 = CalcInfo_->edrc;
  CalcInfo_->e0 = 0.0;
  CalcInfo_->e0_drc = 0.0;
  for (i=0; i<CalcInfo_->num_drc_orbs; i++) {
     outfile->Printf(" orb_energy[%d] = %lf\n", i, CalcInfo_->scfeigval[i]);
     tval = 2.0 * CalcInfo_->scfeigval[i];
     CalcInfo_->e0 += tval;
     CalcInfo_->e0_drc += tval;
     }

  if(Parameters_->zaptn) {
    for (i=0; i<CalcInfo_->num_alp_expl; i++) {
       j = (stralp->occs)[i] + CalcInfo_->num_drc_orbs;
       CalcInfo_->e0 += CalcInfo_->scfeigvala[j];
       }

    for (i=0; i<CalcInfo_->num_bet_expl; i++) {
       j = (strbet->occs)[i] + CalcInfo_->num_drc_orbs;
       CalcInfo_->e0 += CalcInfo_->scfeigvalb[j];
       }
    } else {
    for (i=0; i<CalcInfo_->num_alp_expl; i++) {
       j = (stralp->occs)[i] + CalcInfo_->num_drc_orbs;
       CalcInfo_->e0 += CalcInfo_->scfeigval[j];
       }

    for (i=0; i<CalcInfo_->num_bet_expl; i++) {
       j = (strbet->occs)[i] + CalcInfo_->num_drc_orbs;
       CalcInfo_->e0 += CalcInfo_->scfeigval[j];
       }
    }

   /* prepare the H0 block */

   Hd.diag_mat_els(alplist_, betlist_, CalcInfo_->onel_ints->pointer(),
          CalcInfo_->twoel_ints->pointer(), CalcInfo_->e0_drc, CalcInfo_->num_alp_expl,
          CalcInfo_->num_bet_expl, CalcInfo_->num_ci_orbs, Parameters_->hd_ave);

   H0block_setup(CIblks_->num_blocks, CIblks_->Ia_code, CIblks_->Ib_code);

   mpn_generator(Hd);
}


void CIWavefunction::mpn_generator(CIvect &Hd)
{
  double *mpk_energy, *mp2k_energy, *oei, *tei, *buffer1, *buffer2;
  double tval, Empn = 0.0, **wfn_overlap, Empn2 = 0.0, **cvec_coeff;
  double Empn2a = 0.0;
  double *cvec_norm, *tmp_coeff, tmp_norm, max_overlap = 1.0;
  int i, k, did_vec=0;
  int kvec_offset; /* offset if c_0 is not stored on disk */


  CIvect Cvec(Parameters_->icore, Parameters_->maxnvect, 1,
           Parameters_->c_filenum, CIblks_, CalcInfo_, Parameters_, H0block_, false);
  CIvect Sigma(Parameters_->icore, 1, 1,
            Parameters_->s_filenum, CIblks_, CalcInfo_, Parameters_, H0block_, false);
  CIvect Cvec2(Parameters_->icore, Parameters_->maxnvect, 1,
            Parameters_->c_filenum, CIblks_, CalcInfo_, Parameters_, H0block_, false);


   // setup I/O files, don't open old versions of these files
   Cvec.init_io_files(false);
   Sigma.init_io_files(false);

  /* set up the vector pointers/info */
  if (Cvec.read_new_first_buf() == -1) Cvec.write_new_first_buf();
  if (Sigma.read_new_first_buf() == -1) Sigma.write_new_first_buf();


  Cvec.h0block_buf_init();
  buffer1 = *(Hd.blockptr(0));
  buffer2 = Hd.buf_malloc();
  Hd.buf_unlock();

  wfn_overlap = init_matrix(Parameters_->maxnvect+1, Parameters_->maxnvect+1);
  cvec_coeff = init_matrix(Parameters_->maxnvect+1, Parameters_->maxnvect+1);
  tmp_coeff = init_array(Parameters_->maxnvect+1);
  cvec_norm = init_array(Parameters_->maxnvect+1);
  mpk_energy = init_array(Parameters_->maxnvect+1);
  mp2k_energy = init_array(2*(Parameters_->maxnvect+1)+1);
  mpk_energy[0] = mp2k_energy[0] = CalcInfo_->e0;
  mpk_energy[1] = CalcInfo_->e1 = mp2k_energy[1] =
       CalcInfo_->escf - CalcInfo_->e0 - CalcInfo_->enuc;
  if (CalcInfo_->iopen) kvec_offset = 0;
  else kvec_offset = 0;

  for (i=0; i<=Parameters_->maxnvect; i++) cvec_coeff[i][i] = 1.0;
  cvec_norm[0] = 1.0;
  wfn_overlap[0][0] = 1.0;

  /* oei = CalcInfo_->tf_onel_ints; */
  if (Parameters_->fci) oei = CalcInfo_->tf_onel_ints->pointer();
  else oei = CalcInfo_->gmat->pointer();
  tei = CalcInfo_->twoel_ints->pointer();

  outfile->Printf("   CalcInfo_->escf = %25.15f\n", CalcInfo_->escf);
  outfile->Printf("   CalcInfo_->e0   = %25.15f\n", CalcInfo_->e0);
  outfile->Printf("   CalcInfo_->enuc = %25.15f\n", CalcInfo_->enuc);
  outfile->Printf("   CalcInfo_->e1   = %25.15f\n\n", CalcInfo_->e1);
  if(Parameters_->zaptn) {
    outfile->Printf("   n         Corr. Energy                  E(ZAPTn)          "
        "   n         Corr. Energy                  E(ZAPTn)\n\n");
    } else {
    outfile->Printf("   n         Corr. Energy                  E(MPn)            "
        "   n         Corr. Energy                  E(MPn)\n\n");
    }
  outfile->Printf("   0  %25.15f %25.15f\n", 0.0000000000,
      CalcInfo_->e0+CalcInfo_->enuc);
  outfile->Printf("   1  %25.15f %25.15f\n",mpk_energy[1],
      mpk_energy[0]+mpk_energy[1]+CalcInfo_->enuc);
  Empn = mpk_energy[0]+mpk_energy[1]+CalcInfo_->enuc;
  Empn2 = Empn;


  Cvec.buf_lock(buffer1);

  Cvec.h0block_buf_init();
  tval = 1.0;
  Cvec.init_vals(0, 1, &(CalcInfo_->ref_alp_list), &(CalcInfo_->ref_alp_rel),
      &(CalcInfo_->ref_bet_list), &(CalcInfo_->ref_bet_rel), H0block_->blknum,
      &tval);
  if (print_ >= 5) {
    outfile->Printf("Zeroth-order wavefunction\n");
    Cvec.print();
  }

  Sigma.buf_lock(buffer2);
  Cvec.read(0, 0);  /* Set Cvec up correctly ? */

  //outfile->Printf("Cvec zero_blocks:\n");
  //Cvec.print_zero_blocks();
  Sigma.set_zero_blocks_all();

  //outfile->Printf("Sigma zero_blocks after set_zero_blocks_all.\n");
  //Sigma.print_zero_blocks();

  sigma(Cvec, Sigma, oei, tei, 0);
  if (print_ >= 5) {
    outfile->Printf("Sigma vector for 0 C vector\n");
    Sigma.read(0,0);
    Sigma.print();
    }
  //outfile->Printf("Sigma zero_blocks after sigma call.\n");
  //Sigma.print_zero_blocks();

  Cvec.read(0,0);  /* Set Cvec up correctly ? */
  Sigma.read(0,0); /* Set Sigma up correctly ? */
  tval = Cvec * Sigma;

  //outfile->Printf(" CalcInfo_->enuc = %25.15f\n", CalcInfo_->enuc);
  //outfile->Printf(" <psi0|Hc|psi0> = %25.15f\n", tval);
  //outfile->Printf(" CalcInfo_->edrc = %25.15f\n", CalcInfo_->edrc);
  //outfile->Printf(" mpk_energy[0] = %25.15f\n", mpk_energy[0]);

  tval += CalcInfo_->edrc - mpk_energy[0];

  outfile->Printf("   1  %25.15f %25.15f\n", tval,
   tval+mpk_energy[0]+CalcInfo_->enuc);
  if (tval - mpk_energy[1] > ZERO)
    outfile->Printf( "First-order energies do not agree!\n");


  Cvec.copy_zero_blocks(Sigma); /* Probably don't need this anymore */
  Cvec.copy(Sigma, (1-kvec_offset), 0);
  if (print_ >= 5) {
    outfile->Printf( "Cvec copying Sigma.\n");
    Cvec.print();
  }

  Sigma.buf_unlock();
  Hd.buf_lock(buffer2);

  tval = Cvec.dcalc2(1, CalcInfo_->e0, Hd, 1, alplist_, betlist_);
  Hd.buf_unlock();

  tval = 0.0;
  if (print_ >= 5) {
    outfile->Printf( "Cvec after dcalc2.\n");
    Cvec.print();
  }
  Cvec.read((1-kvec_offset),0);
  Cvec.set_vals((1-kvec_offset), 1, &(CalcInfo_->ref_alp_list),
          &(CalcInfo_->ref_alp_rel), &(CalcInfo_->ref_bet_list),
          &(CalcInfo_->ref_bet_rel), H0block_->blknum, &tval);
  Sigma.buf_lock(buffer2);

  if (print_ >= 5) {
    outfile->Printf( "Cvec after set_vals.\n");
    Cvec.print();
    }

  /* Here buffer1 = Cvec and buffer2 = Sigma */
  k=1;
  while (k<Parameters_->maxnvect) {

     /* Form Sigma */
     Sigma.read(0, 0);
     Cvec.read((k-kvec_offset), 0);  /* Set Cvec up correctly ? */
     sigma(Cvec, Sigma, oei, tei, 0);
     //outfile->Printf("Sigma zero_blocks after sigma call.\n");
     //Sigma.print_zero_blocks();

     /* Compute k+1, 2k, and 2k+1 th order energies from kth order wavefunction */
     Sigma.read(0,0); /* S_k is located in first Sigma space */
     //Sigma.print(outfile);

     #ifdef DEBUG
       tval = Sigma * Sigma;
       outfile->Printf(" Sigma * Sigma = %20.10f\n", tval);
       tval = sqrt(tval);
       outfile->Printf(" norm of Sigma = %20.10f\n", tval);
     #endif

     if (CalcInfo_->iopen) {
       Cvec.read(0,0); /* E_k = C_0 x S_k */
       tval = Cvec * Sigma;
       }
     else Sigma.extract_vals(0, 1, &(CalcInfo_->ref_alp_list),
          &(CalcInfo_->ref_alp_rel), &(CalcInfo_->ref_bet_list),
          &(CalcInfo_->ref_bet_rel), H0block_->blknum, &tval);
     mpk_energy[k+1] = tval;
     Empn += tval;
     outfile->Printf("  %2d  %25.15f %25.15f", k+1, mpk_energy[k+1], Empn);

     std::string label = (Parameters_->zaptn) ? "ZAPT" : "MP";

     Sigma.buf_unlock();
     if (Parameters_->wigner) {
       Cvec.wigner_E2k_formula(Hd, Sigma, Cvec2, alplist_, betlist_, buffer1,
          buffer2, k, mp2k_energy, wfn_overlap, cvec_coeff, cvec_norm, kvec_offset);
       Empn2 += mp2k_energy[2*k];
       outfile->Printf("      %2d %25.15f %25.15f\n",2*k,mp2k_energy[2*k],Empn2);

       /*- strings so that variable-name psi variables get parsed in docs -*/
       /*- Process::environment.globals["ZAPTn TOTAL ENERGY"] -*/
       /*- Process::environment.globals["MPn TOTAL ENERGY"] -*/
       /*- Process::environment.globals["ZAPTn CORRELATION ENERGY"] -*/
       /*- Process::environment.globals["MPn CORRELATION ENERGY"] -*/
       std::stringstream s;
       s << label << (2*k) << " TOTAL ENERGY";
       Process::environment.globals[s.str()] = Empn2;
       s.str(std::string());
       s << label << (2*k) << " CORRELATION ENERGY";
       Process::environment.globals[s.str()] = Empn2 - CalcInfo_->escf;
       //s.str(std::string());
       //s << label << (2*k) << " CORRECTION ENERGY";
       //Process::environment.globals[s.str()] = mp2k_energy[2*k];

       /* 25 November 2003 - JMT
        * Moified to save MP(2n-2) energy */
       Empn2a = Empn2;

       Empn2 += mp2k_energy[2*k+1];
       outfile->Printf("%62s %2d %25.15f %25.15f\n", "", 2*k+1, mp2k_energy[2*k+1], Empn2);

       s.str(std::string());
       s << label << (2*k+1) << " TOTAL ENERGY";
       Process::environment.globals[s.str()] = Empn2;
       s.str(std::string());
       s << label << (2*k+1) << " CORRELATION ENERGY";
       Process::environment.globals[s.str()] = Empn2 - CalcInfo_->escf;
       //s.str(std::string());
       //s << label << (2*k+1) << " CORRECTION ENERGY";
       //Process::environment.globals[s.str()] = mp2k_energy[2*k+1];

       }
     else {
       outfile->Printf( "\n");

       std::stringstream s;
       s << label << (k+1) << " TOTAL ENERGY";
       Process::environment.globals[s.str()] = Empn;
       s.str(std::string());
       s << label << (k+1) << " CORRELATION ENERGY";
       Process::environment.globals[s.str()] = Empn - CalcInfo_->escf;
       //s.str(std::string());
       //s << label << (k+1) << " CORRECTION ENERGY";
       //Process::environment.globals[s.str()] = mpk_energy[k+1];
     }



     if (k+1 == Parameters_->maxnvect) break;

     /* Schmidt orthogonalize vector space */
     if (Parameters_->mpn_schmidt && k>1) {
       Cvec2.buf_lock(buffer2);
       if (Cvec.schmidt_add2(Cvec2,0,k-2,k-1,k-1,cvec_coeff[k-1],
           (&cvec_norm[k-1]),&max_overlap)) did_vec = 1;
       else {
           std::string str = std::to_string(13);
           str += " vector norm = ";
           char*str2 = new char[25];
           sprintf(str2,"%20.15lf",cvec_norm[k-1]);
           str += str2;
           str += " < ";
           sprintf(str2,"%20.15lf",MPn_NORM_TOL);
           str += str2;
           delete[] str2;
           throw PsiException(str,__FILE__,__LINE__);
           }
       while (max_overlap > MPn_ZERO) {
         outfile->Printf("Second Schmidt-Orthogonalization performed.\n");
         Cvec.read(k-1,0);
         tval = Cvec * Cvec;
         tval = 1.0/sqrt(tval);
         outfile->Printf("Norm constant for %d S.O. Cvec = %20.15f\n",
                 k-1, tval);
         Cvec2.read(k-1,0);
         tval = Cvec2 * Cvec2;
         tval = 1.0/sqrt(tval);
         outfile->Printf("Norm constant for %d S.O. Cvec2 = %20.15f\n",
                 k-1, tval);
         if (Cvec.schmidt_add2(Cvec2,0,k-2,k-1,k-1,tmp_coeff,
            &tmp_norm,&max_overlap)) did_vec = 1;
         else {
           std::string str = std::to_string(13);
           str += " vector norm = ";
           char*str2 = new char[25];
           sprintf(str2,"%20.15lf",cvec_norm[k-1]);
           str += str2;
           str += " < ";
           sprintf(str2,"%20.15lf",MPn_NORM_TOL);
           str += str2;
           delete[] str2;
           throw PsiException(str,__FILE__,__LINE__);
           }
         for (i=0; i<k-1; i++) {
            outfile->Printf( "second coeff[%d] = %20.15f\n",i,tmp_coeff[i]);
            cvec_coeff[k-1][i] += (tmp_coeff[i]/cvec_norm[k-1]);
            }
         outfile->Printf( "max_overlap = %20.15f\n", max_overlap);
         outfile->Printf( "cvec_norm = %20.15f\n",cvec_norm[k-1]);
         outfile->Printf( "tmp_norm = %20.15f\n",tmp_norm);
         cvec_norm[k-1] *= tmp_norm;
         outfile->Printf( "cvec_norm new = %20.15f\n",cvec_norm[k-1]);
         }
       Cvec2.buf_unlock();
       }

     /* Construct k+1th order wavefunction */
     Cvec.construct_kth_order_wf(Hd, Sigma, Cvec2, alplist_, betlist_, buffer1,
        buffer2, k+1, mpk_energy, cvec_coeff, cvec_norm);

     // outfile->Printf( "Cvec %d = \n", k+1);
     // Cvec.print(outfile);

     if (Parameters_->mpn_schmidt) {
       outfile->Printf( "cvec_coeff = \n");
       print_mat(cvec_coeff, k-1, k-1, "outfile");
       }

     tval = 0.0;
     Cvec.set_vals(k+1, 1, &(CalcInfo_->ref_alp_list), &(CalcInfo_->ref_alp_rel),
             &(CalcInfo_->ref_bet_list), &(CalcInfo_->ref_bet_rel), H0block_->blknum,
             &tval);
     Sigma.buf_lock(buffer2);
     k++;
     Cvec.copy_zero_blocks(Sigma);
     }

   /* 22 Nov 2003 - JMT
    * Save the MPn or MP(2n-1) energy
    */
   if (Parameters_->save_mpn2 == 1 && Parameters_->wigner) {
     Process::environment.globals["CURRENT ENERGY"] = Empn2;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = Empn2 - Process::environment.globals["CURRENT REFERENCE ENERGY"];

     if(Parameters_->zaptn)
       outfile->Printf( "\n    ZAPT%d energy saved\n", (Parameters_->maxnvect * 2) - 1);
     else
       outfile->Printf( "\n    MP%d energy saved\n", (Parameters_->maxnvect * 2) - 1);
   }
   else if (Parameters_->save_mpn2 == 2 && Parameters_->wigner) {
     Process::environment.globals["CURRENT ENERGY"] = Empn2a;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = Empn2a - Process::environment.globals["CURRENT REFERENCE ENERGY"];
     if(Parameters_->zaptn)
       outfile->Printf( "\n    ZAPT%d energy saved\n", (Parameters_->maxnvect * 2) - 2);
     else
       outfile->Printf( "\n    MP%d energy saved\n", (Parameters_->maxnvect * 2) - 2);
   }
   else {
     Process::environment.globals["CURRENT ENERGY"] = Empn;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = Empn - Process::environment.globals["CURRENT REFERENCE ENERGY"];
     if(Parameters_->zaptn)
       outfile->Printf( "\n    ZAPT%d energy saved\n", Parameters_->maxnvect);
     else
       outfile->Printf( "\n    MP%d energy saved\n", Parameters_->maxnvect);
   }
   if(Parameters_->zaptn)
     outfile->Printf( "\n    EZAPTn = %17.13lf\n", Process::environment.globals["CURRENT ENERGY"]);
   else
     outfile->Printf( "\n    EMPn = %17.13lf\n", Process::environment.globals["CURRENT ENERGY"]);


   outfile->Printf("\n");
}

}} // namespace psi::detci
