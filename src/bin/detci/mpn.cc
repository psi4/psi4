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
#include <boost/lexical_cast.hpp>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include "structs.h"
#include "ci_tol.h"
#include <libpsio/psio.h>
#define EXTERN
#include "globals.h"
#include "civect.h"

namespace psi { namespace detci {

extern void print_vec(unsigned int nprint, int *Iacode, int *Ibcode,
   int *Iaidx, int *Ibidx, double *coeff,
   struct olsen_graph *AlphaG, struct olsen_graph *BetaG,
   struct stringwr **alplist, struct stringwr **betlist,
   std::string );


void mpn_generator(CIvect &Hd, struct stringwr **alplist, 
      struct stringwr **betlist)
{
  double *mpk_energy, *mp2k_energy, *oei, *tei, *buffer1, *buffer2;
  double tval, Empn = 0.0, **wfn_overlap, Empn2 = 0.0, **cvec_coeff;
  double Empn2a = 0.0;
  double *cvec_norm, norm, *tmp_coeff, tmp_norm, max_overlap = 1.0;
  int i, j, k, order, did_vec=0;
  int kvec_offset; /* offset if c_0 is not stored on disk */

  CIvect Cvec;
  CIvect Cvec2;
  CIvect Sigma;

  Cvec.set(CIblks.vectlen,CIblks.num_blocks,Parameters.icore,Parameters.Ms0,
     CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
     CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
     CalcInfo.nirreps, AlphaG->subgr_per_irrep, Parameters.maxnvect,
     Parameters.num_c_tmp_units, Parameters.first_c_tmp_unit,
     CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
  Sigma.set(CIblks.vectlen,CIblks.num_blocks,Parameters.icore,Parameters.Ms0,
     CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
     CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
     CalcInfo.nirreps, AlphaG->subgr_per_irrep, 1,
     Parameters.num_s_tmp_units, Parameters.first_s_tmp_unit,
     CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
  Cvec2.set(CIblks.vectlen,CIblks.num_blocks,Parameters.icore,Parameters.Ms0,
     CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
     CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
     CalcInfo.nirreps, AlphaG->subgr_per_irrep, Parameters.maxnvect,
     Parameters.num_c_tmp_units, Parameters.first_c_tmp_unit,
     CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);

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
  
  wfn_overlap = init_matrix(Parameters.maxnvect+1, Parameters.maxnvect+1);
  cvec_coeff = init_matrix(Parameters.maxnvect+1, Parameters.maxnvect+1);
  tmp_coeff = init_array(Parameters.maxnvect+1);
  cvec_norm = init_array(Parameters.maxnvect+1);
  mpk_energy = init_array(Parameters.maxnvect+1);
  mp2k_energy = init_array(2*(Parameters.maxnvect+1)+1);
  mpk_energy[0] = mp2k_energy[0] = CalcInfo.e0;
  mpk_energy[1] = CalcInfo.e1 = mp2k_energy[1] =  
       CalcInfo.escf - CalcInfo.e0 - CalcInfo.enuc;
  if (CalcInfo.iopen) kvec_offset = 0;
  else kvec_offset = 0;
 
  for (i=0; i<=Parameters.maxnvect; i++) cvec_coeff[i][i] = 1.0; 
  cvec_norm[0] = 1.0;
  wfn_overlap[0][0] = 1.0;

  /* oei = CalcInfo.tf_onel_ints; */
  if (Parameters.fci) oei = CalcInfo.tf_onel_ints;
  else oei = CalcInfo.gmat[0];
  tei = CalcInfo.twoel_ints;

  outfile->Printf("   CalcInfo.escf = %25.15f\n", CalcInfo.escf);
  outfile->Printf("   CalcInfo.e0   = %25.15f\n", CalcInfo.e0);
  outfile->Printf("   CalcInfo.enuc = %25.15f\n", CalcInfo.enuc);
  outfile->Printf("   CalcInfo.e1   = %25.15f\n\n", CalcInfo.e1);
  if(Parameters.zaptn) {
    outfile->Printf("   n         Corr. Energy \t\t E(ZAPTn) \t\t"
        "   n         Corr. Energy \t\t E(ZAPTn)\n\n");
    } else {
    outfile->Printf("   n         Corr. Energy \t\t E(MPn) \t\t"
        "   n         Corr. Energy \t\t E(MPn)\n\n");
    }
  outfile->Printf("   0  %25.15f %25.15f\n", 0.0000000000,
      CalcInfo.e0+CalcInfo.enuc);
  outfile->Printf("   1  %25.15f %25.15f\n",mpk_energy[1],
      mpk_energy[0]+mpk_energy[1]+CalcInfo.enuc);
  Empn = mpk_energy[0]+mpk_energy[1]+CalcInfo.enuc;
  Empn2 = Empn;
  

  Cvec.buf_lock(buffer1);

  Cvec.h0block_buf_init();
  tval = 1.0;
  Cvec.init_vals(0, 1, &(CalcInfo.ref_alp_list), &(CalcInfo.ref_alp_rel),
      &(CalcInfo.ref_bet_list), &(CalcInfo.ref_bet_rel), H0block.blknum, 
      &tval);
  if (Parameters.print_lvl >= 5) {
    outfile->Printf("Zeroth-order wavefunction\n");
    Cvec.print("outfile");
  }

  Sigma.buf_lock(buffer2);
  Cvec.read(0, 0);  /* Set Cvec up correctly ? */
  
  //outfile->Printf("Cvec zero_blocks:\n");
  //Cvec.print_zero_blocks(); 
  Sigma.set_zero_blocks_all();
  
  //outfile->Printf("Sigma zero_blocks after set_zero_blocks_all.\n");
  //Sigma.print_zero_blocks(); 
  
  sigma(alplist, betlist, Cvec, Sigma, oei, tei, Parameters.fci, 0); 
  if (Parameters.print_lvl >= 5) {
    outfile->Printf("Sigma vector for 0 C vector\n");
    Sigma.read(0,0);
    Sigma.print("outfile");
    }
  //outfile->Printf("Sigma zero_blocks after sigma call.\n");
  //Sigma.print_zero_blocks(); 

  Cvec.read(0,0);  /* Set Cvec up correctly ? */
  Sigma.read(0,0); /* Set Sigma up correctly ? */
  tval = Cvec * Sigma;

  //outfile->Printf(" CalcInfo.enuc = %25.15f\n", CalcInfo.enuc);
  //outfile->Printf(" <psi0|Hc|psi0> = %25.15f\n", tval);
  //outfile->Printf(" CalcInfo.efzc = %25.15f\n", CalcInfo.efzc);
  //outfile->Printf(" mpk_energy[0] = %25.15f\n", mpk_energy[0]);
 
  tval += CalcInfo.efzc - mpk_energy[0]; 
 
  outfile->Printf("   1  %25.15f %25.15f\n", tval,
   tval+mpk_energy[0]+CalcInfo.enuc);
  if (tval - mpk_energy[1] > ZERO) 
    outfile->Printf( "First-order energies do not agree!\n");
  

  Cvec.copy_zero_blocks(Sigma); /* Probably don't need this anymore */ 
  Cvec.copy(Sigma, (1-kvec_offset), 0);
  if (Parameters.print_lvl >= 5) {
    outfile->Printf( "Cvec copying Sigma.\n");
    Cvec.print("outfile");
  }

  Sigma.buf_unlock();
  Hd.buf_lock(buffer2); 

  tval = Cvec.dcalc2(1, CalcInfo.e0, Hd, 1, alplist, betlist);
  Hd.buf_unlock();

  tval = 0.0;
  if (Parameters.print_lvl >= 5) {
    outfile->Printf( "Cvec after dcalc2.\n");
    Cvec.print("outfile");
  }
  Cvec.read((1-kvec_offset),0);
  Cvec.set_vals((1-kvec_offset), 1, &(CalcInfo.ref_alp_list), 
          &(CalcInfo.ref_alp_rel), &(CalcInfo.ref_bet_list), 
          &(CalcInfo.ref_bet_rel), H0block.blknum, &tval);
  Sigma.buf_lock(buffer2);

  if (Parameters.print_lvl >= 5) {
    outfile->Printf( "Cvec after set_vals.\n");
    Cvec.print("outfile");
    }

  /* Here buffer1 = Cvec and buffer2 = Sigma */
  k=1; 
  while (k<Parameters.maxnvect) {

     /* Form Sigma */
     Sigma.read(0, 0);
     Cvec.read((k-kvec_offset), 0);  /* Set Cvec up correctly ? */
     sigma(alplist, betlist, Cvec, Sigma, oei, tei, Parameters.fci, 0);
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

     if (CalcInfo.iopen) { 
       Cvec.read(0,0); /* E_k = C_0 x S_k */
       tval = Cvec * Sigma;
       }
     else Sigma.extract_vals(0, 1, &(CalcInfo.ref_alp_list),
          &(CalcInfo.ref_alp_rel), &(CalcInfo.ref_bet_list),
          &(CalcInfo.ref_bet_rel), H0block.blknum, &tval);
     mpk_energy[k+1] = tval; 
     Empn += tval;
     outfile->Printf("  %2d  %25.15f %25.15f", k+1, mpk_energy[k+1], Empn); 
 
     std::string label = (Parameters.zaptn) ? "ZAPT" : "MP";

     Sigma.buf_unlock();
     if (Parameters.wigner) {
       Cvec.wigner_E2k_formula(Hd, Sigma, Cvec2, alplist, betlist, buffer1,
          buffer2, k, mp2k_energy, wfn_overlap, cvec_coeff, cvec_norm, kvec_offset);
       Empn2 += mp2k_energy[2*k];
       outfile->Printf("\t %2d %25.15f %25.15f\n",2*k,mp2k_energy[2*k],Empn2);

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
       Process::environment.globals[s.str()] = Empn2 - CalcInfo.escf;
       //s.str(std::string());
       //s << label << (2*k) << " CORRECTION ENERGY";
       //Process::environment.globals[s.str()] = mp2k_energy[2*k];

       /* 25 November 2003 - JMT
        * Moified to save MP(2n-2) energy */
       Empn2a = Empn2;

       Empn2 += mp2k_energy[2*k+1];
       outfile->Printf("\t\t\t\t\t\t\t\t"
               " %2d %25.15f %25.15f\n", 2*k+1, mp2k_energy[2*k+1], Empn2);

       s.str(std::string());
       s << label << (2*k+1) << " TOTAL ENERGY";
       Process::environment.globals[s.str()] = Empn2;
       s.str(std::string());
       s << label << (2*k+1) << " CORRELATION ENERGY";
       Process::environment.globals[s.str()] = Empn2 - CalcInfo.escf;
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
       Process::environment.globals[s.str()] = Empn - CalcInfo.escf;
       //s.str(std::string());
       //s << label << (k+1) << " CORRECTION ENERGY";
       //Process::environment.globals[s.str()] = mpk_energy[k+1];
     }

     

     if (k+1 == Parameters.maxnvect) break;

     /* Schmidt orthogonalize vector space */
     if (Parameters.mpn_schmidt && k>1) {
       Cvec2.buf_lock(buffer2);
       if (Cvec.schmidt_add2(Cvec2,0,k-2,k-1,k-1,cvec_coeff[k-1],
           (&cvec_norm[k-1]),&max_overlap)) did_vec = 1;
       else {
           std::string str = boost::lexical_cast<std::string>( 13);
           str += " vector norm = "; 
           char*str2 = new char[25];
           sprintf(str2,"%20.15lf",cvec_norm[k-1]);
           str += str2;
           str += " < ";
           sprintf(str2,"%20.15lf",MPn_NORM_TOL);
           str += str2;
           delete str2;
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
           std::string str = boost::lexical_cast<std::string>(13);
           str += " vector norm = "; 
           char*str2 = new char[25];
           sprintf(str2,"%20.15lf",cvec_norm[k-1]);
           str += str2;
           str += " < ";
           sprintf(str2,"%20.15lf",MPn_NORM_TOL);
           str += str2;
           delete str2;
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
     Cvec.construct_kth_order_wf(Hd, Sigma, Cvec2, alplist, betlist, buffer1,
        buffer2, k+1, mpk_energy, cvec_coeff, cvec_norm);
     
     // outfile->Printf( "Cvec %d = \n", k+1);
     // Cvec.print(outfile);

     if (Parameters.mpn_schmidt) {
       outfile->Printf( "cvec_coeff = \n");
       print_mat(cvec_coeff, k-1, k-1, "outfile");
       }

     tval = 0.0;
     Cvec.set_vals(k+1, 1, &(CalcInfo.ref_alp_list), &(CalcInfo.ref_alp_rel),
             &(CalcInfo.ref_bet_list), &(CalcInfo.ref_bet_rel), H0block.blknum,
             &tval);
     Sigma.buf_lock(buffer2);
     k++;
     Cvec.copy_zero_blocks(Sigma); 
     }

   /* 22 Nov 2003 - JMT 
    * Save the MPn or MP(2n-1) energy
    */
   chkpt_init(PSIO_OPEN_OLD);
   if (Parameters.save_mpn2 == 1 && Parameters.wigner) {
     chkpt_wt_etot(Empn2);
     Process::environment.globals["CURRENT ENERGY"] = Empn2;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = Empn2 - Process::environment.globals["CURRENT REFERENCE ENERGY"];

     if(Parameters.zaptn) 
       outfile->Printf( "\nZAPT%d energy saved\n", (Parameters.maxnvect * 2) - 1);
     else
       outfile->Printf( "\nMP%d energy saved\n", (Parameters.maxnvect * 2) - 1);
   }
   else if (Parameters.save_mpn2 == 2 && Parameters.wigner) {
     chkpt_wt_etot(Empn2a);
     Process::environment.globals["CURRENT ENERGY"] = Empn2a;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = Empn2a - Process::environment.globals["CURRENT REFERENCE ENERGY"];
     if(Parameters.zaptn)
       outfile->Printf( "\nZAPT%d energy saved\n", (Parameters.maxnvect * 2) - 2);
     else
       outfile->Printf( "\nMP%d energy saved\n", (Parameters.maxnvect * 2) - 2);
   }
   else {
     chkpt_wt_etot(Empn);
     Process::environment.globals["CURRENT ENERGY"] = Empn;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = Empn - Process::environment.globals["CURRENT REFERENCE ENERGY"];
     if(Parameters.zaptn)
       outfile->Printf( "\nZAPT%d energy saved\n", Parameters.maxnvect);
     else
       outfile->Printf( "\nMP%d energy saved\n", Parameters.maxnvect);
   }
   if(Parameters.zaptn)
     outfile->Printf( "\nEZAPTn = %17.13lf\n", chkpt_rd_etot());
   else
     outfile->Printf( "\nEMPn = %17.13lf\n", chkpt_rd_etot());

   chkpt_close();

   outfile->Printf("\n");
}

}} // namespace psi::detci

