/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

/*! 
  \file
  \ingroup DETCI
  \brief Check the SCF energy
*/
#include <cstdio>
#include <cmath>
#include "psi4-dec.h"
#include "libparallel/ParallelPrinter.h"
namespace psi { namespace detci {

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

extern int *ioff ;

void scf_energy(double *H, double *TE, double *energy_1, double *energy_2, 
      double *energy_e, int *docc, int *dropped_docc, int drc_flag, 
      int nirreps, int *reorder, int *opi);

/*!
** check_energy(): check the SCF energy by calculating it from the two-electr.
**    integrals in the MO basis
** 
** \param H            =  lwr tri of one-electron integrals matrix (MO basis)
** \param twoel_ints   =  two electron integrals (lexically indexed, MO basis)
** \param docc         =  doubly occupied orbitals per irrep
** \param dropped_docc =  dropped occupied orbitals per irrep 
** \param drc_flag     =  1 if we drop core orbitals, 0 otherwise
** \param escf         =  scf energy to compare to
** \param enuc         =  nuclear repulsion energy
** \param edrc         =  energy of dropped core orbitals
** \param nirreps      =  number of irreps 
** \param eorder       =  reordering array for Pitzer->CI ordering
** \param opi          =  orbs per irrep in Pitzer ordering
** \param print_lvl    =  integer describing how much to print
** \param outfile      =  file to write output to
**
** Returns: the computed SCF energy
** \ingroup DETCI
*/
double check_energy(double *H, double *twoel_ints, int *docc, 
  int *dropped_docc, int drc_flag, double escf, double enuc, double edrc, 
      int nirreps, int *reorder, int *opi, int print_lvl, std::string out)
{
   double energy_1 ;     /* one-electron energy */
   double energy_2 ;     /* two-electron energy */
   double energy_e ;     /* total electronic energy */
   boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         boost::shared_ptr<OutFile>(new OutFile(out)));
   scf_energy(H, twoel_ints, &energy_1, &energy_2, &energy_e, docc,
      dropped_docc, drc_flag, nirreps, reorder, opi);

   if (print_lvl) {
     printer->Printf("\nCheck SCF Energy from 1- and 2-electron integrals\n\n");
     printer->Printf("SCF Energy (ref):          %16.10lf\n", escf) ;
     printer->Printf("Nuclear repulsion energy:  %16.10lf\n", enuc) ;
     printer->Printf("One-electron energy:       %16.10lf\n", energy_1) ;
     printer->Printf("Two-electron energy:       %16.10lf\n", energy_2) ;
     printer->Printf("Dropped core energy:       %16.10lf\n", edrc) ;
     printer->Printf("Total electronic energy:   %16.10lf\n", energy_e+edrc) ;
     printer->Printf("Total SCF energy:          %16.10lf\n", enuc +
        energy_e + edrc) ;
    
     if (fabs(enuc + edrc + energy_e - escf) > 0.00000001) {
        printer->Printf(
           "\n*** Calculated Energy Differs from SCF Energy in CHKPT ! ***\n") ;
        }
   }

   return(enuc+edrc+energy_e); 
}   


/*!
** scf_energy(): Function calculates the SCF energy from the one- and
**   two-electron integrals in MO form (closed-shell case).
**
** David Sherrill, Sept 1993
**
** \param H            = Matrix of one-electron integrals in MO basis 
**                       (lower triangle)
** \param TE           = Two-electron integrals in MO basis, stored in 
**                       ijkl-indexed array
** \param energy_1     = pointer to hold one-electron energy
** \param energy_2     = pointer to hold two-electron energy
** \param energy_e     = pointer to hold total electronic energy (sum of two 
**                       terms above)
** \param docc         = array of doubly-occupied orbitals per irrep
** \param dropped_docc = dropped occupied orbitals per irrep 
** \param drc_flag     = 1 if we drop core orbitals, 0 otherwise
** \param nirreps      = number of irreps 
** \param reorder      = reordering array Pitzer->CI order
** \param opi          = orbitals per irrep
**
** Returns: none
**
** \ingroup DETCI
*/

void scf_energy(double *H, double *TE, double *energy_1, double *energy_2, 
      double *energy_e, int *docc, int *dropped_docc, int drc_flag, 
      int nirreps, int *reorder, int *opi)
{
   int irrep, irrep2, d, d2, offset, offset2, ndoc, ndoc2, ndrc, ndrc2, totdrc;
   int i, j;
   int ii, jj, iijj, ij, ijij, iiii;

   *energy_1 = *energy_2 = *energy_e = 0.0;

   totdrc=0;
   if (drc_flag) {
     for (irrep=0; irrep<nirreps; irrep++) {
       totdrc += dropped_docc[irrep];
     }
   }

   for (irrep=0,offset=0; irrep<nirreps; irrep++) {
     if (irrep>0) offset += opi[irrep-1];
     ndoc = docc[irrep];
     if (drc_flag) {
       ndrc = dropped_docc[irrep];
       ndoc -= ndrc;
     }
     else ndrc=0;
     for (d=offset+ndrc; d<ndoc+offset+ndrc; d++) {
       i = reorder[d]-totdrc;
       ii = ioff[i] + i;
       iiii = ioff[ii] + ii;
       *energy_1 += 2.0 * H[ii];
       *energy_2 += TE[iiii];
       
       for (irrep2=0,offset2=0; irrep2<=irrep; irrep2++) {
	 if (irrep2>0) offset2 += opi[irrep2-1];
	 ndoc2 = docc[irrep2];
	 if (drc_flag) {
	   ndrc2 = dropped_docc[irrep2];
	   ndoc2 -= ndrc2;
	 }
	 else ndrc2=0;
	 
	 for (d2=offset2+ndrc2; d2<ndoc2+offset2+ndrc2 && d2<d; d2++) {
	   j = reorder[d2]-totdrc;
	   jj = ioff[j] + j;
	   iijj = INDEX(ii,jj);
	   ij = INDEX(i,j);
	   ijij = ioff[ij] + ij;
	   *energy_2 += 4.0 * TE[iijj] - 2.0 * TE[ijij];
	 }
       }
     }
   }
   
   *energy_e = *energy_1 + *energy_2;
}

}} // namespace psi::detci