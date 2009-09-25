/*! 
  \file
  \ingroup DETCI
  \brief Check the SCF energy
*/
#include <cstdio>
#include <cmath>

namespace psi { namespace detci {

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

extern int *ioff ;

void scf_energy(double *H, double *TE, double *energy_1, double *energy_2, 
      double *energy_e, int *docc, int *frozen_docc, int fzc_flag, 
      int nirreps, int *reorder, int *opi);

/*!
** check_energy(): check the SCF energy by calculating it from the two-electr.
**    integrals in the MO basis
** 
** \param  H         =  lwr tri of one-electron integrals matrix (MO basis)
** \param twoel_ints =  two electron integrals (lexically indexed, MO basis)
** \param   nocc     =  num occupied orbitals (assume closed shell case) and
**                      exclude frozen core
** \param   escf     =  scf energy to compare to
** \param   enuc     =  nuclear repulsion energy
** \param   efzc     =  frozen core energy
** \param   nirreps  =  number of irreps 
** \param   reorder  =  reordering array for Pitzer->CI ordering
** \param   opi      =  orbs per irrep in Pitzer ordering
** \param   outfile  =  file to write output to
**
** Returns: the computed SCF energy
** \ingroup DETCI
*/
double check_energy(double *H, double *twoel_ints, int *docc, int *frozen_docc,
      int fzc_flag, double escf, double enuc, double efzc, 
      int nirreps, int *reorder, int *opi, int print_lvl, FILE *outfile)
{
   double energy_1 ;     /* one-electron energy */
   double energy_2 ;     /* two-electron energy */
   double energy_e ;     /* total electronic energy */

   scf_energy(H, twoel_ints, &energy_1, &energy_2, &energy_e, docc,
      frozen_docc, fzc_flag, nirreps, reorder, opi);

   if (print_lvl) {
     fprintf(outfile,"\nCheck SCF Energy from 1- and 2-electron integrals\n\n");
     fprintf(outfile,"SCF Energy (ref):          %16.10lf\n", escf) ;
     fprintf(outfile,"Nuclear repulsion energy:  %16.10lf\n", enuc) ;
     fprintf(outfile,"One-electron energy:       %16.10lf\n", energy_1) ;
     fprintf(outfile,"Two-electron energy:       %16.10lf\n", energy_2) ;
     fprintf(outfile,"Frozen core energy:        %16.10lf\n", efzc) ;
     fprintf(outfile,"Total electronic energy:   %16.10lf\n", energy_e+efzc) ;
     fprintf(outfile,"Total SCF energy:          %16.10lf\n", enuc + 
        energy_e + efzc) ;
    
     if (fabs(enuc + efzc + energy_e - escf) > 0.00000001) {
        fprintf(outfile, 
           "\n*** Calculated Energy Differs from SCF Energy in CHKPT ! ***\n") ;
        }
   }

   return(enuc+efzc+energy_e); 
}   


/*!
** scf_energy(): Function calculates the SCF energy from the one- and
**   two-electron integrals in MO form (closed-shell case).
**
** David Sherrill, Sept 1993
**
** \param H        = Matrix of one-electron integrals in MO basis (lwr triangle)
** \param TE       = Two-electron integrals in MO basis, stored in 
**                   ijkl-indexed array
** \param energy_1 = pointer to hold one-electron energy
** \param energy_2 = pointer to hold two-electron energy
** \param energy_e = pointer to hold total electronic energy (sum of two 
**                   terms above)
** \param docc     = array of doubly-occupied orbitals per irrep
** \param frozen_docc = array of frozen doubly-occupied orbitals per irrep
** \param fzc_flag = remove explicit consideration of frozen core orbitals ?
** \param nirreps  = number of irreps 
** \param reorder  = reordering array Pitzer->CI order
** \param opi      = orbitals per irrep
**
** Returns: none
**
** \ingroup DETCI
*/

void scf_energy(double *H, double *TE, double *energy_1, double *energy_2, 
      double *energy_e, int *docc, int *frozen_docc, int fzc_flag, 
      int nirreps, int *reorder, int *opi)
{
   int irrep, irrep2, d, d2, offset, offset2, ndoc, ndoc2, nfzc, nfzc2, totfzc;
   int i, j;
   int ii, jj, iijj, ij, ijij, iiii;

   *energy_1 = *energy_2 = *energy_e = 0.0;

   totfzc=0;
   if (fzc_flag) {
     for (irrep=0; irrep<nirreps; irrep++) {
       totfzc += frozen_docc[irrep];
     }
   }

   for (irrep=0,offset=0; irrep<nirreps; irrep++) {
     if (irrep>0) offset += opi[irrep-1];
     ndoc = docc[irrep];
     if (fzc_flag) {
       nfzc = frozen_docc[irrep];
       ndoc -= nfzc;
     }
     else nfzc=0;
     for (d=offset+nfzc; d<ndoc+offset+nfzc; d++) {
       i = reorder[d]-totfzc;
       ii = ioff[i] + i;
       iiii = ioff[ii] + ii;
       *energy_1 += 2.0 * H[ii];
       *energy_2 += TE[iiii];
       
       for (irrep2=0,offset2=0; irrep2<=irrep; irrep2++) {
	 if (irrep2>0) offset2 += opi[irrep2-1];
	 ndoc2 = docc[irrep2];
	 if (fzc_flag) {
	   nfzc2 = frozen_docc[irrep2];
	   ndoc2 -= nfzc2;
	 }
	 else nfzc2=0;
	 
	 for (d2=offset2+nfzc2; d2<ndoc2+offset2+nfzc2 && d2<d; d2++) {
	   j = reorder[d2]-totfzc;
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

