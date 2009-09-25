/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here 
*/
/*
** FZC.C: Source code to handle a few things associated with transforming
**    integrals with frozen core orbitals.
**
**    C. David Sherrill
**    Center for Computational Quantum Chemistry
**    March, 1995
**
*/

#include <cstdio>

extern FILE *outfile;

namespace psi { namespace transqt {

/*
** FZC_DENSITY(): Form the frozen core density matrix, i.e. the density
**    matrix in symmetry orbitals which arises from the electrons occupying
**    frozen core orbitals.
**
**       Pc[mu][nu] = SUM_{i}^{n_core} C[mu][i] * C[nu][i]
**
** Arguments:
**    nirreps =  number of irreducible representations
**    docc    =  array of number of doubly occupied orbitals per irrep
**    Pc      =  frozen core density matrix (lower triangle)
**    C       =  SCF eigenvector matrix
**    first   =  first molecular orbital for each irrep
**    first_so = first SO index for each irrep
**    last_so =  last SO index for each irrep
**    ioff    =  the standard offset array
**
** Returns: none
**
** Notes: The routine assumes Pc has been allocated already.  The matrix
**    multiply uses sparse matrix techniques, but it might be rearranged
**    more efficiently, because C is accessed down rows instead of across
**    them.  Accessing of Pc is fine.  (Since this routine is only called
**    once, it may not be worth the effort to improve its efficiency).
*/
void fzc_density(int nirreps, int *docc, double *Pc, double **C,
      int *first, int *first_so, int *last_so, int *ioff)
{
   int p, q, pq;
   int ir, i, it, j;

   for (ir=0; ir<nirreps; ir++) {  /* loop over irreps */
      for (i=0,it=first[ir]; i<docc[ir]; i++,it++) { /* loop over fzc */
        for (p=first_so[ir]; p <= last_so[ir]; p++) {
           pq = ioff[p] + first_so[ir];
           for (q=first_so[ir]; q <= p; q++,pq++) {
              Pc[pq] += C[p][it] * C[q][it];
              }
           }   
        } /* end loop over fzc */
      } /* end loop over irreps */

}


/*
** IVO_DENSITY(): This is an attempt to form a N-1 electron density to 
**    be used to make an effective N-1 electron Fock matrix, using the
**    normal frozen core operator subroutine found in yosh_rdwto().
**    We take the modified density as the core density (formed by 
**    fzc_density()) and add the valence density, but scaled by a
**    factor of (N-1)/N.  I'm not sure about open-shells.
**
**       P[mu][nu] = SUM_{i}^{n_core} C[mu][i] * C[nu][i]
**                 + (N-1)/N * SUM_{i}^{valence} C[mu][i] * C[nu][i]
**
** Arguments:
**    nirreps =  number of irreducible representations
**    fzdocc  =  number of frozen core doubly occupied orbitals per irrep
**    docc    =  array of number of doubly occupied orbitals per irrep
**    socc    =  array of number of singly occupied orbitals per irrep
**    P       =  our special density matrix (lower triangle)
**    C       =  SCF eigenvector matrix
**    first   =  first molecular orbital for each irrep
**    first_so = first SO index for each irrep
**    last_so =  last SO index for each irrep
**    ioff    =  the standard offset array
**
** Returns: none
**
** Notes: The routine assumes P has been allocated already.
*/
void ivo_density(int nirreps, int *fzdocc, int *docc, int *socc, 
                 double *P, double **C, int *first, int *first_so, 
                 int *last_so, int *ioff)
{
   int p, q, pq;
   int ir, i, it, j;
   double N;
   double tval;

   N = 0.0;
   for (ir=0; ir<nirreps; ir++) {
     N += 2.0 * (double) (docc[ir] - fzdocc[ir]);
     N += (double) socc[ir];
   }


   fprintf(outfile, "%f valence electrons", N);

   for (ir=0; ir<nirreps; ir++) {  /* loop over irreps */
      for (i=0,it=first[ir]; i<fzdocc[ir]; i++,it++) { /* loop over fzc */
        for (p=first_so[ir]; p <= last_so[ir]; p++) {
           pq = ioff[p] + first_so[ir];
           for (q=first_so[ir]; q <= p; q++,pq++) {
              P[pq] += C[p][it] * C[q][it];
              }
           }   
        } /* end loop over fzc */

      /* loop over docc */
      for (i=0,it=first[ir]+fzdocc[ir]; i<docc[ir]-fzdocc[ir]; i++,it++) { 
        fprintf(outfile, "Doing docc orbital ir=%d, orb=%d\n", ir, i);
        for (p=first_so[ir]; p <= last_so[ir]; p++) {
           pq = ioff[p] + first_so[ir];
           for (q=first_so[ir]; q <= p; q++,pq++) {
              tval = (N-1)/N * C[p][it] * C[q][it];
              fprintf(outfile, "%lf\n", tval);
              P[pq] += (N-1)/N * C[p][it] * C[q][it];
              }
           }     
        } /* end loop over docc */

      /* I am very confused about socc, maybe this works? 
       * I suspect that the construction of the frozen core 
       * operator ONLY works correctly for closed shells (which
       * a frozen core always is closed shell...).
       * Hence, it is pretty much irrelevant to trying to 
       * get this density right, but here's something anyway --CDS
       */   
      for (i=0,it=first[ir]+docc[ir]; i<socc[ir]; i++,it++) { 
        fprintf(outfile, "Doing socc orbital ir=%d, orb=%d\n", ir, i);
        for (p=first_so[ir]; p <= last_so[ir]; p++) {
           pq = ioff[p] + first_so[ir];
           for (q=first_so[ir]; q <= p; q++,pq++) {
              P[pq] += 0.5 * (N-1)/N * C[p][it] * C[q][it];
              }
           }
        } /* end loop over socc */

      } /* end loop over irreps */

}



/*
** FZC_ENERGY(): This function calculates the frozen core energy, which is 
**    the expectation value of the determinant made up of the frozen core 
**    orbitals only.
**
** Arguments:
**    nbfso   = number of symmetry orbitals
**    orbsym  = array containing symmetry of each orbital
**    P       = frozen core density matrix
**    Hc      = frozen core operator in SO basis
**    H       = one-electron operator in SO basis
**    first_so = first SO index for each irrep
**    ioff    = the usual offset array
**
** Returns: the frozen core energy
*/
double fzc_energy(int nbfso, int *orbsym, double *P, double *Hc, double *H, 
      int *first_so, int *ioff)
{
   int i,j,ij,isym;
   double fzc = 0.0;

   for (i=0; i<nbfso; i++) {
      isym = orbsym[i];
 
      /* Off-diagonal elements of P */
      for (j=first_so[isym],ij=ioff[i]+first_so[isym]; j<i; j++,ij++) {
         fzc += 2.0 * P[ij] * (Hc[ij] + H[ij]);
         }

      /* Diagonal elements of P */
      ij = ioff[i] + i;
      fzc += P[ij] * (Hc[ij] + H[ij]);
      }

   return(fzc);
}

/*
** FZC_ENERGY_UHF(): This function calculates the UHF frozen core energy.
**
** Arguments:
**    nbfso   = number of symmetry orbitals
**    orbsym  = array containing symmetry of each orbital
**    Pa      = alpha frozen core density matrix
**    Pb      = beta frozen core density matrix
**    Hca     = alpha frozen core operator in SO basis
**    Hcb     = beta frozen core operator in SO basis
**    H       = one-electron operator in SO basis
**    first_so = first SO index for each irrep
**    ioff    = the usual offset array
**
** Returns: the frozen core energy
*/
double fzc_energy_uhf(int nbfso, int *orbsym, double *Pa, double *Pb,
		      double *Hca, double *Hcb, double *H, 
		      int *first_so, int *ioff)
{
   int i,j,ij,isym;
   double fzc = 0.0;

   for (i=0; i<nbfso; i++) {
      isym = orbsym[i];
 
      /* Off-diagonal elements of P */
      for (j=first_so[isym],ij=ioff[i]+first_so[isym]; j<i; j++,ij++) {
         fzc += Pa[ij] * (Hca[ij] + H[ij]);
	 fzc += Pb[ij] * (Hcb[ij] + H[ij]);
         }

      /* Diagonal elements of P */
      ij = ioff[i] + i;
      fzc += 0.5 * Pa[ij] * (Hca[ij] + H[ij]);
      fzc += 0.5 * Pb[ij] * (Hcb[ij] + H[ij]);
      }

   return(fzc);
}

}} // end namespace psi::transqt
