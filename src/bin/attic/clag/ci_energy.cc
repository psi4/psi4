/*! \file
    \ingroup CLAG
    \brief Compute the CI Lagrangian
*/

/*! \defgroup CLAG clag: Compute the CI Lagrangian */

/*****************************************************************************/
/*ci_energy - This is a program to calculate the CI energy for a CI          */
/*         or MCSCF wavefunction given the one- and two-particle density     */
/*         matrix, the one-electron MO integrals Him, and the two-electron   */
/*         MO integrals (im,kl)                                              */ 
/*                                                                           */
/* Brian Hoffman                                                             */
/* Matt Leininger                                                            */
/*****************************************************************************/

#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

extern "C" {
  extern FILE *infile;
  extern FILE *outfile;
}

namespace psi { namespace clag {

#define INDEX(x,y) ((x>y) ? ioff[x] + y : ioff[y] + x)
#define CI_DIFF 1.0E-10

extern int *ioff;
extern int print_lvl;

/*****************************************************************************/
/* The main function                                                         */
/*****************************************************************************/

void ci_energy(double **OPDM, double *TPDM, double *h, double *TwoElec,
             int nbf, double enuc, double eci_chkpt, double lagtr)

{
  int i,j;                /* indecies of lagrangian element  */
  int m,k,l;              /* indecies of integrals needed    */ 
  int ij,kl,ijkl;         /* integral indecies combined      */
  int Tij, Tkl, Tijkl;    /* TPDM indecies combined          */
  double OEsum = 0.0;     /* QjmHim sumed over MO index m    */
  double TEsum = 0.0;     /* Gjmkl(im,kl) summed over m,k,l  */
  double e_ci = 0.0;      /* the CI energy                   */
  double diff;            /* diff between e_ci and eci_chkpt */
 
    for (i=0; i<nbf; i++)
       for (j=0; j<nbf; j++)
         {
           /*
           ** calculate the one-electron contribution
           */
           OEsum += (OPDM[i][j] * h[INDEX(i,j)]);             

           /*
           ** calculate two-elec. contribution
           */
           for (k=0; k<nbf; k++)
              for (l=0; l<nbf; l++) 
                {
                  /*
                  ** calculate integral composite index
                  */
                  ij = INDEX(i,j);
                  kl = INDEX(k,l);
                  ijkl = INDEX(ij,kl);
              
                  /*
                  ** calculate TPDM composite index 
                  */
                  Tij = i * nbf + j;
                  Tkl = k * nbf + l;
                  Tijkl = INDEX(Tij,Tkl);
                  TEsum += (TPDM[Tijkl] * TwoElec[ijkl]); 
                } /* end l loop */ 
         } /* end j loop */

  e_ci = OEsum + TEsum + enuc;
  diff = e_ci - eci_chkpt; 
  if (fabs(diff) > CI_DIFF) {
    fprintf(outfile,
      "Calculated CI Energy differs from the CI Energy in checkpoint file\n");
    fprintf(outfile,"ECI Calc. = %lf\n", e_ci); 
    fprintf(outfile,"ECI Chkpt = %lf\n", eci_chkpt); 
    }
 
  if (print_lvl > 0) {
    fprintf(outfile,"\nCheck CI Energy\n\n");
    fprintf(outfile,"One-electron contribution = %20.10lf\n", OEsum); 
    fprintf(outfile,"Two-electron contribution = %20.10lf\n", TEsum);
    fprintf(outfile,"Total electronic energy   = %20.10lf\n", e_ci);
    fprintf(outfile,"Trace of lagrangian       = %20.10lf\n", lagtr);
    fprintf(outfile,"Nuclear repulsion energy  = %20.10f\n", enuc); 
    fprintf(outfile,"Total CI Energy           = %20.10lf\n", e_ci); 
    fprintf(outfile,"CI Energy from chkpt file = %20.10lf\n", eci_chkpt);
    }          
}

}} // end namespace psi::clag

