/*! \file
    \ingroup CLAG
    \brief Routine that does the work for the CI Lagrangian computation
*/
/*****************************************************************************/
/*lagcalc - This is a program to calculate the lagrangian matrix for a CI    */
/*         or MCSCF wavefunction given the one- and two-particle density     */
/*         matrix, the one-electron MO integrals Him, and the two-electron   */
/*         MO integrals (im, kl)                                             */
/*                                                                           */ 
/* Brian Hoffman                                                             */
/* Matt Leininger                                                            */
/*****************************************************************************/

#include <libqt/qt.h>         
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include "MCSCF.h"
#include "globaldefs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

// #define INDEX(x,y) ((x>y) ? ioff[x] + y : ioff[y] + x)

// Compute the Lagrangian matrix
//
// this code assumes h is the bare one-electron ints, not the frozen
// core operator [if desired, it is possible to write one that would take the 
// fzc operator, along the lines of equations 4.3 in Chaban, Schmidt, and
// Gordon, Theor. Chem. Acc. 97, 88-95 (1997), but the computational cost
// of this step is pretty negligible anyway] 
//
double** MCSCF::lagcalc(double **OPDM, double *TPDM, double *h, double *TwoElec,
                 int nmo, int npop, int print_lvl, int lag_file)
{
  int i,j;                      /* indecies of lagrangian element          */
  int m,k,l;                    /* indecies of integrals needed            */ 
  int im,kl,imkl;               /* integral indecies combined              */
  int Tjm,Tmj,Tkl,Tjmkl,Tmjkl;  /* TPDM indices combined                   */
  int mj, mjkl;                 /* revised TPDM indices                    */
  double **oe_lag;              /* One-electron part of lagrangian         */
  double **te_lag;              /* Two-electron part of lagrangian         */ 
  double OEsum;                 /* QjmHim sumed over MO index m            */
  double TEsum;                 /* Gjmkl(im,kl) summed over m,k,l          */
  double lagtrace;              /* Trace of the Lagrangian (as a check)    */
  double **lag;                 /* Lagrangian                              */



  lag = block_matrix(nmo, nmo);
  oe_lag = block_matrix(nmo, nmo);
  te_lag = block_matrix(nmo, nmo);

  for (i=0; i<nmo; i++)
    for (j=0; j<npop; j++)
       {
          /*
          ** zero locations of intermediate sums 
          */
          OEsum = 0.0;
          TEsum = 0.0;

          for (m=0; m<npop; m++)
            {
              /*
              ** calculate the one-electron contribution
              */
              OEsum += (OPDM[j][m] * h[INDEX(i,m)]);             
              for (k=0; k<npop; k++)
                for (l=0; l<npop; l++) 
                  {
                    /*
                    ** calculate two-elec. contribution and
                    ** calculate integral composite index
                    */ 
                    im = INDEX(i,m);
                    kl = INDEX(k,l);
                    imkl = INDEX(im,kl);
              
                    /*
                    ** Calculate TPDM composite index and
                    **  total sum.  If TPDM[ijkl] != TPDM[jikl] we
                    **  need two terms, unlike in Dr. Yamaguchi's book.
                    ** We use two terms below.
                    */
                    mj = INDEX(m,j);
                    mjkl = INDEX(mj,kl);
                    TEsum += TPDM[mjkl] * TwoElec[imkl]; 
                  } /* end l loop */
            } /* end m loop */

          oe_lag[i][j] = OEsum;
          te_lag[i][j] = TEsum;
          lag[i][j] = OEsum + TEsum;
      
        } /* end j loop */ 

  /*
  ** check the trace of the Lagrangian (supposedly = energy)
  */

  if (print_lvl > 1) {
    outfile->Printf("\n\n One-electron part of the Lagrangian");
    print_mat(oe_lag, nmo, nmo, "outfile");
    outfile->Printf("\n\n Two-electron part of the Lagrangian");
    print_mat(te_lag, nmo, nmo, "outfile");
    outfile->Printf("\nLagrangian Matrix\n\n");
    print_mat(lag, nmo, nmo, "outfile");
    }

  for (lagtrace=0.0,i=0; i<nmo; i++)
    lagtrace += oe_lag[i][i] + 0.5 * te_lag[i][i];

  outfile->Printf("Lagrangian Trace is %6.10f\n", lagtrace);

  free_block(oe_lag);
  free_block(te_lag);

  return(lag);
}

}} // end namespace psi::detci

