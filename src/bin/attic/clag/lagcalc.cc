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

#include <cstdio>     
#include <cmath>      
#include <libciomr/libciomr.h>
#include <libqt/qt.h>         
#include <libpsio/psio.h>

extern "C" {
  extern FILE *infile;
  extern FILE *outfile;
}

namespace psi { namespace clag {

#define INDEX(x,y) ((x>y) ? ioff[x] + y : ioff[y] + x)
extern int *ioff;

/***************************************************************************/
/* The main function                                                       */
/***************************************************************************/

double lagcalc(double **OPDM, double *TPDM, double *h, double *TwoElec,
             double **lag, int nmo, int npop, int print_lvl, int lag_file)

{
  int i,j;                      /* indecies of lagrangian element          */
  int m,k,l;                    /* indecies of integrals needed            */ 
  int im,kl,imkl;               /* integral indecies combined              */
  int Tjm,Tmj,Tkl,Tjmkl,Tmjkl;  /* TPDM indices combined                   */
  double **oe_lag;              /* One-electron part of lagrangian         */
  double **te_lag;              /* Two-electron part of lagrangian         */ 
  double OEsum;                 /* QjmHim sumed over MO index m            */
  double TEsum;                 /* Gjmkl(im,kl) summed over m,k,l          */
  double lagtrace;              /* Trace of the Lagrangian (as a check)    */

  psio_open(lag_file, PSIO_OPEN_OLD);

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
                    Tjm = j * npop + m;
                    Tmj = m * npop + j;
                    Tkl = k * npop + l;
                    Tjmkl = INDEX(Tjm,Tkl);
                    Tmjkl = INDEX(Tmj,Tkl);
                    TEsum += (TPDM[Tmjkl] + TPDM[Tjmkl]) * TwoElec[imkl]; 
                  } /* end l loop */
            } /* end m loop */

          oe_lag[i][j] = OEsum;
          te_lag[i][j] = TEsum;
          lag[i][j] = OEsum + TEsum;
      
        } /* end j loop */ 

  /*
  ** check the trace of the Lagrangian (supposedly = energy)
  */
  for (lagtrace=0.0,i=0; i<nmo; i++)
    lagtrace += oe_lag[i][i] + 0.5 * te_lag[i][i];

  psio_write_entry(lag_file, "MO-basis Lagrangian", (char *) lag[0],
    nmo*nmo*sizeof(double));

  if (print_lvl > 1) {
    psi::fprintf(outfile,"\n\n One-electron part of the Lagrangian");
    print_mat(oe_lag, nmo, nmo, outfile);
    psi::fprintf(outfile,"\n\n Two-electron part of the Lagrangian");
    print_mat(te_lag, nmo, nmo, outfile);
    psi::fprintf(outfile,"\nLagrangian Matrix\n\n");
    print_mat(lag, nmo, nmo, outfile);
    }

  psio_close(lag_file, 1);

  free_block(oe_lag);
  free_block(te_lag);

  return(lagtrace);
}

}} // end namespace psi::clag

