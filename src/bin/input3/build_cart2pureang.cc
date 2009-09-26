/*! \file
 \ingroup INPUT
 \brief Enter brief description of file here 
 */
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <psiconfig.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace {
  double xyz2lm_Coeff(int, int, int, int, int);
  double **bc;
  double *fac;
}

namespace psi {
  namespace input {
    
    /*-----------------------------------------------------------------------------------------------------------------
     This function builds cartesian to pure angular momentum transformation matrix
     -----------------------------------------------------------------------------------------------------------------*/

    void build_cart2pureang() {
      int i, j, k, m;
      int n;
      int lmax, irr, ao, atom, l;
      int shell, shell_first, shell_last;
      int ao_max;
      int cnt, so_cnt;
      int ao_offset;
      int n_max;
      double min, tmp;
      
      /*------------------------
       Initialize global arrays
       -------------------------*/

      num_so_per_irrep = init_int_array(nirreps);
      cart2pureang = (double ***) malloc(sizeof(double **)*(max_angmom+1));
      for (l=0; l<=max_angmom; l++)
        cart2pureang[l] = init_matrix(2*l+1, ioff[l+1]);
      fac = init_array(2*max_angmom+1);
      fac[0] = 1.0;
      for (i=1; i<=2*max_angmom; i++)
        fac[i] = i*fac[i-1];
      bc = init_matrix(max_angmom+1, max_angmom+1);
      for (i=0; i<=max_angmom; i++)
        for (j=0; j<=i; j++)
          bc[i][j] = combinations(i, j);
      
      /*------------------------
       Initialize local arrays
       -------------------------*/

      /*--------------------------------------------------------------------------------------------------
       Compute number of pure angular momentum SOs in each symmetry block for each angular momentum type
       --------------------------------------------------------------------------------------------------*/

      num_pureang_so = init_int_matrix(max_angmom+1, nirreps);
      num_redun_so = init_int_matrix(max_angmom+1, nirreps);
      for (l=0; l<=max_angmom; l++)
        for (irr=0; irr<nirreps; irr++)
          num_pureang_so[l][irr] = num_cart_so[l][irr];
      for (i=2; i<=max_angmom; i++)
        for (l=i%2; l<i; l+=2)
          for (irr=0; irr<nirreps; irr++) {
            num_redun_so[i][irr] += num_pureang_so[l][irr];
            num_pureang_so[i][irr] -= num_pureang_so[l][irr];
          }
      
      /*--------------------------------------------------------------------------------------------------
       Compute cartesian to pure angular momentum transformation matrices for each angular momentum type
       --------------------------------------------------------------------------------------------------*/

      for (l=0; l<=max_angmom; l++) {
        ao_max = ioff[l+1];
        for (m=-l; m<=l; m++)
          for (ao=0; ao<ao_max; ao++)
            cart2pureang[l][m+l][ao] = xyz2lm_Coeff(l, m, xexp_ao[l][ao],
                                                    yexp_ao[l][ao],
                                                    zexp_ao[l][ao]);
      }
      
      /*-----------------------------
       Remove after testing is over
       -----------------------------*/
      if (print_lvl >= DEBUGPRINT)
        for (l=0; l<=max_angmom; l++) {
          fprintf(outfile, "    -Cart2PureAng matrix for l=%d:\n", l);
          print_mat(cart2pureang[l], 2*l+1, ioff[l+1], outfile);
          fprintf(outfile, "\n");
        }
      
      free(fac);
      free_matrix(bc, max_angmom+1);
      return;
    }
  
  }
} // namespace psi::input

namespace {
  /*---------------------------------------------------------------------------------------------
   Computes transformation coefficients from cartesian to real pure angular momentum functions.
   See IJQC 54, 83 (1995), eqn (15). If m is negative, imaginary part is computed, whereas
   a positive m indicates that the real part of spherical harmonic Ylm is requested.
   ---------------------------------------------------------------------------------------------*/
  double xyz2lm_Coeff(int l, int m, int lx, int ly, int lz) {
    using namespace psi::input;
    static int use_cca_integrals_standard = (PSI_INTEGRALS_STANDARD == 1);
    int i, j, k, i_max;
    int k_min, k_max;
    int abs_m;
    int comp;
    double pfac, pfac1, sum, sum1;
    
    abs_m = abs(m);
    if ((lx + ly - abs(m))%2)
      return 0.0;
    else
      j = (lx + ly - abs(m))/2;
    
    if (j < 0)
      return 0.0;
    
    /*----------------------------------------------------------------------------------------
     Checking whether the cartesian polynomial contributes to the requested component of Ylm
     ----------------------------------------------------------------------------------------*/
    comp = (m >= 0) ? 1 : -1;
    i = abs_m-lx;
    if (comp != parity(i))
      return 0.0;
    
    pfac = sqrt(fac[2*lx]*fac[2*ly]*fac[2*lz]*fac[l-abs_m]/(fac[2*l]*fac[l]
        *fac[lx]*fac[ly]*fac[lz]*fac[l+abs_m]));
    pfac /= (1 << l);
    
    if (m < 0)
      pfac *= parity((i-1)/2);
    else
      pfac *= parity(i/2);
    
    i_max = (l-abs_m)/2;
    sum = 0.0;
    for (i=0; i<=i_max; i++) {
      pfac1 = bc[l][i]*bc[i][j];
      if (pfac1 == 0.0)
        continue;
      else
        pfac1 *= (parity(i)*fac[2*(l-i)]/fac[l-abs_m-2*i]);
      sum1 = 0.0;
      k_min = MAX((lx-abs_m)/2,0);
      k_max = MIN(j,lx/2);
      for (k=k_min; k<=k_max; k++)
        sum1 += bc[j][k]*bc[abs_m][lx-2*k]*parity(k);
      sum += pfac1*sum1;
    }
    
    if (use_cca_integrals_standard)
      sum *= sqrt(df[2*l]/(df[2*lx]*df[2*ly]*df[2*lz]));
    
    if (m == 0)
      return pfac*sum;
    else
      return M_SQRT2*pfac*sum;
  }
} // namespace
