/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"small_fns.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"fjt.h"
#endif


namespace psi { namespace CINTS {
/*! Recurrence relation are from the same paper - pp. 3971-3972 */

void AI_OSrecurs(double ***AI0, struct coordinates PA, struct coordinates PB,
		 struct coordinates PC, double gamma, int iang, int jang)
{
  int a,b,m;
  int izm = 1;
  int iym = iang + 1;
  int ixm = iym * iym;
  int jzm = 1;
  int jym = jang + 1;
  int jxm = jym * jym;
  int ix,iy,iz,jx,jy,jz;
  int iind,jind;
  double pp = 1/(2*gamma);
  int mmax = iang+jang;
  double tmp = sqrt(gamma)*M_2_SQRTPI;
  double u = gamma*(PC.x*PC.x + PC.y*PC.y + PC.z*PC.z);
  static double F[2*CINTS_MAX_AM+1];

#ifdef USE_TAYLOR_FM
  taylor_compute_fm(F,u,mmax);
#else
  calc_f(F,mmax,u);
#endif


	/* Computing starting integrals for recursion */

  for(m=0;m<=mmax;m++)
    AI0[0][0][m] = tmp*F[m];

	/* Upward recursion in j with i=0 */
  
  for(b=1;b<=jang;b++)
    for(jx=0;jx<=b;jx++)
    for(jy=0;jy<=b-jx;jy++) {
      jz = b-jx-jy;
      jind = jx*jxm+jy*jym+jz*jzm;
      if (jz > 0) {
        for(m=0;m<=mmax-b;m++)	/* Electrostatic potential integrals */
          AI0[0][jind][m] = PB.z*AI0[0][jind-jzm][m] - 
                            PC.z*AI0[0][jind-jzm][m+1];
	if (jz > 1) {
          for(m=0;m<=mmax-b;m++)
            AI0[0][jind][m] += pp*(jz-1)*(AI0[0][jind-2*jzm][m] -
                                          AI0[0][jind-2*jzm][m+1]);
        }
      }
      else 
      if (jy > 0) {
        for(m=0;m<=mmax-b;m++)
          AI0[0][jind][m] = PB.y*AI0[0][jind-jym][m] -
                            PC.y*AI0[0][jind-jym][m+1];
        if (jy > 1) {
          for(m=0;m<=mmax-b;m++)  
            AI0[0][jind][m] += pp*(jy-1)*(AI0[0][jind-2*jym][m] -
                                          AI0[0][jind-2*jym][m+1]);
        }
      }
      else
      if (jx > 0) {
        for(m=0;m<=mmax-b;m++)
          AI0[0][jind][m] = PB.x*AI0[0][jind-jxm][m] -
                            PC.x*AI0[0][jind-jxm][m+1];
        if (jx > 1) {
          for(m=0;m<=mmax-b;m++)  
            AI0[0][jind][m] += pp*(jx-1)*(AI0[0][jind-2*jxm][m] -
                                          AI0[0][jind-2*jxm][m+1]);
        }
      }
      else
        throw std::domain_error("  There's some error in the AI_OSrecurs algorithm\n\n");
    }
 


  /* The following fragment cannot be vectorized easily, I guess :-) */
	/* Upward recursion in i with all possible j's */

  for(b=0;b<=jang;b++)
    for(jx=0;jx<=b;jx++)
    for(jy=0;jy<=b-jx;jy++) {
    jz = b-jx-jy;
    jind = jx*jxm + jy*jym + jz*jzm;
    for(a=1;a<=iang;a++)
      for(ix=0;ix<=a;ix++)
      for(iy=0;iy<=a-ix;iy++) {
        iz = a-ix-iy;
        iind = ix*ixm + iy*iym + iz*izm;
        if (iz > 0) {
          for(m=0;m<=mmax-a-b;m++)
            AI0[iind][jind][m] = PA.z*AI0[iind-izm][jind][m] - 
                                 PC.z*AI0[iind-izm][jind][m+1];
          if (iz > 1) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*(iz-1)*
               (AI0[iind-2*izm][jind][m] - AI0[iind-2*izm][jind][m+1]);
          }
          if (jz > 0) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*jz*
               (AI0[iind-izm][jind-jzm][m] - AI0[iind-izm][jind-jzm][m+1]);
          }
        }
        else
	if (iy > 0) {
          for(m=0;m<=mmax-a-b;m++)
            AI0[iind][jind][m] = PA.y*AI0[iind-iym][jind][m] -
                                 PC.y*AI0[iind-iym][jind][m+1];
	  if (iy > 1) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*(iy-1)*
              (AI0[iind-2*iym][jind][m] - AI0[iind-2*iym][jind][m+1]);
          }
	  if (jy > 0) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*jy*
               (AI0[iind-iym][jind-jym][m] - AI0[iind-iym][jind-jym][m+1]);
          }
        }
        else
	if (ix > 0) {
          for(m=0;m<=mmax-a-b;m++)
            AI0[iind][jind][m] = PA.x*AI0[iind-ixm][jind][m] -
                                 PC.x*AI0[iind-ixm][jind][m+1];
          if (ix > 1) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*(ix-1)*
               (AI0[iind-2*ixm][jind][m] - AI0[iind-2*ixm][jind][m+1]);
          }
          if (jx > 0) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*jx*
               (AI0[iind-ixm][jind-jxm][m] - AI0[iind-ixm][jind-jxm][m+1]);  
          }
        }
        else
          throw std::domain_error("  There's some error in the AI_OSrecurs algorithm\n\n");
      }
    }

  return;
}


void OI_OSrecurs(double **OIX, double **OIY, double **OIZ, struct coordinates PA, struct coordinates PB,
		 double gamma, int lmaxi, int lmaxj)
{
  int i,j,k;
  double pp = 1/(2*gamma);

  OIX[0][0] = OIY[0][0] = OIZ[0][0] = 1.0;

	/* Upward recursion in j for i=0 */

  OIX[0][1] = PB.x;
  OIY[0][1] = PB.y;
  OIZ[0][1] = PB.z;

  for(j=1;j<lmaxj;j++) {
    OIX[0][j+1] = PB.x*OIX[0][j];
    OIY[0][j+1] = PB.y*OIY[0][j];
    OIZ[0][j+1] = PB.z*OIZ[0][j];
    OIX[0][j+1] += j*pp*OIX[0][j-1];
    OIY[0][j+1] += j*pp*OIY[0][j-1];
    OIZ[0][j+1] += j*pp*OIZ[0][j-1];
  }

	/* Upward recursion in i for all j's */

  OIX[1][0] = PA.x;
  OIY[1][0] = PA.y;
  OIZ[1][0] = PA.z;
  for(j=1;j<=lmaxj;j++) {
    OIX[1][j] = PA.x*OIX[0][j];
    OIY[1][j] = PA.y*OIY[0][j];
    OIZ[1][j] = PA.z*OIZ[0][j];
    OIX[1][j] += j*pp*OIX[0][j-1];
    OIY[1][j] += j*pp*OIY[0][j-1];
    OIZ[1][j] += j*pp*OIZ[0][j-1];
  }
  for(i=1;i<lmaxi;i++) {
    OIX[i+1][0] = PA.x*OIX[i][0];
    OIY[i+1][0] = PA.y*OIY[i][0];
    OIZ[i+1][0] = PA.z*OIZ[i][0];
    OIX[i+1][0] += i*pp*OIX[i-1][0];
    OIY[i+1][0] += i*pp*OIY[i-1][0];
    OIZ[i+1][0] += i*pp*OIZ[i-1][0];
    for(j=1;j<=lmaxj;j++) {
      OIX[i+1][j] = PA.x*OIX[i][j];
      OIY[i+1][j] = PA.y*OIY[i][j];
      OIZ[i+1][j] = PA.z*OIZ[i][j];
      OIX[i+1][j] += i*pp*OIX[i-1][j];
      OIY[i+1][j] += i*pp*OIY[i-1][j];
      OIZ[i+1][j] += i*pp*OIZ[i-1][j];
      OIX[i+1][j] += j*pp*OIX[i][j-1];
      OIY[i+1][j] += j*pp*OIY[i][j-1];
      OIZ[i+1][j] += j*pp*OIZ[i][j-1];
    }
  }

  return;
}
};};
