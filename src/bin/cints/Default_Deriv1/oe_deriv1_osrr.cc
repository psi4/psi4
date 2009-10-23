/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libint/libint.h>
#include <libderiv/libderiv.h>
#include <psifiles.h>
#include "defines.h"
#define EXTERN
#include "global.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include "fjt.h"
#endif

namespace psi { namespace cints {
void AI_Deriv1_OSrecurs(double ***AI0, double ***AIX, double ***AIY, double ***AIZ, struct coordinates PA, struct coordinates PB,
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
  for(m=0;m<=mmax-1;m++) {
    AIX[0][0][m] = 2*gamma*PC.x*AI0[0][0][m+1];
    AIY[0][0][m] = 2*gamma*PC.y*AI0[0][0][m+1];
    AIZ[0][0][m] = 2*gamma*PC.z*AI0[0][0][m+1];
  }

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
	for(m=0;m<=mmax-b-1;m++) {	/* Electric field integrals */
          AIX[0][jind][m] = PB.z*AIX[0][jind-jzm][m] -
                            PC.z*AIX[0][jind-jzm][m+1];
          AIY[0][jind][m] = PB.z*AIY[0][jind-jzm][m] -
                            PC.z*AIY[0][jind-jzm][m+1];
          AIZ[0][jind][m] = PB.z*AIZ[0][jind-jzm][m] -
                            PC.z*AIZ[0][jind-jzm][m+1] +
                                AI0[0][jind-jzm][m+1];
        }
	if (jz > 1) {
          for(m=0;m<=mmax-b;m++)
            AI0[0][jind][m] += pp*(jz-1)*(AI0[0][jind-2*jzm][m] -
                                          AI0[0][jind-2*jzm][m+1]);
	  for(m=0;m<=mmax-b-1;m++) {
            AIX[0][jind][m] += pp*(jz-1)*(AIX[0][jind-2*jzm][m] -
                                          AIX[0][jind-2*jzm][m+1]);
            AIY[0][jind][m] += pp*(jz-1)*(AIY[0][jind-2*jzm][m] -
                                          AIY[0][jind-2*jzm][m+1]);
            AIZ[0][jind][m] += pp*(jz-1)*(AIZ[0][jind-2*jzm][m] -
                                          AIZ[0][jind-2*jzm][m+1]);
          }
	}
      }
      else 
      if (jy > 0) {
        for(m=0;m<=mmax-b;m++)
          AI0[0][jind][m] = PB.y*AI0[0][jind-jym][m] -
                            PC.y*AI0[0][jind-jym][m+1];
	for(m=0;m<=mmax-b-1;m++) {
          AIX[0][jind][m] = PB.y*AIX[0][jind-jym][m] -
                            PC.y*AIX[0][jind-jym][m+1];
          AIY[0][jind][m] = PB.y*AIY[0][jind-jym][m] -
                            PC.y*AIY[0][jind-jym][m+1] +
                                AI0[0][jind-jym][m+1];
          AIZ[0][jind][m] = PB.y*AIZ[0][jind-jym][m] -
                            PC.y*AIZ[0][jind-jym][m+1];
        }
        if (jy > 1) {
          for(m=0;m<=mmax-b;m++)  
            AI0[0][jind][m] += pp*(jy-1)*(AI0[0][jind-2*jym][m] -
                                          AI0[0][jind-2*jym][m+1]);
	  for(m=0;m<=mmax-b-1;m++) {
            AIX[0][jind][m] += pp*(jy-1)*(AIX[0][jind-2*jym][m] -
                                          AIX[0][jind-2*jym][m+1]);
            AIY[0][jind][m] += pp*(jy-1)*(AIY[0][jind-2*jym][m] -
                                          AIY[0][jind-2*jym][m+1]);
            AIZ[0][jind][m] += pp*(jy-1)*(AIZ[0][jind-2*jym][m] -
                                          AIZ[0][jind-2*jym][m+1]);
          }
        }
      }
      else
      if (jx > 0) {
        for(m=0;m<=mmax-b;m++)
          AI0[0][jind][m] = PB.x*AI0[0][jind-jxm][m] -
                            PC.x*AI0[0][jind-jxm][m+1];
	for(m=0;m<=mmax-b-1;m++) {
          AIX[0][jind][m] = PB.x*AIX[0][jind-jxm][m] -
                            PC.x*AIX[0][jind-jxm][m+1] +
                                AI0[0][jind-jxm][m+1];
          AIY[0][jind][m] = PB.x*AIY[0][jind-jxm][m] -
                            PC.x*AIY[0][jind-jxm][m+1];
          AIZ[0][jind][m] = PB.x*AIZ[0][jind-jxm][m] -
                            PC.x*AIZ[0][jind-jxm][m+1];
        }
        if (jx > 1) {
          for(m=0;m<=mmax-b;m++)  
            AI0[0][jind][m] += pp*(jx-1)*(AI0[0][jind-2*jxm][m] -
                                          AI0[0][jind-2*jxm][m+1]);
	  for(m=0;m<=mmax-b-1;m++) {
            AIX[0][jind][m] += pp*(jx-1)*(AIX[0][jind-2*jxm][m] -
                                          AIX[0][jind-2*jxm][m+1]);
            AIY[0][jind][m] += pp*(jx-1)*(AIY[0][jind-2*jxm][m] -
                                          AIY[0][jind-2*jxm][m+1]);
            AIZ[0][jind][m] += pp*(jx-1)*(AIZ[0][jind-2*jxm][m] -
                                          AIZ[0][jind-2*jxm][m+1]);
          }
        }
      }
      else /*--- Should never happen ---*/
	abort();
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
	  for(m=0;m<=mmax-a-b-1;m++) {	/* Electric field integrals */
            AIX[iind][jind][m] = PA.z*AIX[iind-izm][jind][m] -
                                 PC.z*AIX[iind-izm][jind][m+1];
            AIY[iind][jind][m] = PA.z*AIY[iind-izm][jind][m] -
                                 PC.z*AIY[iind-izm][jind][m+1];
            AIZ[iind][jind][m] = PA.z*AIZ[iind-izm][jind][m] -
                                 PC.z*AIZ[iind-izm][jind][m+1] +
                                     AI0[iind-izm][jind][m+1];
          }
          if (iz > 1) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*(iz-1)*
               (AI0[iind-2*izm][jind][m] - AI0[iind-2*izm][jind][m+1]);
	    for(m=0;m<=mmax-a-b-1;m++) {
              AIX[iind][jind][m] += pp*(iz-1)*(AIX[iind-2*izm][jind][m] -
                                               AIX[iind-2*izm][jind][m+1]);
              AIY[iind][jind][m] += pp*(iz-1)*(AIY[iind-2*izm][jind][m] -
                                               AIY[iind-2*izm][jind][m+1]);
              AIZ[iind][jind][m] += pp*(iz-1)*(AIZ[iind-2*izm][jind][m] -
                                               AIZ[iind-2*izm][jind][m+1]);
            }
          }
          if (jz > 0) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*jz*
               (AI0[iind-izm][jind-jzm][m] - AI0[iind-izm][jind-jzm][m+1]);
	    for(m=0;m<=mmax-a-b-1;m++) {
              AIX[iind][jind][m] += pp*jz*(AIX[iind-izm][jind-jzm][m] -
                                           AIX[iind-izm][jind-jzm][m+1]);
              AIY[iind][jind][m] += pp*jz*(AIY[iind-izm][jind-jzm][m] -
                                           AIY[iind-izm][jind-jzm][m+1]);
              AIZ[iind][jind][m] += pp*jz*(AIZ[iind-izm][jind-jzm][m] -
                                           AIZ[iind-izm][jind-jzm][m+1]);
            }
          }
        }
        else
	if (iy > 0) {
          for(m=0;m<=mmax-a-b;m++)
            AI0[iind][jind][m] = PA.y*AI0[iind-iym][jind][m] -
                                 PC.y*AI0[iind-iym][jind][m+1];
	  for(m=0;m<=mmax-a-b-1;m++) {
            AIX[iind][jind][m] = PA.y*AIX[iind-iym][jind][m] -
                                 PC.y*AIX[iind-iym][jind][m+1];
            AIY[iind][jind][m] = PA.y*AIY[iind-iym][jind][m] -
                                 PC.y*AIY[iind-iym][jind][m+1] +
                                     AI0[iind-iym][jind][m+1];
            AIZ[iind][jind][m] = PA.y*AIZ[iind-iym][jind][m] -
                                 PC.y*AIZ[iind-iym][jind][m+1];
          }
	  if (iy > 1) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*(iy-1)*
              (AI0[iind-2*iym][jind][m] - AI0[iind-2*iym][jind][m+1]);
	    for(m=0;m<=mmax-a-b-1;m++) {
              AIX[iind][jind][m] += pp*(iy-1)*(AIX[iind-2*iym][jind][m] -
                                               AIX[iind-2*iym][jind][m+1]);
              AIY[iind][jind][m] += pp*(iy-1)*(AIY[iind-2*iym][jind][m] -
                                               AIY[iind-2*iym][jind][m+1]);
              AIZ[iind][jind][m] += pp*(iy-1)*(AIZ[iind-2*iym][jind][m] -
                                               AIZ[iind-2*iym][jind][m+1]);
            }
          }
	  if (jy > 0) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*jy*
               (AI0[iind-iym][jind-jym][m] - AI0[iind-iym][jind-jym][m+1]);
	    for(m=0;m<=mmax-a-b-1;m++) {
              AIX[iind][jind][m] += pp*jy*(AIX[iind-iym][jind-jym][m] -
                                           AIX[iind-iym][jind-jym][m+1]);
              AIY[iind][jind][m] += pp*jy*(AIY[iind-iym][jind-jym][m] -
                                           AIY[iind-iym][jind-jym][m+1]);
              AIZ[iind][jind][m] += pp*jy*(AIZ[iind-iym][jind-jym][m] -
                                           AIZ[iind-iym][jind-jym][m+1]);
            }
          }
        }
        else
	if (ix > 0) {
          for(m=0;m<=mmax-a-b;m++)
            AI0[iind][jind][m] = PA.x*AI0[iind-ixm][jind][m] -
                                 PC.x*AI0[iind-ixm][jind][m+1];
	  for(m=0;m<=mmax-a-b-1;m++) {	/* Electric field integrals */
            AIX[iind][jind][m] = PA.x*AIX[iind-ixm][jind][m] -
                                 PC.x*AIX[iind-ixm][jind][m+1] +
                                     AI0[iind-ixm][jind][m+1];
            AIY[iind][jind][m] = PA.x*AIY[iind-ixm][jind][m] -
                                 PC.x*AIY[iind-ixm][jind][m+1];
            AIZ[iind][jind][m] = PA.x*AIZ[iind-ixm][jind][m] -
                                 PC.x*AIZ[iind-ixm][jind][m+1];
          }
          if (ix > 1) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*(ix-1)*
               (AI0[iind-2*ixm][jind][m] - AI0[iind-2*ixm][jind][m+1]);
	    for(m=0;m<=mmax-a-b-1;m++) {
              AIX[iind][jind][m] += pp*(ix-1)*(AIX[iind-2*ixm][jind][m] -
                                               AIX[iind-2*ixm][jind][m+1]);
              AIY[iind][jind][m] += pp*(ix-1)*(AIY[iind-2*ixm][jind][m] -
                                               AIY[iind-2*ixm][jind][m+1]);
              AIZ[iind][jind][m] += pp*(ix-1)*(AIZ[iind-2*ixm][jind][m] -
                                               AIZ[iind-2*ixm][jind][m+1]);
            }
          }
          if (jx > 0) {
            for(m=0;m<=mmax-a-b;m++)
              AI0[iind][jind][m] += pp*jx*
               (AI0[iind-ixm][jind-jxm][m] - AI0[iind-ixm][jind-jxm][m+1]);
	    for(m=0;m<=mmax-a-b-1;m++) {
              AIX[iind][jind][m] += pp*jx*(AIX[iind-ixm][jind-jxm][m] -
                                           AIX[iind-ixm][jind-jxm][m+1]);
              AIY[iind][jind][m] += pp*jx*(AIY[iind-ixm][jind-jxm][m] -
                                           AIY[iind-ixm][jind-jxm][m+1]);
              AIZ[iind][jind][m] += pp*jx*(AIZ[iind-ixm][jind-jxm][m] -
                                           AIZ[iind-ixm][jind-jxm][m+1]);
            }
          }
        }
        else  /*--- Should never happen ---*/
	  abort();
      }
    }

  return;
}
}}
