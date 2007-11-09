/*! \file oe_deriv2_osrr.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
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

namespace psi { namespace CINTS {

void AI_Deriv2_OSrecurs(double ***AI0, double ***AIX, double ***AIY, double ***AIZ,
			double ***AIXX, double ***AIXY, double ***AIXZ,
			double ***AIYY, double ***AIYZ, double ***AIZZ,
			struct coordinates PA, struct coordinates PB,
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
  for(m=0;m<=mmax-2;m++) {
    AIXX[0][0][m] = 4*gamma*gamma*PC.x*PC.x*AI0[0][0][m+2] -
                    2*gamma*AI0[0][0][m+1];
    AIYY[0][0][m] = 4*gamma*gamma*PC.y*PC.y*AI0[0][0][m+2] -  
                    2*gamma*AI0[0][0][m+1];
    AIZZ[0][0][m] = 4*gamma*gamma*PC.z*PC.z*AI0[0][0][m+2] -  
                    2*gamma*AI0[0][0][m+1];
    AIXY[0][0][m] = 4*gamma*gamma*PC.x*PC.y*AI0[0][0][m+2];
    AIXZ[0][0][m] = 4*gamma*gamma*PC.x*PC.z*AI0[0][0][m+2];   
    AIYZ[0][0][m] = 4*gamma*gamma*PC.y*PC.z*AI0[0][0][m+2];   
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
        for(m=0;m<=mmax-b-2;m++) {	/* Gradients of the electric field */
          AIXX[0][jind][m] = PB.z*AIXX[0][jind-jzm][m] -
                             PC.z*AIXX[0][jind-jzm][m+1];
          AIYY[0][jind][m] = PB.z*AIYY[0][jind-jzm][m] -
                             PC.z*AIYY[0][jind-jzm][m+1];
          AIZZ[0][jind][m] = PB.z*AIZZ[0][jind-jzm][m] -
                             PC.z*AIZZ[0][jind-jzm][m+1] +
                               2*AIZ[0][jind-jzm][m+1];
          AIXY[0][jind][m] = PB.z*AIXY[0][jind-jzm][m] -
                             PC.z*AIXY[0][jind-jzm][m+1];
          AIXZ[0][jind][m] = PB.z*AIXZ[0][jind-jzm][m] -
                             PC.z*AIXZ[0][jind-jzm][m+1] +
                                 AIX[0][jind-jzm][m+1];
          AIYZ[0][jind][m] = PB.z*AIYZ[0][jind-jzm][m] -
                             PC.z*AIYZ[0][jind-jzm][m+1] +
                                 AIY[0][jind-jzm][m+1];
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
          for(m=0;m<=mmax-b-2;m++) {
            AIXX[0][jind][m] += pp*(jz-1)*(AIXX[0][jind-2*jzm][m] -
                                           AIXX[0][jind-2*jzm][m+1]);
            AIYY[0][jind][m] += pp*(jz-1)*(AIYY[0][jind-2*jzm][m] -
                                           AIYY[0][jind-2*jzm][m+1]);
            AIZZ[0][jind][m] += pp*(jz-1)*(AIZZ[0][jind-2*jzm][m] -             
                                           AIZZ[0][jind-2*jzm][m+1]);
            AIXY[0][jind][m] += pp*(jz-1)*(AIXY[0][jind-2*jzm][m] -             
                                           AIXY[0][jind-2*jzm][m+1]);
            AIXZ[0][jind][m] += pp*(jz-1)*(AIXZ[0][jind-2*jzm][m] -             
                                           AIXZ[0][jind-2*jzm][m+1]);
            AIYZ[0][jind][m] += pp*(jz-1)*(AIYZ[0][jind-2*jzm][m] -             
                                           AIYZ[0][jind-2*jzm][m+1]);
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
        for(m=0;m<=mmax-b-2;m++) {
          AIXX[0][jind][m] = PB.y*AIXX[0][jind-jym][m] -
                             PC.y*AIXX[0][jind-jym][m+1];
          AIYY[0][jind][m] = PB.y*AIYY[0][jind-jym][m] -
                             PC.y*AIYY[0][jind-jym][m+1] +
                               2*AIY[0][jind-jym][m+1];
          AIZZ[0][jind][m] = PB.y*AIZZ[0][jind-jym][m] -
                             PC.y*AIZZ[0][jind-jym][m+1];
          AIXY[0][jind][m] = PB.y*AIXY[0][jind-jym][m] -
                             PC.y*AIXY[0][jind-jym][m+1] +
                                 AIX[0][jind-jym][m+1];
          AIXZ[0][jind][m] = PB.y*AIXZ[0][jind-jym][m] -
                             PC.y*AIXZ[0][jind-jym][m+1];
          AIYZ[0][jind][m] = PB.y*AIYZ[0][jind-jym][m] -
                             PC.y*AIYZ[0][jind-jym][m+1] +
                                 AIZ[0][jind-jym][m+1];
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
          for(m=0;m<=mmax-b-2;m++) {
            AIXX[0][jind][m] += pp*(jy-1)*(AIXX[0][jind-2*jym][m] -
                                           AIXX[0][jind-2*jym][m+1]);
            AIYY[0][jind][m] += pp*(jy-1)*(AIYY[0][jind-2*jym][m] -
                                           AIYY[0][jind-2*jym][m+1]);
            AIZZ[0][jind][m] += pp*(jy-1)*(AIZZ[0][jind-2*jym][m] -             
                                           AIZZ[0][jind-2*jym][m+1]);
            AIXY[0][jind][m] += pp*(jy-1)*(AIXY[0][jind-2*jym][m] -             
                                           AIXY[0][jind-2*jym][m+1]);
            AIXZ[0][jind][m] += pp*(jy-1)*(AIXZ[0][jind-2*jym][m] -             
                                           AIXZ[0][jind-2*jym][m+1]);
            AIYZ[0][jind][m] += pp*(jy-1)*(AIYZ[0][jind-2*jym][m] -             
                                           AIYZ[0][jind-2*jym][m+1]);
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
        for(m=0;m<=mmax-b-2;m++) {
          AIXX[0][jind][m] = PB.x*AIXX[0][jind-jxm][m] -
                             PC.x*AIXX[0][jind-jxm][m+1] +
                               2*AIX[0][jind-jxm][m+1];
          AIYY[0][jind][m] = PB.x*AIYY[0][jind-jxm][m] -
                             PC.x*AIYY[0][jind-jxm][m+1];
          AIZZ[0][jind][m] = PB.x*AIZZ[0][jind-jxm][m] -
                             PC.x*AIZZ[0][jind-jxm][m+1];
          AIXY[0][jind][m] = PB.x*AIXY[0][jind-jxm][m] -
                             PC.x*AIXY[0][jind-jxm][m+1] +
                                 AIY[0][jind-jxm][m+1];
          AIXZ[0][jind][m] = PB.x*AIXZ[0][jind-jxm][m] -
                             PC.x*AIXZ[0][jind-jxm][m+1] +
                                 AIZ[0][jind-jxm][m+1];
          AIYZ[0][jind][m] = PB.x*AIYZ[0][jind-jxm][m] -
                             PC.x*AIYZ[0][jind-jxm][m+1];
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
          for(m=0;m<=mmax-b-2;m++) {
            AIXX[0][jind][m] += pp*(jx-1)*(AIXX[0][jind-2*jxm][m] -
                                           AIXX[0][jind-2*jxm][m+1]);
            AIYY[0][jind][m] += pp*(jx-1)*(AIYY[0][jind-2*jxm][m] -
                                           AIYY[0][jind-2*jxm][m+1]);
            AIZZ[0][jind][m] += pp*(jx-1)*(AIZZ[0][jind-2*jxm][m] -             
                                           AIZZ[0][jind-2*jxm][m+1]);
            AIXY[0][jind][m] += pp*(jx-1)*(AIXY[0][jind-2*jxm][m] -             
                                           AIXY[0][jind-2*jxm][m+1]);
            AIXZ[0][jind][m] += pp*(jx-1)*(AIXZ[0][jind-2*jxm][m] -             
                                           AIXZ[0][jind-2*jxm][m+1]);
            AIYZ[0][jind][m] += pp*(jx-1)*(AIYZ[0][jind-2*jxm][m] -             
                                           AIYZ[0][jind-2*jxm][m+1]);
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
          for(m=0;m<=mmax-a-b-2;m++) {	/* Gradients of the electric field */
            AIXX[iind][jind][m] = PA.z*AIXX[iind-izm][jind][m] -
                                  PC.z*AIXX[iind-izm][jind][m+1];
            AIYY[iind][jind][m] = PA.z*AIYY[iind-izm][jind][m] -
                                  PC.z*AIYY[iind-izm][jind][m+1];
            AIZZ[iind][jind][m] = PA.z*AIZZ[iind-izm][jind][m] -
                                  PC.z*AIZZ[iind-izm][jind][m+1] +
                                    2*AIZ[iind-izm][jind][m+1];
            AIXY[iind][jind][m] = PA.z*AIXY[iind-izm][jind][m] -
                                  PC.z*AIXY[iind-izm][jind][m+1];
            AIXZ[iind][jind][m] = PA.z*AIXZ[iind-izm][jind][m] -
                                  PC.z*AIXZ[iind-izm][jind][m+1] +
                                      AIX[iind-izm][jind][m+1];
            AIYZ[iind][jind][m] = PA.z*AIYZ[iind-izm][jind][m] -
                                  PC.z*AIYZ[iind-izm][jind][m+1] +
                                      AIY[iind-izm][jind][m+1];
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
            for(m=0;m<=mmax-a-b-2;m++) {
              AIXX[iind][jind][m] += pp*(iz-1)*(AIXX[iind-2*izm][jind][m] -
                                                AIXX[iind-2*izm][jind][m+1]);
              AIYY[iind][jind][m] += pp*(iz-1)*(AIYY[iind-2*izm][jind][m] -
                                                AIYY[iind-2*izm][jind][m+1]);
              AIZZ[iind][jind][m] += pp*(iz-1)*(AIZZ[iind-2*izm][jind][m] -             
                                                AIZZ[iind-2*izm][jind][m+1]);
              AIXY[iind][jind][m] += pp*(iz-1)*(AIXY[iind-2*izm][jind][m] -             
                                                AIXY[iind-2*izm][jind][m+1]);
              AIXZ[iind][jind][m] += pp*(iz-1)*(AIXZ[iind-2*izm][jind][m] -             
                                                AIXZ[iind-2*izm][jind][m+1]);
              AIYZ[iind][jind][m] += pp*(iz-1)*(AIYZ[iind-2*izm][jind][m] -             
                                                AIYZ[iind-2*izm][jind][m+1]);
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
            for(m=0;m<=mmax-a-b-2;m++) {
              AIXX[iind][jind][m] += pp*jz*(AIXX[iind-izm][jind-jzm][m] -
                                            AIXX[iind-izm][jind-jzm][m+1]);
              AIYY[iind][jind][m] += pp*jz*(AIYY[iind-izm][jind-jzm][m] -
                                            AIYY[iind-izm][jind-jzm][m+1]);
              AIZZ[iind][jind][m] += pp*jz*(AIZZ[iind-izm][jind-jzm][m] -             
                                            AIZZ[iind-izm][jind-jzm][m+1]);
              AIXY[iind][jind][m] += pp*jz*(AIXY[iind-izm][jind-jzm][m] -             
                                            AIXY[iind-izm][jind-jzm][m+1]);
              AIXZ[iind][jind][m] += pp*jz*(AIXZ[iind-izm][jind-jzm][m] -             
                                            AIXZ[iind-izm][jind-jzm][m+1]);
              AIYZ[iind][jind][m] += pp*jz*(AIYZ[iind-izm][jind-jzm][m] -             
                                            AIYZ[iind-izm][jind-jzm][m+1]);
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
          for(m=0;m<=mmax-a-b-2;m++) {
            AIXX[iind][jind][m] = PA.y*AIXX[iind-iym][jind][m] -
                                  PC.y*AIXX[iind-iym][jind][m+1];
            AIYY[iind][jind][m] = PA.y*AIYY[iind-iym][jind][m] -
                                  PC.y*AIYY[iind-iym][jind][m+1] +
                                    2*AIY[iind-iym][jind][m+1];
            AIZZ[iind][jind][m] = PA.y*AIZZ[iind-iym][jind][m] -
                                  PC.y*AIZZ[iind-iym][jind][m+1];
            AIXY[iind][jind][m] = PA.y*AIXY[iind-iym][jind][m] -
                                  PC.y*AIXY[iind-iym][jind][m+1] +
                                      AIX[iind-iym][jind][m+1];
            AIXZ[iind][jind][m] = PA.y*AIXZ[iind-iym][jind][m] -
                                  PC.y*AIXZ[iind-iym][jind][m+1];
            AIYZ[iind][jind][m] = PA.y*AIYZ[iind-iym][jind][m] -
                                  PC.y*AIYZ[iind-iym][jind][m+1] +
                                      AIZ[iind-iym][jind][m+1];
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
            for(m=0;m<=mmax-a-b-2;m++) {
              AIXX[iind][jind][m] += pp*(iy-1)*(AIXX[iind-2*iym][jind][m] -
                                                AIXX[iind-2*iym][jind][m+1]);
              AIYY[iind][jind][m] += pp*(iy-1)*(AIYY[iind-2*iym][jind][m] -
                                                AIYY[iind-2*iym][jind][m+1]);
              AIZZ[iind][jind][m] += pp*(iy-1)*(AIZZ[iind-2*iym][jind][m] -             
                                                AIZZ[iind-2*iym][jind][m+1]);
              AIXY[iind][jind][m] += pp*(iy-1)*(AIXY[iind-2*iym][jind][m] -             
                                                AIXY[iind-2*iym][jind][m+1]);
              AIXZ[iind][jind][m] += pp*(iy-1)*(AIXZ[iind-2*iym][jind][m] -             
                                                AIXZ[iind-2*iym][jind][m+1]);
              AIYZ[iind][jind][m] += pp*(iy-1)*(AIYZ[iind-2*iym][jind][m] -             
                                                AIYZ[iind-2*iym][jind][m+1]);
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
            for(m=0;m<=mmax-a-b-2;m++) {
              AIXX[iind][jind][m] += pp*jy*(AIXX[iind-iym][jind-jym][m] -
                                            AIXX[iind-iym][jind-jym][m+1]);
              AIYY[iind][jind][m] += pp*jy*(AIYY[iind-iym][jind-jym][m] -
                                            AIYY[iind-iym][jind-jym][m+1]);
              AIZZ[iind][jind][m] += pp*jy*(AIZZ[iind-iym][jind-jym][m] -             
                                            AIZZ[iind-iym][jind-jym][m+1]);
              AIXY[iind][jind][m] += pp*jy*(AIXY[iind-iym][jind-jym][m] -             
                                            AIXY[iind-iym][jind-jym][m+1]);
              AIXZ[iind][jind][m] += pp*jy*(AIXZ[iind-iym][jind-jym][m] -             
                                            AIXZ[iind-iym][jind-jym][m+1]);
              AIYZ[iind][jind][m] += pp*jy*(AIYZ[iind-iym][jind-jym][m] -             
                                            AIYZ[iind-iym][jind-jym][m+1]);
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
          for(m=0;m<=mmax-a-b-2;m++) {	/* Gradients of the electric field */
            AIXX[iind][jind][m] = PA.x*AIXX[iind-ixm][jind][m] -
                                  PC.x*AIXX[iind-ixm][jind][m+1] +
                                    2*AIX[iind-ixm][jind][m+1];
            AIYY[iind][jind][m] = PA.x*AIYY[iind-ixm][jind][m] -
                                  PC.x*AIYY[iind-ixm][jind][m+1];
            AIZZ[iind][jind][m] = PA.x*AIZZ[iind-ixm][jind][m] -
                                  PC.x*AIZZ[iind-ixm][jind][m+1];
            AIXY[iind][jind][m] = PA.x*AIXY[iind-ixm][jind][m] -
                                  PC.x*AIXY[iind-ixm][jind][m+1] +
                                      AIY[iind-ixm][jind][m+1];
            AIXZ[iind][jind][m] = PA.x*AIXZ[iind-ixm][jind][m] -
                                  PC.x*AIXZ[iind-ixm][jind][m+1] +
                                      AIZ[iind-ixm][jind][m+1];
            AIYZ[iind][jind][m] = PA.x*AIYZ[iind-ixm][jind][m] -
                                  PC.x*AIYZ[iind-ixm][jind][m+1];
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
            for(m=0;m<=mmax-a-b-2;m++) {
              AIXX[iind][jind][m] += pp*(ix-1)*(AIXX[iind-2*ixm][jind][m] -
                                                AIXX[iind-2*ixm][jind][m+1]);
              AIYY[iind][jind][m] += pp*(ix-1)*(AIYY[iind-2*ixm][jind][m] -
                                                AIYY[iind-2*ixm][jind][m+1]);
              AIZZ[iind][jind][m] += pp*(ix-1)*(AIZZ[iind-2*ixm][jind][m] -             
                                                AIZZ[iind-2*ixm][jind][m+1]);
              AIXY[iind][jind][m] += pp*(ix-1)*(AIXY[iind-2*ixm][jind][m] -             
                                                AIXY[iind-2*ixm][jind][m+1]);
              AIXZ[iind][jind][m] += pp*(ix-1)*(AIXZ[iind-2*ixm][jind][m] -             
                                                AIXZ[iind-2*ixm][jind][m+1]);
              AIYZ[iind][jind][m] += pp*(ix-1)*(AIYZ[iind-2*ixm][jind][m] -             
                                                AIYZ[iind-2*ixm][jind][m+1]);
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
            for(m=0;m<=mmax-a-b-2;m++) {
              AIXX[iind][jind][m] += pp*jx*(AIXX[iind-ixm][jind-jxm][m] -
                                            AIXX[iind-ixm][jind-jxm][m+1]);
              AIYY[iind][jind][m] += pp*jx*(AIYY[iind-ixm][jind-jxm][m] -
                                            AIYY[iind-ixm][jind-jxm][m+1]);
              AIZZ[iind][jind][m] += pp*jx*(AIZZ[iind-ixm][jind-jxm][m] -             
                                            AIZZ[iind-ixm][jind-jxm][m+1]);
              AIXY[iind][jind][m] += pp*jx*(AIXY[iind-ixm][jind-jxm][m] -             
                                            AIXY[iind-ixm][jind-jxm][m+1]);
              AIXZ[iind][jind][m] += pp*jx*(AIXZ[iind-ixm][jind-jxm][m] -             
                                            AIXZ[iind-ixm][jind-jxm][m+1]);
              AIYZ[iind][jind][m] += pp*jx*(AIYZ[iind-ixm][jind-jxm][m] -             
                                            AIYZ[iind-ixm][jind-jxm][m+1]);
            }
          }
        }
        else  /*--- Should never happen ---*/
	  abort();
      }
    }

  return;
}
};};
