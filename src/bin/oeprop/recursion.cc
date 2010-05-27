/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "globals.h"
#include "prototypes.h"

namespace {
  void calc_f(double *, int, double);
}

namespace psi { namespace oeprop {

	/* Recursion relations are taken from Obara and Saika paper
	   JCP 84, 3963, 1986. */

void MI_OSrecurs(double pax, double pay, double paz, 
                 double pbx, double pby, double pbz, double gamma,
                 int lmaxi, int lmaxj, int maxm)
{
  int i,j,k;
  double pp = 1/(2*gamma);

	/* Computing starting integrals for recursive procedure */

  if (maxm > 1)
    MIX[0][0][2] = MIY[0][0][2] = MIZ[0][0][2] = pp;

	/* Upward recursion in j for i=0 */

  for(j=0;j<lmaxj;j++)
    for(k=0;k<=maxm;k++) {
      MIX[0][j+1][k] = pbx*MIX[0][j][k];
      MIY[0][j+1][k] = pby*MIY[0][j][k];
      MIZ[0][j+1][k] = pbz*MIZ[0][j][k];
      if (j>0) {
        MIX[0][j+1][k] += j*pp*MIX[0][j-1][k];
        MIY[0][j+1][k] += j*pp*MIY[0][j-1][k];
        MIZ[0][j+1][k] += j*pp*MIZ[0][j-1][k];
      }
      if (k>0) {
        MIX[0][j+1][k] += k*pp*MIX[0][j][k-1];
        MIY[0][j+1][k] += k*pp*MIY[0][j][k-1];
        MIZ[0][j+1][k] += k*pp*MIZ[0][j][k-1];
      }
    }

	/* Upward recursion in i for all j's */

  for(i=0;i<lmaxi;i++)
    for(j=0;j<=lmaxj;j++)
      for(k=0;k<=maxm;k++) {
        MIX[i+1][j][k] = pax*MIX[i][j][k];
        MIY[i+1][j][k] = pay*MIY[i][j][k];
        MIZ[i+1][j][k] = paz*MIZ[i][j][k];
        if (i>0) {
          MIX[i+1][j][k] += i*pp*MIX[i-1][j][k];
          MIY[i+1][j][k] += i*pp*MIY[i-1][j][k];
          MIZ[i+1][j][k] += i*pp*MIZ[i-1][j][k];
        }
        if (j>0) {
          MIX[i+1][j][k] += j*pp*MIX[i][j-1][k];
          MIY[i+1][j][k] += j*pp*MIY[i][j-1][k];
          MIZ[i+1][j][k] += j*pp*MIZ[i][j-1][k];
        }
        if (k>0) {
          MIX[i+1][j][k] += k*pp*MIX[i][j][k-1];
          MIY[i+1][j][k] += k*pp*MIY[i][j][k-1];
          MIZ[i+1][j][k] += k*pp*MIZ[i][j][k-1];
        }
      }
}


/* Recurrence relation are from the same paper - pp. 3971-3972 */

void AI_OSrecurs(double pax, double pay, double paz,
                 double pbx, double pby, double pbz, 
                 double pcx, double pcy, double pcz, 
                 double gamma, int iang, int jang)
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
  double mmax = 2*lmax+2;
  double tmp = sqrt(gamma)*M_2_SQRTPI;
  double u = gamma*(pcx*pcx + pcy*pcy + pcz*pcz);
  double *F;
  
  F = init_array((int) (mmax+1));
  calc_f(F,(int) mmax,u);

	/* Computing starting integrals for recursion */
  fprintf(outfile, "setup:\n");
  for(m=0;m<=mmax;m++)
    AI0[0][0][m] = tmp*F[m];
  for(m=0;m<=mmax-1;m++) {
    AIX[0][0][m] = 2*gamma*pcx*AI0[0][0][m+1];
    AIY[0][0][m] = 2*gamma*pcy*AI0[0][0][m+1];
    AIZ[0][0][m] = 2*gamma*pcz*AI0[0][0][m+1];

    fprintf(outfile, "AIX[0][0][%d] = %lf\tAIY[0][0][%d] = %lf\tAIZ[0][0][%d] = %lf\n", m, AIX[0][0][m], m, AIY[0][0][m], m, AIZ[0][0][m]);
  }
  for(m=0;m<=mmax-2;m++) {
    AIXX[0][0][m] = 4*gamma*gamma*pcx*pcx*AI0[0][0][m+2] -
                    2*gamma*AI0[0][0][m+1];
    AIYY[0][0][m] = 4*gamma*gamma*pcy*pcy*AI0[0][0][m+2] -  
                    2*gamma*AI0[0][0][m+1];
    AIZZ[0][0][m] = 4*gamma*gamma*pcz*pcz*AI0[0][0][m+2] -  
                    2*gamma*AI0[0][0][m+1];
    AIXY[0][0][m] = 4*gamma*gamma*pcx*pcy*AI0[0][0][m+2];
    AIXZ[0][0][m] = 4*gamma*gamma*pcx*pcz*AI0[0][0][m+2];   
    AIYZ[0][0][m] = 4*gamma*gamma*pcy*pcz*AI0[0][0][m+2];   
  }
  

	/* Upward recursion in j with i=0 */
  fprintf(outfile, "upward in j with i=0:\n");
  for(b=1;b<=jang;b++)
    for(jx=0;jx<=b;jx++)
    for(jy=0;jy<=b-jx;jy++) {
      jz = b-jx-jy;
      jind = jx*jxm+jy*jym+jz*jzm;
      if (jz > 0) {
        for(m=0;m<=mmax-b;m++)	/* Electrostatic potential integrals */
          AI0[0][jind][m] = pbz*AI0[0][jind-jzm][m] - 
                            pcz*AI0[0][jind-jzm][m+1];
        for(m=0;m<=mmax-b-1;m++) {	/* Electric field integrals */
          AIX[0][jind][m] = pbz*AIX[0][jind-jzm][m] -
                            pcz*AIX[0][jind-jzm][m+1];
          AIY[0][jind][m] = pbz*AIY[0][jind-jzm][m] -
                            pcz*AIY[0][jind-jzm][m+1];
          AIZ[0][jind][m] = pbz*AIZ[0][jind-jzm][m] -
                            pcz*AIZ[0][jind-jzm][m+1] +
                                AI0[0][jind-jzm][m+1];
          fprintf(outfile, "AIX[0][%d][%d] = %lf\tAIY[0][%d][%d] = %lf\tAIZ[0][%d][%d] = %lf\n", jind, m, AIX[0][jind][m], jind, m, AIY[0][jind][m], jind, m, AIZ[0][jind][m]);
        }
        for(m=0;m<=mmax-b-2;m++) {	/* Gradients of the electric field */
          AIXX[0][jind][m] = pbz*AIXX[0][jind-jzm][m] -
                             pcz*AIXX[0][jind-jzm][m+1];
          AIYY[0][jind][m] = pbz*AIYY[0][jind-jzm][m] -
                             pcz*AIYY[0][jind-jzm][m+1];
          AIZZ[0][jind][m] = pbz*AIZZ[0][jind-jzm][m] -
                             pcz*AIZZ[0][jind-jzm][m+1] +
                               2*AIZ[0][jind-jzm][m+1];
          AIXY[0][jind][m] = pbz*AIXY[0][jind-jzm][m] -
                             pcz*AIXY[0][jind-jzm][m+1];
          AIXZ[0][jind][m] = pbz*AIXZ[0][jind-jzm][m] -
                             pcz*AIXZ[0][jind-jzm][m+1] +
                                 AIX[0][jind-jzm][m+1];
          AIYZ[0][jind][m] = pbz*AIYZ[0][jind-jzm][m] -
                             pcz*AIYZ[0][jind-jzm][m+1] +
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
          fprintf(outfile, "AIX[0][%d][%d] = %lf\tAIY[0][%d][%d] = %lf\tAIZ[0][%d][%d] = %lf\n", jind, m, AIX[0][jind][m], jind, m, AIY[0][jind][m], jind, m, AIZ[0][jind][m]);
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
          AI0[0][jind][m] = pby*AI0[0][jind-jym][m] -
                            pcy*AI0[0][jind-jym][m+1];
        for(m=0;m<=mmax-b-1;m++) {
          AIX[0][jind][m] = pby*AIX[0][jind-jym][m] -
                            pcy*AIX[0][jind-jym][m+1];
          AIY[0][jind][m] = pby*AIY[0][jind-jym][m] -
                            pcy*AIY[0][jind-jym][m+1] +
                                AI0[0][jind-jym][m+1];
          AIZ[0][jind][m] = pby*AIZ[0][jind-jym][m] -
                            pcy*AIZ[0][jind-jym][m+1];
          fprintf(outfile, "AIX[0][%d][%d] = %lf\tAIY[0][%d][%d] = %lf\tAIZ[0][%d][%d] = %lf\n", jind, m, AIX[0][jind][m], jind, m, AIY[0][jind][m], jind, m, AIZ[0][jind][m]);
        }
        for(m=0;m<=mmax-b-2;m++) {
          AIXX[0][jind][m] = pby*AIXX[0][jind-jym][m] -
                             pcy*AIXX[0][jind-jym][m+1];
          AIYY[0][jind][m] = pby*AIYY[0][jind-jym][m] -
                             pcy*AIYY[0][jind-jym][m+1] +
                               2*AIY[0][jind-jym][m+1];
          AIZZ[0][jind][m] = pby*AIZZ[0][jind-jym][m] -
                             pcy*AIZZ[0][jind-jym][m+1];
          AIXY[0][jind][m] = pby*AIXY[0][jind-jym][m] -
                             pcy*AIXY[0][jind-jym][m+1] +
                                 AIX[0][jind-jym][m+1];
          AIXZ[0][jind][m] = pby*AIXZ[0][jind-jym][m] -
                             pcy*AIXZ[0][jind-jym][m+1];
          AIYZ[0][jind][m] = pby*AIYZ[0][jind-jym][m] -
                             pcy*AIYZ[0][jind-jym][m+1] +
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
          fprintf(outfile, "AIX[0][%d][%d] = %lf\tAIY[0][%d][%d] = %lf\tAIZ[0][%d][%d] = %lf\n", jind, m, AIX[0][jind][m], jind, m, AIY[0][jind][m], jind, m, AIZ[0][jind][m]);
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
          AI0[0][jind][m] = pbx*AI0[0][jind-jxm][m] -
                            pcx*AI0[0][jind-jxm][m+1];
        for(m=0;m<=mmax-b-1;m++) {
          AIX[0][jind][m] = pbx*AIX[0][jind-jxm][m] -
                            pcx*AIX[0][jind-jxm][m+1] +
                                AI0[0][jind-jxm][m+1];
          AIY[0][jind][m] = pbx*AIY[0][jind-jxm][m] -
                            pcx*AIY[0][jind-jxm][m+1];
          AIZ[0][jind][m] = pbx*AIZ[0][jind-jxm][m] -
                            pcx*AIZ[0][jind-jxm][m+1];
          fprintf(outfile, "AIX[0][%d][%d] = %lf\tAIY[0][%d][%d] = %lf\tAIZ[0][%d][%d] = %lf\n", jind, m, AIX[0][jind][m], jind, m, AIY[0][jind][m], jind, m, AIZ[0][jind][m]);
        }
        for(m=0;m<=mmax-b-2;m++) {
          AIXX[0][jind][m] = pbx*AIXX[0][jind-jxm][m] -
                             pcx*AIXX[0][jind-jxm][m+1] +
                               2*AIX[0][jind-jxm][m+1];
          AIYY[0][jind][m] = pbx*AIYY[0][jind-jxm][m] -
                             pcx*AIYY[0][jind-jxm][m+1];
          AIZZ[0][jind][m] = pbx*AIZZ[0][jind-jxm][m] -
                             pcx*AIZZ[0][jind-jxm][m+1];
          AIXY[0][jind][m] = pbx*AIXY[0][jind-jxm][m] -
                             pcx*AIXY[0][jind-jxm][m+1] +
                                 AIY[0][jind-jxm][m+1];
          AIXZ[0][jind][m] = pbx*AIXZ[0][jind-jxm][m] -
                             pcx*AIXZ[0][jind-jxm][m+1] +
                                 AIZ[0][jind-jxm][m+1];
          AIYZ[0][jind][m] = pbx*AIYZ[0][jind-jxm][m] -
                             pcx*AIYZ[0][jind-jxm][m+1];
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
          fprintf(outfile, "AIX[0][%d][%d] = %lf\tAIY[0][%d][%d] = %lf\tAIZ[0][%d][%d] = %lf\n", jind, m, AIX[0][jind][m], jind, m, AIY[0][jind][m], jind, m, AIZ[0][jind][m]);
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
      else /* This should and will never happen */
        throw PsiException("There's some error in the AI_OSrecurs algorithm", __FILE__, __LINE__);
//        punt("There's some error in the AI_OSrecurs algorithm");
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
            AI0[iind][jind][m] = paz*AI0[iind-izm][jind][m] - 
                                 pcz*AI0[iind-izm][jind][m+1];
          for(m=0;m<=mmax-a-b-1;m++) {	/* Electric field integrals */
            AIX[iind][jind][m] = paz*AIX[iind-izm][jind][m] -
                                 pcz*AIX[iind-izm][jind][m+1];
            AIY[iind][jind][m] = paz*AIY[iind-izm][jind][m] -
                                 pcz*AIY[iind-izm][jind][m+1];
            AIZ[iind][jind][m] = paz*AIZ[iind-izm][jind][m] -
                                 pcz*AIZ[iind-izm][jind][m+1] +
                                     AI0[iind-izm][jind][m+1];
          fprintf(outfile, "AIX[%d][%d][%d] = %lf\tAIY[%d][%d][%d] = %lf\tAIZ[%d][%d][%d] = %lf\n", iind, jind, m, AIX[iind][jind][m], iind, jind, m, AIY[iind][jind][m], iind, jind, m, AIZ[iind][jind][m]);
          }
          for(m=0;m<=mmax-a-b-2;m++) {	/* Gradients of the electric field */
            AIXX[iind][jind][m] = paz*AIXX[iind-izm][jind][m] -
                                  pcz*AIXX[iind-izm][jind][m+1];
            AIYY[iind][jind][m] = paz*AIYY[iind-izm][jind][m] -
                                  pcz*AIYY[iind-izm][jind][m+1];
            AIZZ[iind][jind][m] = paz*AIZZ[iind-izm][jind][m] -
                                  pcz*AIZZ[iind-izm][jind][m+1] +
                                    2*AIZ[iind-izm][jind][m+1];
            AIXY[iind][jind][m] = paz*AIXY[iind-izm][jind][m] -
                                  pcz*AIXY[iind-izm][jind][m+1];
            AIXZ[iind][jind][m] = paz*AIXZ[iind-izm][jind][m] -
                                  pcz*AIXZ[iind-izm][jind][m+1] +
                                      AIX[iind-izm][jind][m+1];
            AIYZ[iind][jind][m] = paz*AIYZ[iind-izm][jind][m] -
                                  pcz*AIYZ[iind-izm][jind][m+1] +
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
          fprintf(outfile, "AIX[%d][%d][%d] = %lf\tAIY[%d][%d][%d] = %lf\tAIZ[%d][%d][%d] = %lf\n", iind, jind, m, AIX[iind][jind][m], iind, jind, m, AIY[iind][jind][m], iind, jind, m, AIZ[iind][jind][m]);
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
          fprintf(outfile, "AIX[%d][%d][%d] = %lf\tAIY[%d][%d][%d] = %lf\tAIZ[%d][%d][%d] = %lf\n", iind, jind, m, AIX[iind][jind][m], iind, jind, m, AIY[iind][jind][m], iind, jind, m, AIZ[iind][jind][m]);
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
            AI0[iind][jind][m] = pay*AI0[iind-iym][jind][m] -
                                 pcy*AI0[iind-iym][jind][m+1];
          for(m=0;m<=mmax-a-b-1;m++) {
            AIX[iind][jind][m] = pay*AIX[iind-iym][jind][m] -
                                 pcy*AIX[iind-iym][jind][m+1];
            AIY[iind][jind][m] = pay*AIY[iind-iym][jind][m] -
                                 pcy*AIY[iind-iym][jind][m+1] +
                                     AI0[iind-iym][jind][m+1];
            AIZ[iind][jind][m] = pay*AIZ[iind-iym][jind][m] -
                                 pcy*AIZ[iind-iym][jind][m+1];
          fprintf(outfile, "AIX[%d][%d][%d] = %lf\tAIY[%d][%d][%d] = %lf\tAIZ[%d][%d][%d] = %lf\n", iind, jind, m, AIX[iind][jind][m], iind, jind, m, AIY[iind][jind][m], iind, jind, m, AIZ[iind][jind][m]);
          }
          for(m=0;m<=mmax-a-b-2;m++) {
            AIXX[iind][jind][m] = pay*AIXX[iind-iym][jind][m] -
                                  pcy*AIXX[iind-iym][jind][m+1];
            AIYY[iind][jind][m] = pay*AIYY[iind-iym][jind][m] -
                                  pcy*AIYY[iind-iym][jind][m+1] +
                                    2*AIY[iind-iym][jind][m+1];
            AIZZ[iind][jind][m] = pay*AIZZ[iind-iym][jind][m] -
                                  pcy*AIZZ[iind-iym][jind][m+1];
            AIXY[iind][jind][m] = pay*AIXY[iind-iym][jind][m] -
                                  pcy*AIXY[iind-iym][jind][m+1] +
                                      AIX[iind-iym][jind][m+1];
            AIXZ[iind][jind][m] = pay*AIXZ[iind-iym][jind][m] -
                                  pcy*AIXZ[iind-iym][jind][m+1];
            AIYZ[iind][jind][m] = pay*AIYZ[iind-iym][jind][m] -
                                  pcy*AIYZ[iind-iym][jind][m+1] +
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
          fprintf(outfile, "AIX[%d][%d][%d] = %lf\tAIY[%d][%d][%d] = %lf\tAIZ[%d][%d][%d] = %lf\n", iind, jind, m, AIX[iind][jind][m], iind, jind, m, AIY[iind][jind][m], iind, jind, m, AIZ[iind][jind][m]);
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
          fprintf(outfile, "AIX[%d][%d][%d] = %lf\tAIY[%d][%d][%d] = %lf\tAIZ[%d][%d][%d] = %lf\n", iind, jind, m, AIX[iind][jind][m], iind, jind, m, AIY[iind][jind][m], iind, jind, m, AIZ[iind][jind][m]);
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
            AI0[iind][jind][m] = pax*AI0[iind-ixm][jind][m] -
                                 pcx*AI0[iind-ixm][jind][m+1];
          for(m=0;m<=mmax-a-b-1;m++) {	/* Electric field integrals */
            AIX[iind][jind][m] = pax*AIX[iind-ixm][jind][m] -
                                 pcx*AIX[iind-ixm][jind][m+1] +
                                     AI0[iind-ixm][jind][m+1];
            AIY[iind][jind][m] = pax*AIY[iind-ixm][jind][m] -
                                 pcx*AIY[iind-ixm][jind][m+1];
            AIZ[iind][jind][m] = pax*AIZ[iind-ixm][jind][m] -
                                 pcx*AIZ[iind-ixm][jind][m+1];
          fprintf(outfile, "AIX[%d][%d][%d] = %lf\tAIY[%d][%d][%d] = %lf\tAIZ[%d][%d][%d] = %lf\n", iind, jind, m, AIX[iind][jind][m], iind, jind, m, AIY[iind][jind][m], iind, jind, m, AIZ[iind][jind][m]);
          }
          for(m=0;m<=mmax-a-b-2;m++) {	/* Gradients of the electric field */
            AIXX[iind][jind][m] = pax*AIXX[iind-ixm][jind][m] -
                                  pcx*AIXX[iind-ixm][jind][m+1] +
                                    2*AIX[iind-ixm][jind][m+1];
            AIYY[iind][jind][m] = pax*AIYY[iind-ixm][jind][m] -
                                  pcx*AIYY[iind-ixm][jind][m+1];
            AIZZ[iind][jind][m] = pax*AIZZ[iind-ixm][jind][m] -
                                  pcx*AIZZ[iind-ixm][jind][m+1];
            AIXY[iind][jind][m] = pax*AIXY[iind-ixm][jind][m] -
                                  pcx*AIXY[iind-ixm][jind][m+1] +
                                      AIY[iind-ixm][jind][m+1];
            AIXZ[iind][jind][m] = pax*AIXZ[iind-ixm][jind][m] -
                                  pcx*AIXZ[iind-ixm][jind][m+1] +
                                      AIZ[iind-ixm][jind][m+1];
            AIYZ[iind][jind][m] = pax*AIYZ[iind-ixm][jind][m] -
                                  pcx*AIYZ[iind-ixm][jind][m+1];
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
          fprintf(outfile, "AIX[%d][%d][%d] = %lf\tAIY[%d][%d][%d] = %lf\tAIZ[%d][%d][%d] = %lf\n", iind, jind, m, AIX[iind][jind][m], iind, jind, m, AIY[iind][jind][m], iind, jind, m, AIZ[iind][jind][m]);
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
          fprintf(outfile, "AIX[%d][%d][%d] = %lf\tAIY[%d][%d][%d] = %lf\tAIZ[%d][%d][%d] = %lf\n", iind, jind, m, AIX[iind][jind][m], iind, jind, m, AIY[iind][jind][m], iind, jind, m, AIZ[iind][jind][m]);
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
        else /* This should and will never happen */
          throw PsiException("There's some error in the AI_OSrecurs algorithm", __FILE__, __LINE__);
//          punt("There's some error in the AI_OSrecurs algorithm");
      }
    }
  free(F);
}

}} // namespace psi::oeprop

namespace {
using namespace psi::oeprop;

	/* This function computes infamous integral Fn(t). For its definition 
	   see Obara and Saika paper, or Shavitt's chapter in the 
	   "Methods in Computational Physics" book (see reference below).
	   This piece of code is from Dr. Justin Fermann's program CINTS */

void calc_f(double *F, int n, double t)
{
  int i, m, k;
  int m2;
  double t2;
  double num;
  double sum;
  double term1, term2;
  static double K = 1.0/M_2_SQRTPI;
  double et;


  if (t>20.0){   /* For big t's do upward recursion */
    t2 = 2*t;
    et = exp(-t);
    t = sqrt(t);
    F[0] = K*erf(t)/t;
    for(m=0; m<=n-1; m++){
      F[m+1] = ((2*m + 1)*F[m] - et)/(t2);
    }
  }
  else {	/* For smaller t's compute F with highest n using 
  		   asymptotic series (see I. Shavitt in 
  		   "Methods in Computational Physics", ed. B. Alder etal, 
  		   vol 2, 1963, page 8 */
    et = exp(-t);
    t2 = 2*t;
    m2 = 2*n;
    num = df[m2];
    i=0;
    sum = 1.0/(m2+1);
    do{
      i++;
      num = num*t2;
      term1 = num/df[m2+2*i+2];
      sum += term1;
    } while (fabs(term1) > EPS && i < MAXFACT);
    F[n] = sum*et;
    for(m=n-1;m>=0;m--){	/* And then do downward recursion */
      F[m] = (t2*F[m+1] + et)/(2*m+1);
    }
  }
}

} // namespace
