/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.2  2002/12/06 15:50:32  crawdad
 * Changed all exit values to PSI_RETURN_SUCCESS or PSI_RETURN_FAILURE as
 * necessary.  This is new for the PSI3 execution driver.
 * -TDC
 *
/* Revision 1.1.1.1  2000/02/04 22:52:30  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.3  1999/08/17 19:04:14  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.2  1999/07/24 18:13:51  crawdad
/* Renamed variable "nint" to "cscf_nint" to avoid DEC compiler type conflict.
/* -Daniel
/*
 * Revision 1.1.1.1  1999/04/12  16:59:26  evaleev
 * Added a version of CSCF that can work with CINTS.
 * -Ed
 * */

static char *rcsid = "$Id: findit.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void findit(int ii, int jj, int kk, int ll, int ism, int ksm, double value, int iab)
{
   register int i,j;
   unsigned int *ijtmp, *kltmp;
   int *nxtmp;
   int p1,p2,p3;
   int noi,nok;
   int next,start;
   int lij,lkl;
   int keep=128;
   int keep2=127;
   int d2i = sizeof(double)/sizeof(int);
   double *patmp, *pbtmp;

   if(nbasis > 150) {
     keep=1024;
     keep2=1023;
     }

   /* if(!inext) inext = (int *) init_array((int) (keep+intmx)/d2i); */
   if (inext == NULL) inext = (int *) init_int_array(keep+intmx);

   noi = scf_info[ism].nopen+scf_info[ism].nhalf;
   nok = scf_info[ksm].nopen+scf_info[ksm].nhalf;

   lij = ioff[ii]+jj;
   lkl = ioff[kk]+ll;

   if(!cscf_nint) {
      bzero(inext,sizeof(int)*old_nint);
      bzero(&inext[intmx],sizeof(int)*keep);
      }

   start = 2*lij + lkl;
   start = (start & keep2) + intmx;

L1:
   next=inext[start];
   if(next) {
      if (lbij[next-1] == lij && lbkl[next-1] == lkl) i=next-1;
      else {
         start = next;
         goto L1;
         }
      }
   else {
      i=cscf_nint;
      if(cscf_nint >= intmx) {
        fprintf(outfile,"\n  increasing size of buffers in findit\n");
        fprintf(outfile,"  intmx was %d, is %d\n",intmx,intmx*2);
        fflush(outfile);
        intmx*=2;

 /* i don't use realloc because strange things were happening */

        /* nxtmp = (int *) init_array((int) (keep+intmx)/d2i); */
        nxtmp = (int *) init_int_array(keep+intmx);
        bcopy(inext,nxtmp,(int)sizeof(int)*(intmx/2));
        for(j=0; j < keep ; j++) nxtmp[j+intmx]=inext[j+intmx/2];
        free(inext);
        inext=nxtmp;

        /* ijtmp = (unsigned int *) init_array(intmx/d2i); */
        ijtmp = (unsigned int *) init_int_array(intmx);
        bcopy(lbij,ijtmp,sizeof(int)*(intmx/2));
        free(lbij);
        lbij = ijtmp;

        /* kltmp = (unsigned int *) init_array(intmx/d2i); */
        kltmp = (unsigned int *) init_int_array(intmx);
        bcopy(lbkl,kltmp,sizeof(int)*(intmx/2));
        free(lbkl);
        lbkl = kltmp;

        patmp = (double *) init_array(intmx);
        bcopy(pa,patmp,sizeof(double)*(intmx/2));
        free(pa);
        pa = patmp;

        pbtmp = (double *) init_array(intmx);
        bcopy(pb,pbtmp,sizeof(double)*(intmx/2));
        free(pb);
        pb = pbtmp;

        if(inext==NULL || lbij==NULL || lbkl==NULL) {
          fprintf(outfile,"\n pathological problems with realloc in findit\n");
          fprintf(outfile," try upping intmx to %d\n",intmx);
          fflush(outfile);
          exit(PSI_RETURN_FAILURE);
          }

        start = 2*lij + lkl;
        start = (start & keep2) + intmx;
        }
      inext[start] = ++cscf_nint;
      lbij[i] = lij;
      lbkl[i] = lkl;
      pa[i] = pb[i] = 0.0;
      }

   value = (lij == lkl) ? value*0.5 : value;

   /* EFV 10/27/98
   fprintf(outfile,"%3d %3d  %3d  %3d  %3d  %3d  %lf\n",lbij[i],lbkl[i],ii,jj,kk,ll,value); */

   switch(iab) {
      case 1:
         pa[i] += value;
         if(noi && nok) {
            if(special) {
               p1 = MAX0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p2 = MIN0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p3 = ioff[p1]+p2;
               pb[i] += alpha[p3]*value;
               }
            }
         break;
      case 2:
         pa[i] -= 0.25*value;
         if(noi && nok) {
            if(hsos) pb[i] += 0.25*value;
            else if(singlet) {
               if (ism != ksm) pb[i] -= 0.75*value;
               else pb[i] += 0.25*value;
               }
            else if(twocon) {
               if (ism != ksm) pb[i] += 0.25*value;
               }
            else {
               p1 = MAX0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p2 = MIN0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p3 = ioff[p1]+p2;
               pb[i] += 0.25*beta[p3]*value;
               }
            }
         break;
      case 3:
         pa[i] += 0.75*value;
         if(noi && nok) {
            if(hsos) pb[i] += 0.25*value;
            else if(singlet) {
               if (ism != ksm) pb[i] -= 0.75*value;
               else pb[i] += 0.25*value;
               }
            else if(twocon) {
               if (ism != ksm) pb[i] += 0.25*value;
               }
            else {
               p1 = MAX0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p2 = MIN0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p3 = ioff[p1]+p2;
               pb[i] += (alpha[p3] + 0.25*beta[p3])*value;
               }
            }
         break;
      case 4:
         pa[i] -= 0.5*value;
         if(noi && nok) {
            if(hsos) pb[i] = 0.5*value;
            else if(singlet) {
               if (ism != ksm) pb[i] -= 1.5*value;
               else pb[i] += 0.5*value;
               }
            else if(twocon) {
               if (ism != ksm) pb[i] += 0.5*value;
               }
            else {
               p1 = MAX0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p2 = MIN0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p3 = ioff[p1]+p2;
               pb[i] += 0.5*beta[p3]*value;
               }
            }
         break;
      case 5:
         pa[i] = 0.5*value;
         if(noi && nok) {
            if(hsos) pb[i] = 0.5*value;
            else if(singlet) {
               if (ism != ksm) pb[i] = -1.5*value;
               else pb[i] = 0.5*value;
               }
            else if(twocon) {
               if (ism != ksm) pb[i] += 0.5*value;
               }
            else {
               p1 = MAX0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p2 = MIN0(scf_info[ism].os_num,scf_info[ksm].os_num);
               p3 = ioff[p1]+p2;
               pb[i] = (alpha[p3] + 0.5*beta[p3])*value;
               }
            }
      }
   old_nint=cscf_nint;
   }

}} // namespace psi::cscf
