/*! \file 
    \ingroup (CSCF)
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.2  2002/12/06 15:50:32  crawdad
 * Changed all exit values to PSI_RETURN_SUCCESS or PSI_RETURN_FAILURE as
 * necessary.  This is new for the PSI3 execution driver.
 * -TDC
 *
/* Revision 1.1.1.1  2000/02/04 22:52:33  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.1  1999/11/02 23:55:57  localpsi
/* Shawn Brown - (11/2/99) Modified to the code in a few major ways.
/*
/* 1.  Added the capability to do UHF.  All of the features available with the
/* other refrences have been added for UHF.
/*
/* 2.  For UHF, I had to alter the structure of file30. (See cleanup.c for a
/* map)  This entailed adding a pointer array right after the header in the SCF
/* section of file30 that pointed to all of the data for the SCF caclulation.
/* Functions were added to libfile30 to account for this and they are
/* incorporated in this code.
/*
/* 3.  Updated and fixed all of the problems associated with my previous
/* guessing code.  The code no longer uses OPENTYPE to specify the type of
/* occupation.  The keword REFERENCE and MULTP can now be used to indicate any
/* type of calculation.  (e.g. ROHF with MULTP of 1 is an open shell singlet
/* ROHF calculation)  This code was moved to occ_fun.c.  The code can also
/* guess at any multplicity in a highspin case, provided enough electrons.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:26  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
/*
 * Revision 1.3  1996/11/26  05:20:52  sherrill
 * Add casts in front of init_int_array() calls to avoid compiler warnings
 * (really should include libciomr everywhere instead but that causes more
 * warnings because wwritw's need casts to char *).  Also fixed problem
 * where phase() was trying to free unallocated memory.
 *
 * Revision 1.2  1995/11/07  19:24:05  sherrill
 * Replace all initializations of integer arrays with calls to init_int_array().
 * There had been problems with using init_array() because the division by
 * the doubles-to-int factor was rounding off (integer division).
 *
 * Revision 1.1  1991/06/15  20:22:25  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: findit_uhf.cc 3592 2007-09-28 13:01:33Z evaleev $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

extern double *pa, *pb;
extern int *inext;
extern unsigned int *lbij,*lbkl;
extern int intmx;
static int old_nint=25920;

void findit_uhf(int ii, int jj, int kk, int ll, int ism, int ksm, double value, int iab)
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
   if (!inext) inext = (int *) init_int_array(keep+intmx);

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
   
   /*value = (lij == lkl) ? value*0.5 : value;*/
   
   switch(iab) {
   case 1:
       pa[i] += value;
       break;
   case 2:
       pb[i] += 0.50*value;
       break;
   case 3:
       pa[i] += value;
       pb[i] += 0.50*value;
       break;
   case 4:
       pa[i] +=value;
       pb[i] +=value;
       break;
   case 5:
       pb[i] +=value;
       break;
   }
   old_nint=cscf_nint;
}

}} // namespace psi::cscf
