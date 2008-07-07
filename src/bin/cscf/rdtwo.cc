/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.5  2004/01/02 06:37:12  crawdad
 * Merging psi-3-2 branch (from release tag psi-3-2-rc-2 to head at psi-3-2-0)
 * into main trunk.  This code compiles and runs correctly on sirius.
 * -TDC
 *
/* Revision 1.4.4.1  2003/12/31 01:59:54  crawdad
/* Removed use_iwl option, which was deprecated anyway.  cscf now depends only
/* on libpsio for its I/O functions.
/* -TDC
/*
/* Revision 1.4  2002/12/06 15:50:32  crawdad
/* Changed all exit values to PSI_RETURN_SUCCESS or PSI_RETURN_FAILURE as
/* necessary.  This is new for the PSI3 execution driver.
/* -TDC
/*
/* Revision 1.3  2002/04/03 02:06:01  janssen
/* Finish changes to use new include paths for libraries.
/*
/* Revision 1.2  2002/03/25 02:17:36  janssen
/* Get rid of tmpl.  Use new naming scheme for libipv1 includes.
/*
/* Revision 1.1.1.1  2000/02/04 22:52:29  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.6  1999/11/02 23:55:59  localpsi
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
/* Revision 1.5  1999/08/17 19:04:16  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.4  1999/07/24 18:13:54  crawdad
/* Renamed variable "nint" to "cscf_nint" to avoid DEC compiler type conflict.
/* -Daniel
/*
 * Revision 1.3  1999/07/14  18:50:40  evaleev
 * Very important fix of a careless mistake I introduced before. (intmx has to
 * be the same as old_nint).
 *
 * Revision 1.2  1999/06/13  20:30:04  evaleev
 * Fixed problems with Makefile.
 *
/* Revision 1.1.1.1  1999/04/12 16:59:27  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id: rdtwo.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>

namespace psi { namespace cscf {

double *pa, *pb;
int *inext;
unsigned int *lbij,*lbkl;
int intmx=25920;
/*int intmx=51840;*/

extern void findit(int ii, int jj, int kk, int ll, int ism, int ksm, double value, int iab);
extern void findit_uhf(int ii, int jj, int kk, int ll, int ism, int ksm, double value, int iab);
extern void packit_closed(unsigned int* lbij, unsigned int* lbkl, int endflg);
extern void packit_open(unsigned int* lbij, unsigned int* lbkl, int endflg);

void rdtwo()
{
   int ilsti, nbuf;
   int ibufsz = 8942;
   int ibufqt = 2236;
   int i;
   int iii;
   int fi;
   int errcod;
   int ior, ism, jor, jsm;
   int kor, ksm, lor, lsm;
   int ii,jj,kk,ll;
   int p1,p2;
   int d2i = sizeof(double)/sizeof(int);
   double pki_int;
   int *tmp;
   union psi_buffer inbuf;
   struct iwlbuf ERIIN;
   int pk_flush;
   double *ints;

   cscf_nint=0;

   if(nbasis > 150) intmx *= 4;

   iwl_buf_init(&ERIIN,itap33,0.0,1,1);
   fprintf(outfile,"  reading integrals in the IWL format from files 33,35,36,37\n");
   pa = (double *) init_array(intmx);
   pb = (double *) init_array(intmx);

   /* lbij = (unsigned int *) init_array((int) intmx/d2i); */
   /* lbkl = (unsigned int *) init_array((int) intmx/d2i); */
   lbij = (unsigned int *) init_int_array(intmx);
   lbkl = (unsigned int *) init_int_array(intmx);

   do {
	ilsti = ERIIN.lastbuf;
	nbuf = ERIIN.inbuf;

      if (print & 8) fprintf(outfile,"%5d\n",nbuf);

      fi = 0;
      for (i=0 ; i < nbuf ; i++,tmp += 2) {
	ii = ERIIN.labels[fi];
	jj = ERIIN.labels[fi+1];
	kk = ERIIN.labels[fi+2];
	ll = ERIIN.labels[fi+3];
	pki_int = ERIIN.values[i];
	pk_flush = 0;
	/* If the first index is negative - it's the last integral necessary for a pk-block
	   Set the flushing flag.. Hey, and we need the positive index back */
	if (ii < 0) {
	  pk_flush = 1;
	  ii = abs(ii);
	}
	ism = so2symblk[ii];
	jsm = so2symblk[jj];
	ksm = so2symblk[kk];
	lsm = so2symblk[ll];
	ior = ii - scf_info[ism].ideg;
	jor = jj - scf_info[jsm].ideg;
	kor = kk - scf_info[ksm].ideg;
	lor = ll - scf_info[lsm].ideg;
	fi += 4;

	if (print & 128)
	  fprintf(outfile," i = %d j = %d k = %d l = %d  %lf\n",ii,jj,kk,ll,pki_int);
	if(!uhf){
	    if (ism == jsm && ksm == lsm && ism == ksm) {
		if (ior == jor && ior == kor || jor == kor && jor == lor) {
		    findit(ii,jj,kk,ll,ism,ksm,pki_int,5);
		}
		else if (ior == kor || jor == lor) {
		    findit(ii,jj,kk,ll,ism,ksm,pki_int,3);
		    
		    p1 = MAX0(jj,ll);
		    p2 = MIN0(jj,ll);
		    
		    findit(ii,kk,p1,p2,ism,ksm,pki_int,4);
		}
		else if (jor == kor) {
		    findit(ii,jj,kk,ll,ism,ksm,pki_int,3);
		    
		    findit(ii,ll,jj,kk,ism,ksm,pki_int,4);
		}
		else if (ior == jor || kor == lor) {
		    findit(ii,jj,kk,ll,ism,ksm,pki_int,1);
		    
		    p1 = MAX0(jj,ll);
		    p2 = MIN0(jj,ll);
		    
		    findit(ii,kk,p1,p2,ism,ksm,pki_int,2);
		}
		else {
		    findit(ii,jj,kk,ll,ism,ksm,pki_int,1);
		    
		    p1 = MAX0(jj,ll);
		    p2 = MIN0(jj,ll);
		    
		    findit(ii,kk,p1,p2,ism,ksm,pki_int,2);
		    
		    p1 = MAX0(jj,kk);
		    p2 = MIN0(jj,kk);
		    
		    findit(ii,ll,p1,p2,ism,ksm,pki_int,2);
		}
	    }
	    else if (ism == jsm) {
		findit(ii,jj,kk,ll,ism,ksm,pki_int,1);
	    }
	    else if (ism == ksm) {
		if (ior == kor || jor == lor) {
		    p1 = MAX0(jj,ll);
		    p2 = MIN0(jj,ll);
		    
		    findit(ii,kk,p1,p2,ism,jsm,pki_int,4);
		}
		else {
		    p1 = MAX0(jj,ll);
		    p2 = MIN0(jj,ll);
		    
		    findit(ii,kk,p1,p2,ism,jsm,pki_int,2);
		}
	    }
	}
	else{
	    if (ism == jsm && ksm == lsm && ism == ksm) {
		if (ior == jor && ior == kor || jor == kor && jor == lor) {
		    findit_uhf(ii,jj,kk,ll,ism,ksm,pki_int,4);
		}
		else if (ior == kor || jor == lor) {
		    findit_uhf(ii,jj,kk,ll,ism,ksm,pki_int,3);
		    
		    p1 = MAX0(jj,ll);
		    p2 = MIN0(jj,ll);
		    
		    findit_uhf(ii,kk,p1,p2,ism,ksm,pki_int,5);
		}
		else if (jor == kor) {
		    findit_uhf(ii,jj,kk,ll,ism,ksm,pki_int,3);
		    
		    findit_uhf(ii,ll,jj,kk,ism,ksm,pki_int,5);
		}
		else if (ior == jor || kor == lor) {
		    findit_uhf(ii,jj,kk,ll,ism,ksm,pki_int,1);
		    
		    p1 = MAX0(jj,ll);
		    p2 = MIN0(jj,ll);
		    
		    findit_uhf(ii,kk,p1,p2,ism,ksm,pki_int,2);
		}
		else {
		    findit_uhf(ii,jj,kk,ll,ism,ksm,pki_int,1);
		    
		    p1 = MAX0(jj,ll);
		    p2 = MIN0(jj,ll);
		    
		    findit_uhf(ii,kk,p1,p2,ism,ksm,pki_int,2);
		    
		    p1 = MAX0(jj,kk);
		    p2 = MIN0(jj,kk);
		    
		    findit_uhf(ii,ll,p1,p2,ism,ksm,pki_int,2);
		}
	    }
	    else if (ism == jsm) {
		findit_uhf(ii,jj,kk,ll,ism,ksm,pki_int,1);
	    }
	    else if (ism == ksm) {
		if (ior == kor || jor == lor) {
		  
		    p1 = MAX0(jj,ll);
		    p2 = MIN0(jj,ll);

		    findit_uhf(ii,kk,p1,p2,ism,ksm,pki_int,5);
		}
		else {
		    p1 = MAX0(jj,ll);
		    p2 = MIN0(jj,ll);
		    
		    findit_uhf(ii,kk,p1,p2,ism,jsm,pki_int,2);
		}
	    }
	}
	if(pk_flush && cscf_nint) {
            if(iopen || uhf) packit_open(lbij,lbkl,0);
            else packit_closed(lbij,lbkl,0);
	}
      }
      if (!ilsti)
	  iwl_buf_fetch(&ERIIN);
   } while(!ilsti);
   
   if(iopen || uhf) packit_open(lbij,lbkl,1);
   else packit_closed(lbij,lbkl,1);
   
   if(!iopen || !uhf) fprintf(outfile,"  wrote %d integrals to file92\n",num_ints);
   free(pa);
   free(pb);
   free(lbij);
   free(lbkl);
   free(inext);
   pa = pb = NULL;
   lbij = lbkl = NULL;
   inext = NULL;
   iwl_buf_close(&ERIIN,!delete_2e);
}

}} // namespace psi::cscf
