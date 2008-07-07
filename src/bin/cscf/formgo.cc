/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.2  2001/06/29 20:39:29  evaleev
 * Modified cscf to use libpsio to store supermatrix files.
 *
/* Revision 1.1.1.1  2000/02/04 22:52:30  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.4  1999/11/17 19:40:46  evaleev
/* Made all the adjustments necessary to have direct UHF working. Still doesn't work though..
/*
/* Revision 1.3  1999/11/02 23:55:57  localpsi
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
/* Revision 1.2  1999/08/17 19:04:15  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:26  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id: formgo.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

extern double *gtmp,*ptmp; // defined in formgc.cc
extern double *gtmpo2,*ptmpo2; // defined in formg2.cc
double *gtmpo,*ptmpo;
extern int num_bufs_c,num_bufs_o,readflgc,readflgo;
extern struct c_pkints {
         int ij;
         int kl;
         double pval;
         } *c_outbuf;
extern struct o_pkints {
         int ij;
         int kl;
         double pval;
         double qval;
         } *o_outbuf;

extern int lasto,lastc;
int wherec=0;
int whereo=0;
int *int_nums_c;
int *int_nums_o;

void formg_open()
{
   register int i,j,k,joff,nn;
   register int ij,kl;
   int ilast,num;
   int tmpsiz;
   double dotest,tmpval;
   struct o_pkints *o_temp;
   struct c_pkints *c_temp;

   tmpsiz=ioff[nbasis];

   if(gtmp == NULL) {
       gtmp = (double *) init_array(tmpsiz);
       gtmpo = (double *) init_array(tmpsiz);
       ptmp = (double *) init_array(tmpsiz);
       ptmpo = (double *) init_array(tmpsiz);
       if(uhf){
	   gtmpo2 = (double *) init_array(tmpsiz);
	   ptmpo2 = (double *) init_array(tmpsiz);
       }
   }
   else {
       bzero(gtmp,sizeof(double)*tmpsiz);
       bzero(gtmpo,sizeof(double)*tmpsiz);
       if(uhf)
	   bzero(gtmpo2,sizeof(double)*tmpsiz); 
   }
   
   for(k=joff=0; k < num_ir ; k++) {
       if(nn=scf_info[k].num_so) {
	   for(i=0; i < nn ; i++) {
	       if(twocon) {
		   for(j=0; j <= i ; j++) {
		       ptmp[ioff[i+joff]+j+joff] = 
			   scf_info[k].pmat2[ioff[i]+j];
		       ptmpo[ioff[i+joff]+j+joff] = 
			   scf_info[k].pmato2[ioff[i]+j];
		   }
               }
	       else if(uhf){
		   for(j=0; j <= i ; j++) {
		       ptmp[ioff[i+joff]+j+joff] = 
			   scf_info[k].dpmat[ioff[i]+j];
		       ptmpo[ioff[i+joff]+j+joff] = 
			   -spin_info[0].scf_spin[k].dpmat[ioff[i]+j];
		       ptmpo2[ioff[i+joff]+j+joff] = 
			   -spin_info[1].scf_spin[k].dpmat[ioff[i]+j];
		   } 
	       }
	       else {
		   for(j=0; j <= i ; j++) {
		       ptmp[ioff[i+joff]+j+joff] = 
			   scf_info[k].dpmat[ioff[i]+j];
		       ptmpo[ioff[i+joff]+j+joff] = 
			   scf_info[k].dpmato[ioff[i]+j];
		   }
               }
	   }
	   joff += nn;
       }
   }

   if(!wherec) {
      /* int_nums_o = (int *) init_array(num_bufs_o+1); */
      int_nums_o = (int *) init_int_array(num_bufs_o+1);
      /* int_nums_c = (int *) init_array(num_bufs_c+1); */
      int_nums_c = (int *) init_int_array(num_bufs_c+1);
      for(i=1; i < num_bufs_o ; i++) int_nums_o[i]=maxbuf;
      for(i=1; i < num_bufs_c ; i++) int_nums_c[i]=maxbuf;
      int_nums_o[num_bufs_o]=lasto;
      int_nums_c[num_bufs_c]=lastc;
      whereo=num_bufs_o;
      wherec=num_bufs_c;
      }

   num=int_nums_o[whereo];
   for (j=0; j < num_bufs_o ; j++) {
      o_temp = o_outbuf;

      for (i=num; i ; i--,o_temp++) {
         ij = (*o_temp).ij;
         kl = (*o_temp).kl;
         tmpval = (*o_temp).pval;
         dotest = (*o_temp).qval;
	 
	 if(uhf){
	     gtmp[ij] += ptmp[kl]*tmpval;
	     gtmpo[ij] += ptmpo[kl]*dotest;
	     gtmpo2[ij] += ptmpo2[kl]*dotest;
	     if(ij!=kl){
		 gtmp[kl] += ptmp[ij]*tmpval;
		 gtmpo[kl] += ptmpo[ij]*dotest;
		 gtmpo2[kl] += ptmpo2[ij]*dotest;
	     }
	 }
	 else{
	     gtmp[ij] += ptmp[kl]*tmpval;
	     gtmp[kl] += ptmp[ij]*tmpval;
	     gtmpo[ij] += ptmpo[kl]*dotest;
	     gtmpo[kl] += ptmpo[ij]*dotest;
         }
      }
   
      if (readflgo && j < num_bufs_o-1) {
         if(whereo==num_bufs_o) {
	    PKmat.bufpos = PSIO_ZERO;
	    whereo=0;
            }
         whereo++;
         num=int_nums_o[whereo];
	 psio_read(PKmat.unit, PKmat.key, (char *) o_outbuf, sizeof(struct o_pkints)*num,
		   PKmat.bufpos, &(PKmat.bufpos));
         }
      }

   num=int_nums_c[wherec];
   for (j=0; j < num_bufs_c ; j++) {
      c_temp = c_outbuf;

      for (i=num; i ; i--,c_temp++) {
         ij = (*c_temp).ij;
         kl = (*c_temp).kl;
         tmpval = (*c_temp).pval;

         gtmp[ij] += ptmp[kl]*tmpval;
         gtmp[kl] += ptmp[ij]*tmpval;
         }
   
      if (readflgc && j < num_bufs_c-1) {
         if(wherec==num_bufs_c) {
	    Pmat.bufpos = PSIO_ZERO;
            wherec=0;
            }
         wherec++;
         num=int_nums_c[wherec];
	 psio_read(Pmat.unit, Pmat.key, (char *) c_outbuf, sizeof(struct c_pkints)*num,
		   Pmat.bufpos, &(Pmat.bufpos));
         }
      }

   for(k=joff=0; k < num_ir ; k++) {
      if(nn=scf_info[k].num_so) {
         for(i=0; i < nn ; i++) {
            if(twocon) {
               for(j=0; j <= i ; j++) {
                  scf_info[k].dpmat[ioff[i]+j] = 
		      gtmp[ioff[i+joff]+j+joff];
                  scf_info[k].dpmato[ioff[i]+j] = 
		      gtmpo[ioff[i+joff]+j+joff];
                  }
               }
	    else if(uhf){
		   for(j=0; j <=i; j++){
		       spin_info[0].scf_spin[k].gmat[ioff[i]+j] +=
			   gtmp[ioff[i+joff]+j+joff]
			   +gtmpo[ioff[i+joff]+j+joff];
		       spin_info[1].scf_spin[k].gmat[ioff[i]+j] +=
			   gtmp[ioff[i+joff]+j+joff]
			   +gtmpo2[ioff[i+joff]+j+joff];
		   }
	    }
            else {
		for(j=0; j <= i ; j++) {
		    scf_info[k].gmat[ioff[i]+j] += 
			gtmp[ioff[i+joff]+j+joff];
		    scf_info[k].gmato[ioff[i]+j] += 
			gtmpo[ioff[i+joff]+j+joff];
		}
	    }
	 }
         joff += nn;
      }
   }
}

void formg_open_free() {
  free(gtmp); gtmp = NULL;
  free(gtmpo); gtmpo = NULL;
  free(gtmpo2); gtmpo2 = NULL;
  free(ptmp); ptmp = NULL;
  free(ptmpo); ptmpo = NULL;
  free(ptmpo2); ptmpo2 = NULL;
  wherec=0;
  whereo=0;
  free(int_nums_c); int_nums_c = NULL;
  free(int_nums_o); int_nums_o = NULL;
}
}} // namespace psi::cscf
