/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.5  2007/04/05 15:45:25  crawdad
 * Fixed a few memory leaks identified by valgrind. -TDC
 *
 * Revision 1.4  2002/12/06 15:50:32  crawdad
 * Changed all exit values to PSI_RETURN_SUCCESS or PSI_RETURN_FAILURE as
 * necessary.  This is new for the PSI3 execution driver.
 * -TDC
 *
 * Revision 1.3  2001/06/29 20:39:29  evaleev
 * Modified cscf to use libpsio to store supermatrix files.
 *
 * Revision 1.2  2000/06/22 22:15:01  evaleev
 * Modifications for KS DFT. Reading in XC Fock matrices and XC energy in formg_direct need to be uncommented (at present those are not produced by CINTS yet).
 *
 * Revision 1.1.1.1  2000/02/04 22:52:31  evaleev
 * Started PSI 3 repository
 *
 * Revision 1.4  1999/11/02 23:55:59  localpsi
 * Shawn Brown - (11/2/99) Modified to the code in a few major ways.
 *
 * 1.  Added the capability to do UHF.  All of the features available with the
 * other refrences have been added for UHF.
 *
 * 2.  For UHF, I had to alter the structure of file30. (See cleanup.c for a
 * map)  This entailed adding a pointer array right after the header in the SCF
 * section of file30 that pointed to all of the data for the SCF caclulation.
 * Functions were added to libfile30 to account for this and they are
 * incorporated in this code.
 *
 * 3.  Updated and fixed all of the problems associated with my previous
 * guessing code.  The code no longer uses OPENTYPE to specify the type of
 * occupation.  The keword REFERENCE and MULTP can now be used to indicate any
 * type of calculation.  (e.g. ROHF with MULTP of 1 is an open shell singlet
 * ROHF calculation)  This code was moved to occ_fun.c.  The code can also
 * guess at any multplicity in a highspin case, provided enough electrons.
 *
 * Revision 1.3  1999/08/17 19:04:16  evaleev
 * Changed the default symmetric orthogonalization to the canonical
 * orthogonalization. Now, if near-linear dependencies in the basis are found,
 * eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
 * left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
 * longer a square matrix. Had to rework some routines in libfile30, and add some.
 * The progrem prints out a warning if near-linear dependencies are found. TRANSQT
 * and a whole bunch of other codes has to be fixed to work with such basis sets.
 *
 * Revision 1.2  1999/07/24 18:13:53  crawdad
 * Renamed variable "nint" to "cscf_nint" to avoid DEC compiler type conflict.
 * -Daniel
 *
 * Revision 1.1.1.1  1999/04/12  16:59:27  evaleev
 * Added a version of CSCF that can work with CINTS.
 * -Ed
 * */

static char *rcsid = "$Id: packit_o.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void packit_open(unsigned int* lbij, unsigned int* lbkl, int endflg)
{
   int i,j,k,joff,ij,kl,l,lmax,ijkl;
   int tmpsiz,nn;
   double pval, qval;
   double tol = 10e-14;

   if(!o_outbuf) {
      if((c_outbuf =(struct c_pkints *) malloc(maxbuf*sizeof(struct c_pkints)))
                                                                   ==NULL) {
         fprintf(stderr,"cannot allocate memory for c_outbuf in packit\n");
         exit(PSI_RETURN_FAILURE);
         }
      if((o_outbuf = (struct o_pkints *) malloc(maxbuf*sizeof(struct o_pkints)))
                                                                   ==NULL) {
         fprintf(stderr,"cannot allocate memory for o_outbuf in packit\n");
         exit(PSI_RETURN_FAILURE);
         }
      }

   tmpsiz = ioff[nbasis];

   if(gtmp==NULL && !twocon) {
      gtmp = (double *) init_array(tmpsiz);
      ptmp = (double *) init_array(tmpsiz);
      gtmpo = (double *) init_array(tmpsiz);
      ptmpo = (double *) init_array(tmpsiz);
      if(uhf){
	  ptmpo2 = (double *) init_array(tmpsiz);
	  gtmpo2 = (double *) init_array(tmpsiz);
	  /*  testj = (double *) init_array(ioff[tmpsiz]);
	      testk = (double *) init_array(ioff[tmpsiz]);*/
      }
      
      for(k=joff=0; k < num_ir ; k++) {
         if(nn=scf_info[k].num_so) {
            for(i=0; i < nn ; i++)
//	      printf("uhf = %d\n", uhf);
		   if(!uhf){
		       for(j=0; j <= i ; j++) {
			   ptmp[ioff[i+joff]+j+joff] = 
			       scf_info[k].pmat[ioff[i]+j];
			   ptmpo[ioff[i+joff]+j+joff] = 
			       scf_info[k].pmato[ioff[i]+j];
		       }
		   }
		   else{
		      for(j=0; j <= i ; j++) { 
			  ptmp[ioff[i+joff]+j+joff] =
			      scf_info[k].pmat[ioff[i]+j];
			  ptmpo[ioff[i+joff]+j+joff] =
			      -spin_info[0].scf_spin[k].pmat[ioff[i]+j];
			  ptmpo2[ioff[i+joff]+j+joff] =
			      -spin_info[1].scf_spin[k].pmat[ioff[i]+j];
		      }
		   }
	    joff += nn;
	 }
      }
   }
   
   if(!endflg) {
       for(i=0; i < cscf_nint ; i++) {
	   pval=pa[i];
	   qval=pb[i];
	   ij = lbij[i];
	   kl = lbkl[i];
	   if(print & 16) fprintf(outfile,"%5d%5d%9.5f%9.5f\n",
				  ij,kl,pval,qval);
	   if (fabs(pval) >= tol || fabs(qval) >= tol) {
	       if(!uhf){
		   if(qval) {
		       o_outbuf[iblo].ij = ij;
		       o_outbuf[iblo].kl = kl;
		       o_outbuf[iblo].pval = pval;
		       o_outbuf[iblo].qval = qval;
		       
		       if(!twocon) {
			   gtmp[ij] += ptmp[kl]*pval;
			   gtmp[kl] += ptmp[ij]*pval;
			   gtmpo[ij] += ptmpo[kl]*qval;
			   gtmpo[kl] += ptmpo[ij]*qval;
		       }
		       
		       iblo++;
		       
		       if (iblo >= maxbuf) {
			   readflgo=1;
			   psio_write(PKmat.unit, PKmat.key, (char *) o_outbuf, sizeof(struct o_pkints)*maxbuf,
				      PKmat.bufpos, &(PKmat.bufpos));
			   num_ints_o += iblo;
			   num_bufs_o++;
			   iblo=0;
		       }
		   }
		   else {
		       c_outbuf[iblc].ij = ij;
		       c_outbuf[iblc].kl = kl;
		       c_outbuf[iblc].pval = pval;
		       
		       if(!twocon) {
			   gtmp[ij] += ptmp[kl]*pval;
			   gtmp[kl] += ptmp[ij]*pval;
		       }
		       
		       iblc++;
		       
		       if (iblc >= maxbuf) {
			   readflgc=1;
			   psio_write(Pmat.unit, Pmat.key, (char *) c_outbuf, sizeof(struct c_pkints)*maxbuf,
				      Pmat.bufpos, &(Pmat.bufpos));
			   num_ints_c += iblc;
			   num_bufs_c++;
			   iblc=0;
		       }
		   }
	       }
	       else{
		   o_outbuf[iblo].ij = ij;
		   o_outbuf[iblo].kl = kl;
		   o_outbuf[iblo].pval = pval;
		   o_outbuf[iblo].qval = qval;
		   
		   /*testj[ioff[ij]+kl] = pval;
		     testk[ioff[ij]+kl] = qval;*/
		   
		   gtmp[ij] += ptmp[kl]*pval;
		   gtmpo[ij] += ptmpo[kl]*qval;
		   gtmpo2[ij] += ptmpo2[kl]*qval;
		   if(ij!=kl){
		       gtmp[kl] += ptmp[ij]*pval;
		       gtmpo[kl] += ptmpo[ij]*qval;
		       gtmpo2[kl] += ptmpo2[ij]*qval;
		       }
		   
		   iblo++;
		   
		   if (iblo >= maxbuf) {
		       readflgo=1;
		       psio_write(PKmat.unit, PKmat.key, (char *) o_outbuf, sizeof(struct o_pkints)*maxbuf,
				  PKmat.bufpos, &(PKmat.bufpos));
		       num_ints_o += iblo;
		       num_bufs_o++;
		       iblo=0;
		   }
	       }
	   }
       }
       cscf_nint=0;
   }
   else {
      num_ints_o += iblo;
      num_ints_c += iblc;
      num_bufs_c++;
      num_bufs_o++;
      fprintf(outfile,"\n%10d integrals written to file92 in %3d buffers\n",
                                       num_ints_o,num_bufs_o);
      fprintf(outfile,"%10d integrals written to file93 in %3d buffers\n",
                                       num_ints_c,num_bufs_c);
      lasto = iblo;
      lastc = iblc;
      if(readflgo)
	psio_write(PKmat.unit, PKmat.key, (char *) o_outbuf, sizeof(struct o_pkints)*iblo,
		   PKmat.bufpos, &(PKmat.bufpos));
      if(readflgc)
	psio_write(Pmat.unit, Pmat.key, (char *) c_outbuf, sizeof(struct c_pkints)*iblc,
		   Pmat.bufpos, &(Pmat.bufpos));
      

      /* testing stuff */

/*      for(i=0;i<nbasis;i++){
	  for(j=0;j<=i;j++){
	      for(k=0;k<=i;k++){
		  lmax = (k==i) ? j : k;
		  for(l=0;l<=lmax;l++){
		      ijkl = ioff[ioff[i]+j]+ioff[k]+l;
		      if(uhf) *//*fprintf(JK,"\n%5d%5d%5d%5d%5d%5d%5d\t%10.15lf\t%10.15lf",
			i,j,k,l,ioff[i]+j,ioff[k]+l
			,ijkl,testj[ijkl]
			,testk[ijkl]);*/
      /* fprintf(JK,"\n%5d%5d%5d%5d\t%10.15lf\t%10.15lf",
			      i,j,k,l,
			      testj[ijkl],testk[ijkl]);
		  }}}}
      */
      if(!twocon) {
         for(k=joff=0; k < num_ir ; k++) {
            if(nn=scf_info[k].num_so) {
               for(i=0; i < nn ; i++)
                  for(j=0; j <= i ; j++) {
		      if(!uhf){
			  scf_info[k].gmat[ioff[i]+j] = 
			      gtmp[ioff[i+joff]+j+joff];
			  scf_info[k].gmato[ioff[i]+j] = 
			      gtmpo[ioff[i+joff]+j+joff];
		      }
		      else{
			  spin_info[0].scf_spin[k].gmat[ioff[i]+j] =
			      gtmp[ioff[i+joff]+j+joff]
			      +gtmpo[ioff[i+joff]+j+joff];
			  spin_info[1].scf_spin[k].gmat[ioff[i]+j] =
			      gtmp[ioff[i+joff]+j+joff]
			      +gtmpo2[ioff[i+joff]+j+joff];
		      }
		  }
               joff += nn;
	    }
	 }
	 /* for(k=0;k<num_ir;k++){
	     for(i=0;i<2;i++){
	     print_array(spin_info[i].scf_spin[k].gmat,scf_info[k].num_so,gmat);
	     }
	     }*/

         free(gtmp);
         free(ptmp);
         free(gtmpo);
         free(ptmpo);
         gtmp = ptmp = gtmpo = ptmpo = NULL;
	 if(uhf){
	     free(gtmpo2);
		free(ptmpo2);
        gtmpo2 = ptmpo2 = NULL;
	 }
      }
   }
}

}} // namespace psi::cscf
