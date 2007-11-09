/*! \file 
    \ingroup (CSCF)
    \brief Enter brief description of file here 
*/
static char *rcsid = "$Id: formg_direct.cc 3592 2007-09-28 13:01:33Z evaleev $";

#define EXTERN
#include <libpsio/psio.h>
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void formg_direct()
{
   double *gtmp;
   int stat;
   int i,j,k,jj,kk,l,off,joff,nn,max;      
   int ntri = ioff[nbasis];
   double den_percent;

   /*-----------------------------------
     Call CINTS to do all the dirty job
     Check if it ran fine
    -----------------------------------*/
   if(exitflag)
       exit(PSI_RETURN_FAILURE);
   
   stat = system("cints --fock");
   
   switch (stat) {
   case 0:
       /* everything is OK */
       break;
       
   default:
       /* Something went wrong */
       fprintf(outfile,"  formg_direct: System call to CINTS failed. Check to see if it's in your PATH\n");
       fprintf(stderr,"System call to CINTS failed. Check to see if it's in your PATH.\n");
       exit(PSI_RETURN_FAILURE);
   }
   
   gtmp = (double *) init_array(ntri);

   /*----------------------------------------------------
     Zero-out ERI-dependent portions of the Fock matrix
     if ERI accuracy has just been switched (see dmat.c)
    ----------------------------------------------------*/
   if(acc_switch) {
       if (uhf) {
	   for(k=0; k < num_ir ; k++) {
               if(nn=scf_info[k].num_so) {
		   for(i=0; i < nn ; i++) {
		       for(j=0; j <= i ; j++) {
			   spin_info[0].scf_spin[k].gmat[ioff[i]+j] = 0.0;
			   spin_info[1].scf_spin[k].gmat[ioff[i]+j] = 0.0;
		       }
		   }
               }
	   }
       }
       else if (iopen) {
	   for(k=joff=0; k < num_ir ; k++) {
               if(nn=scf_info[k].num_so) {
		   for(i=0; i < nn ; i++) {
		       for(j=0; j <= i ; j++) {
			   scf_info[k].gmat[ioff[i]+j] = 0.0;
			   scf_info[k].gmato[ioff[i]+j] = 0.0;
		       }
		   }
               }
	   }
       }
       else {
	   for(i=0;i<num_ir;i++) {
               max = scf_info[i].num_so;
               for(j=0;j<max;j++) {
		   for(k=0;k<=j;k++) {
		       scf_info[i].gmat[ioff[j]+k] = 0.0;
		   }
               }
	   }
       }
   }
   
   psio_open(itapDSCF, PSIO_OPEN_OLD);
   if (ksdft) {
     psio_read_entry(itapDSCF, "DFT XC-energy", 
		     (char *) &(exc), sizeof(double));
    psio_read_entry(itapDSCF, "DFT X-energy", 
		    (char *) &(exch_energy), sizeof(double));
    psio_read_entry(itapDSCF, "DFT C-energy", 
		    (char *) &(corr_energy), sizeof(double));
    psio_read_entry(itapDSCF, "DFT Den", 
		    (char *) &(den_trace), sizeof(double));
    
    /* Check for consistency of the numrical integration */
    
    den_percent = 100*fabs((double) nelec-den_trace)/(double) nelec;
    
    if(den_percent > 1E-3){
	fprintf(outfile,"\nWarning, Numerically Integrated Density Inaccurate");
	fprintf(outfile,"\nNum of Electrons = %d",nelec);
	fprintf(outfile,"\nNum from Trace   = %5.5lf",den_trace);
	fprintf(outfile,"\nAbsolute Error   = %%%5.5lf\n\n",den_percent);
    }
   }
   
   
   if (uhf) {
     psio_read_entry(itapDSCF, "Alpha JX G-matrix", (char *) gtmp, sizeof(double)*ntri);
     for(k=joff=0; k < num_ir ; k++) {
       if(nn=scf_info[k].num_so) {
	 for(i=0; i < nn ; i++) {
	   for(j=0; j <= i ; j++) {
	     spin_info[0].scf_spin[k].gmat[ioff[i]+j] += gtmp[ioff[i+joff]+j+joff];
	   }
	 }
       }
       joff += nn;
     }
     
     psio_read_entry(itapDSCF, "Beta JX G-matrix", (char *) gtmp, sizeof(double)*ntri);
     for(k=joff=0; k < num_ir ; k++) {
       if(nn=scf_info[k].num_so) {
	 for(i=0; i < nn ; i++) {
	   for(j=0; j <= i ; j++) {
	     spin_info[1].scf_spin[k].gmat[ioff[i]+j] += gtmp[ioff[i+joff]+j+joff];
	   }
	 }
       }
       joff += nn;
     }

     if (ksdft) {
	 psio_read_entry(itapDSCF, "Alpha XC G-matrix", (char *) gtmp, sizeof(double)*ntri);
	 for(k=joff=0; k < num_ir ; k++) {
	     if(nn=scf_info[k].num_so) {
		 for(i=0; i < nn ; i++) {
		     for(j=0; j <= i ; j++) {
			 spin_info[0].scf_spin[k].xcmat[ioff[i]+j] = gtmp[ioff[i+joff]+j+joff];
		     }
		 }
	     }
	     joff += nn;
	 }
	 
	 psio_read_entry(itapDSCF, "Beta XC G-matrix", (char *) gtmp, sizeof(double)*ntri);
	 for(k=joff=0; k < num_ir ; k++) {
	     if(nn=scf_info[k].num_so) {
		 for(i=0; i < nn ; i++) {
		     for(j=0; j <= i ; j++) {
			 spin_info[1].scf_spin[k].xcmat[ioff[i]+j] = gtmp[ioff[i+joff]+j+joff];
		     }
		 }
	     }
	     joff += nn;
	 }
     }
   }
   else {
     psio_read_entry(itapDSCF, "Total JX G-matrix", (char *) gtmp, sizeof(double)*ntri);
     for(k=joff=0; k < num_ir ; k++) {
       if(nn=scf_info[k].num_so) {
	 for(i=0; i < nn ; i++) {
	   for(j=0; j <= i ; j++) {
	     scf_info[k].gmat[ioff[i]+j] += gtmp[ioff[i+joff]+j+joff];
	   }
	 }
       }
       joff += nn;
     }

     if (ksdft) {
	 psio_read_entry(itapDSCF, "Total XC G-matrix", (char *) gtmp, sizeof(double)*ntri);
	 for(k=joff=0; k < num_ir ; k++) {
	     if(nn=scf_info[k].num_so) {
		 for(i=0; i < nn ; i++) {
		     for(j=0; j <= i ; j++) {
			 scf_info[k].xcmat[ioff[i]+j] = gtmp[ioff[i+joff]+j+joff];
	     }
		 }
	     }
	     joff += nn;
	 }
     }

     if (iopen) {
       psio_read_entry(itapDSCF, "Open-shell JX G-matrix", (char *) gtmp, sizeof(double)*ntri);
       for(k=joff=0; k < num_ir ; k++) {
	 if(nn=scf_info[k].num_so) {
	   for(i=0; i < nn ; i++) {
	     for(j=0; j <= i ; j++) {
	       scf_info[k].gmato[ioff[i]+j] += gtmp[ioff[i+joff]+j+joff];
	     }
	   }
	 }
	 joff += nn;
       }
     }
   }
   psio_close(itapDSCF, 1);
   
   free(gtmp);
   return;
}

}} // namespace psi::cscf
