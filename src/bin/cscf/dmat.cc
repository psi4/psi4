/*! \file 
    \ingroup (CSCF)
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.8  2002/04/03 02:06:01  janssen
 * Finish changes to use new include paths for libraries.
 *
/* Revision 1.7  2000/12/05 19:40:03  sbrown
/* Added Unrestricted Kohn-Sham DFT.
/*
/* Revision 1.6  2000/07/10 18:03:30  sbrown
/* Enabling cscf to send over just the occupied SCF eigenvector for DFT
/* calculations.  Only done for the RHF case.
/*
/* Revision 1.5  2000/07/06 20:04:01  sbrown
/* Added capabilities to send the eigenvector to cints for DFT
/* calculations.
/*
/* Revision 1.4  2000/06/26 19:04:09  sbrown
/* Added DFT capapbilities to interface with cints using direct scf
/*
/* Revision 1.3  2000/06/22 22:14:59  evaleev
/* Modifications for KS DFT. Reading in XC Fock matrices and XC energy in formg_direct need to be uncommented (at present those are not produced by CINTS yet).
/*
/* Revision 1.2  2000/06/02 13:32:15  kenny
/*
/*
/* Added dynamic integral accuracy cutoffs for direct scf.  Added a few global
/* variables.  Added keyword 'dyn_acc'; true--use dynamic cutoffs.  Use of
/* 'dconv' and 'delta' to keep track of density convergence somewhat awkward,
/* but avoids problems when accuracy is switched and we have to wipe out density
/* matrices.  Also added error message and exit if direct rohf singlet is
/* attempted since it doesn't work.
/* --Joe Kenny
/*
/* Revision 1.1.1.1  2000/02/04 22:52:29  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.6  1999/11/17 19:40:45  evaleev
/* Made all the adjustments necessary to have direct UHF working. Still doesn't work though..
/*
/* Revision 1.5  1999/11/02 23:23:42  evaleev
/* Commented out a line in dmat..
/*
/* Revision 1.4  1999/11/02 18:10:13  evaleev
/* Direct SCF improved
/*
/* Revision 1.3  1999/10/22 19:47:18  evaleev
/* A direct SCF-enabled version (set DIRECT_SCF=TRUE in input.dat).
/*
/* Revision 1.2  1999/08/17 19:04:14  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:25  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
/*
 * Revision 1.1  1991/06/15  20:22:21  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: dmat.cc 3592 2007-09-28 13:01:33Z evaleev $";

#define EXTERN
#include <libpsio/psio.h>
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void dmat()
{
   int i,j,k,l,ij,jj,kk,nn;
   int max, off, ntri, joff;
   int ndocc,nsocc,nhocc;
   double ptempc,ptempo,ctmp;
   double **cmat;
   double *dmat;
   extern double delta;
   struct symm *s;

   for (l=0; l < num_ir ; l++) {
      s = &scf_info[l];
      if (nn=s->num_so) {
         ndocc = s->nclosed;
         nsocc = s->nopen;
         nhocc = 0;
         if(s->nhalf) nhocc=1;
            
         for (i=ij=0; i < nn ; i++ ) {
            for (j=0; j < i; j++,ij++) {
               ptempc=ptempo=0.0;
               for (k=0; k < ndocc ; k++)
                  ptempc += 4.0*s->cmat[i][k]*s->cmat[j][k];

               for (k=ndocc; k < ndocc+nsocc ; k++)
                  ptempo += 2.0*s->occ_num[k]*s->cmat[i][k]*s->cmat[j][k];

               for (k=ndocc+nsocc; k < ndocc+nsocc+nhocc ; k++)
                  ptempo += 2.0*s->occ_num[k]*s->cmat[i][k]*s->cmat[j][k];

               if(iopen) {
                  s->dpmato[ij] = ptempo - s->pmato[ij];
                  s->pmato[ij] = ptempo;
                  }
               s->dpmat[ij] = ptempc+ptempo - s->pmat[ij];
               s->pmat[ij] = ptempc+ptempo;
               }
            ptempc=ptempo=0.0;
            for (k=0; k < ndocc ; k++) {
               ctmp=s->cmat[i][k];
               ptempc += 2.0*ctmp*ctmp;
               }
            for (k=ndocc; k < ndocc+nsocc ; k++) {
               ctmp=s->cmat[i][k];
               ptempo += ctmp*ctmp*s->occ_num[k];
               }
            for (k=ndocc+nsocc; k < ndocc+nsocc+nhocc ; k++) {
               ctmp=s->cmat[i][k];
               ptempo += ctmp*ctmp*s->occ_num[k];
               }
            if(iopen) {
               s->dpmato[ij] = ptempo - s->pmato[ij];
               s->pmato[ij] = ptempo;
               }
            s->dpmat[ij] = ptempc+ptempo - s->pmat[ij];
            s->pmat[ij] = ptempc+ptempo;
            ij++;
            }

         if(print & 4) {
            fprintf(outfile,
                       "\ntotal density matrix for irrep %s",s->irrep_label);
            print_array(s->pmat,nn,outfile);
            print_array(s->dpmat,nn,outfile);
            if(iopen) {
               fprintf(outfile,"\nopen-shell density matrix for irrep %s",
                                                              s->irrep_label);
               print_array(s->pmato,nn,outfile);
               print_array(s->dpmato,nn,outfile);
               }
            }
         }
      }

   /*-----------------------
     Handle direct SCF here
    -----------------------*/
   if(direct_scf) {
     /*decide what accuracy to request for direct_scf*/
     if (dyn_acc) {
       if((iter<30)&&(tight_ints==0)&&(delta>1.0E-5)) {
          eri_cutoff=1.0E-6;
       }
 
       if((tight_ints==0)&&(delta<=1.0E-5)){
          fprintf(outfile,"  Switching to full integral accuracy\n");
          acc_switch=1;
          tight_ints=1;
          eri_cutoff=1.0E-14;
       }
     }
     
     psio_open(itapDSCF,PSIO_OPEN_NEW);
     psio_write_entry(itapDSCF, "Integrals cutoff", (char *) &eri_cutoff,
		      sizeof(double));
     
     ntri = nbasis*(nbasis+1)/2;

     /*--- Get full dpmat ---*/
     dmat = init_array(ntri);
     for(i=0;i<num_ir;i++) {
       max = scf_info[i].num_so;
       off = scf_info[i].ideg;
       for(j=0;j<max;j++) {
	 jj = j + off;
	 for(k=0;k<=j;k++) {
	   kk = k + off;
	   if(acc_switch) {
	     dmat[ioff[jj]+kk] = scf_info[i].pmat[ioff[j]+k];
	     scf_info[i].dpmat[ioff[j]+k] = 0.0;
	   }
	   else
	     dmat[ioff[jj]+kk] = scf_info[i].dpmat[ioff[j]+k];
	 }
       }
     }
     psio_write_entry(itapDSCF, "Difference Density", (char *) dmat,
		      sizeof(double)*ntri);

     /* ----- Get Occupied Eigenvector matrix for DFT ---- */
     if (ksdft) {
	 ntri = nbfso*n_closed;
	 cmat = block_matrix(nbfso,n_closed);
	 for(i=l=0;i<num_ir;i++){
	     max = scf_info[i].nclosed;
	     off = scf_info[i].ideg;
	     for(j=0;j<max;j++){
		 for(k=0;k<scf_info[i].num_mo;k++) {
		     kk = k + off;
		     cmat[kk][l] = scf_info[i].cmat[k][j];
		 }
		 l++;
	     }
	 }
	 
	 /*fprintf(outfile,"\nOccupied Eigenvector from CSCF");
	   print_mat(cmat,nbfso,n_closed,outfile);*/
	 
	 psio_write_entry(itapDSCF, "Number of MOs", 
			  (char *) &(nmo),sizeof(int));   
	 psio_write_entry(itapDSCF, "Number of DOCC",
			  (char *) &(n_closed),sizeof(int));
	 psio_write_entry(itapDSCF, "Occupied SCF Eigenvector", 
			  (char *) &(cmat[0][0]),sizeof(double)*ntri);
	 free_block(cmat);
	 
	 /*--- Get full dpmat ---*/
	 for(i=0;i<num_ir;i++) {
	     max = scf_info[i].num_so;
	     off = scf_info[i].ideg;
	     for(j=0;j<max;j++) {
		 jj = j + off;
		 for(k=0;k<=j;k++) {
		     kk = k + off;
		     dmat[ioff[jj]+kk] = scf_info[i].pmat[ioff[j]+k];
		 }
	     }
	 }
	 psio_write_entry(itapDSCF, "Total Density", 
			  (char *) dmat, sizeof(double)*ntri);
     }
   

     /*--- Get full dpmato ---*/
     if (iopen) {
       for(i=0;i<num_ir;i++) {
	 max = scf_info[i].num_so;
	 off = scf_info[i].ideg;
	 for(j=0;j<max;j++) {
	   jj = j + off;
	   for(k=0;k<=j;k++) {
	     kk = k + off;
	     if(acc_switch) {
	       dmat[ioff[jj]+kk] = scf_info[i].pmato[ioff[j]+k];
	       scf_info[i].dpmato[ioff[j]+k] = 0.0;
	     }
	     else
	       dmat[ioff[jj]+kk] = scf_info[i].dpmato[ioff[j]+k];
	   }
	 }
       }
       psio_write_entry(itapDSCF, "Difference Open-Shell Density", 
			(char *) dmat,sizeof(double)*ntri);
     }
     free(dmat);
     psio_close(itapDSCF, 1);
   }
   
   return;
}                                            
                                                                                
}} // namespace psi::cscf
