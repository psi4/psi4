/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.6  2004/05/03 04:32:40  crawdad
 * Major mods based on merge with stable psi-3-2-1 release.  Note that this
 * version has not been fully tested and some scf-optn test cases do not run
 * correctly beccause of changes in mid-March 2004 to optking.
 * -TDC
 *
 * Revision 1.5.8.1  2004/04/10 19:41:32  crawdad
 * Fixed the DIIS code for UHF cases.  The new version uses the Pulay scheme of
 * building the error vector in the AO basis as FDS-SDF, followed by xformation
 * into the orthogonal AO basis.   This code converges faster for test cases
 * like cc8, but fails for linearly dependent basis sets for unknown reasons.
 * -TDC
 *
 * Revision 1.5  2002/04/03 02:06:01  janssen
 * Finish changes to use new include paths for libraries.
 *
 * Revision 1.4  2000/12/05 19:40:03  sbrown
 * Added Unrestricted Kohn-Sham DFT.
 *
 * Revision 1.3  2000/06/22 22:15:00  evaleev
 * Modifications for KS DFT. Reading in XC Fock matrices and XC energy in formg_direct need to be uncommented (at present those are not produced by CINTS yet).
 *
 * Revision 1.2  2000/06/02 13:32:15  kenny
 *
 *
 * Added dynamic integral accuracy cutoffs for direct scf.  Added a few global
 * variables.  Added keyword 'dyn_acc'; true--use dynamic cutoffs.  Use of
 * 'dconv' and 'delta' to keep track of density convergence somewhat awkward,
 * but avoids problems when accuracy is switched and we have to wipe out density
 * matrices.  Also added error message and exit if direct rohf singlet is
 * attempted since it doesn't work.
 * --Joe Kenny
 *
 * Revision 1.1.1.1  2000/02/04 22:52:33  evaleev
 * Started PSI 3 repository
 *
 * Revision 1.2  1999/11/17 19:40:46  evaleev
 * Made all the adjustments necessary to have direct UHF working. Still doesn't work though..
 *
 * Revision 1.1  1999/11/02 23:55:56  sbrown
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
 * Revision 1.1.1.1  1999/04/12 16:59:25  evaleev
 * Added a version of CSCF that can work with CINTS.
 * -Ed
 *
 * Revision 1.1  1991/06/15  20:22:21  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: dmatuhf.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include <libpsio/psio.h>
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void dmatuhf()
{
  int i,j,k,l,ij,n,m,jj,kk;
  int max, off, ntri;
  int ndocc;
  int nn;
  double ptemp,ctmp;
  struct symm *s;
  struct spin *sp;
  double *dmat;
  double **cmat;
   
  for(m = 0; m < 2; m++){
    for (l=0; l < num_ir ; l++) {
      sp = &spin_info[m];
      if (n=scf_info[l].num_so) {
	ndocc = sp->scf_spin[l].noccup;
	for (i=ij=0; i < n ; i++ ) {
	  /*--------------------------------------

	  OFF DIAGONAL ELEMENTS

	  -------------------------------------*/

	  for (j=0; j <i; j++,ij++) {
	    ptemp=0.0;
	    for (k=0; k < ndocc ; k++)
	      ptemp += 2.0*sp->scf_spin[l].cmat[i][k]
		*sp->scf_spin[l].cmat[j][k];
	    sp->scf_spin[l].dpmat[ij] = ptemp - sp->scf_spin[l].pmat[ij];
	    sp->scf_spin[l].pmat[ij] = ptemp; 
	  }
	  /*----------------------------------
 
	  DIAGONAL ELEMENTS

	  --------------------------------*/
	  ptemp = 0.0;
	  for (k=0; k < ndocc ; k++) {
	    ctmp=sp->scf_spin[l].cmat[i][k];
	    ptemp += ctmp*ctmp;
	  }
	  sp->scf_spin[l].dpmat[ij] = ptemp - sp->scf_spin[l].pmat[ij];
	  sp->scf_spin[l].pmat[ij] = ptemp;
	  ij++;
	} 
	if(print & 4) {
	  fprintf(outfile,
		  "\nSpin case %d density matrix for irrep %s",m,scf_info[l].irrep_label);
	  print_array(spin_info[m].scf_spin[l].pmat,
		      scf_info[l].num_so,outfile);
	}
      }
    }
  }
   
  for(l=0;l < num_ir; l++){
    if(nn=scf_info[l].num_so) {
      for(ij=0;ij<ioff[nn];ij++) {
	ptemp = spin_info[0].scf_spin[l].pmat[ij] +
	  spin_info[1].scf_spin[l].pmat[ij];
	scf_info[l].dpmat[ij] = ptemp - scf_info[l].pmat[ij];
	scf_info[l].pmat[ij] = ptemp;
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

    psio_open(itapDSCF, PSIO_OPEN_NEW);
    psio_write_entry(itapDSCF, "Integrals cutoff", (char *) &eri_cutoff, sizeof(double));

    ntri = nbasis*(nbasis+1)/2;
    dmat = init_array(ntri);

    /*--- Get full dpmata ---*/
    for(i=0;i<num_ir;i++) {
      max = scf_info[i].num_so;
      off = scf_info[i].ideg;
      for(j=0;j<max;j++) {
	jj = j + off;
	for(k=0;k<=j;k++) {
	  kk = k + off;
	  if(acc_switch) {
	    dmat[ioff[jj]+kk] = spin_info[0].scf_spin[i].pmat[ioff[j]+k];
	    spin_info[0].scf_spin[i].dpmat[ioff[j]+k] = 0.0;
	  }
	  else {
	    dmat[ioff[jj]+kk] = spin_info[0].scf_spin[i].dpmat[ioff[j]+k];
	  }
	}
      }
    }
    psio_write_entry(itapDSCF, "Difference Alpha Density", (char *) dmat, sizeof(double)*ntri);

    /*--- Get full dpmatb ---*/
    for(i=0;i<num_ir;i++) {
      max = scf_info[i].num_so;
      off = scf_info[i].ideg;
      for(j=0;j<max;j++) {
	jj = j + off;
	for(k=0;k<=j;k++) {
	  kk = k + off;
	  if(acc_switch) {
	    dmat[ioff[jj]+kk] = spin_info[1].scf_spin[i].pmat[ioff[j]+k];
	    spin_info[1].scf_spin[i].dpmat[ioff[j]+k] = 0.0;
	  }
	  else {
	    dmat[ioff[jj]+kk] = spin_info[1].scf_spin[i].dpmat[ioff[j]+k];
	  }
	}
      }
    }
    psio_write_entry(itapDSCF, "Difference Beta Density", (char *) dmat, sizeof(double)*ntri);
     
    if (ksdft) {
      /* ---- Get Occupied Eigenvector matrix for DFT ----*/
	 
      /*------ Get the Alpha Occupied Eigenvector --------*/
      ntri = nbfso*a_elec;
      cmat = block_matrix(nbfso,a_elec);
      for(i=l=0;i<num_ir;i++){
	max = spin_info[0].scf_spin[i].nclosed;
	off = scf_info[i].ideg;
	for(j=0;j<max;j++){
	  for(k=0;k<scf_info[i].num_mo;k++) {
	    kk = k + off;
	    cmat[kk][l] = spin_info[0].scf_spin[i].cmat[k][j];
	  }
	  l++;
	}
      }
      /*fprintf(outfile,"\na_elec = %d",a_elec);
	fprintf(outfile,"\nAlpha Occupied Eigenvector from CSCF");
	print_mat(cmat,nbfso,a_elec,outfile);*/
	 
      psio_write_entry(itapDSCF, "Number of MOs", 
		       (char *) &(nmo),sizeof(int)); 
	
      psio_write_entry(itapDSCF, "Number of Alpha DOCC",
		       (char *) &(a_elec),sizeof(int));
       
      psio_write_entry(itapDSCF, "Alpha Occupied SCF Eigenvector", 
		       (char *) &(cmat[0][0]),sizeof(double)*ntri);
      /*fprintf(outfile,"\na_elec = %d",a_elec);*/
      free_block(cmat);
      cmat = NULL;
      /*------ Get the Alpha Occupied Eigenvector --------*/
      ntri = nbfso*b_elec;
      cmat = block_matrix(nbfso,b_elec);
      for(i=l=0;i<num_ir;i++){
	max = spin_info[1].scf_spin[i].nclosed;
	off = scf_info[i].ideg;
	for(j=0;j<max;j++){
	  for(k=0;k<scf_info[i].num_mo;k++) {
	    kk = k + off;
	    cmat[kk][l] = spin_info[1].scf_spin[i].cmat[k][j];
	  }
	  l++;
	}
      }
      /*fprintf(outfile,"\nBeta Occupied Eigenvector from CSCF");
	print_mat(cmat,nbfso,b_elec,outfile);*/

      psio_write_entry(itapDSCF, "Number of Beta DOCC",
		       (char *) &(b_elec),sizeof(int));
      psio_write_entry(itapDSCF, "Beta Occupied SCF Eigenvector", 
		       (char *) &(cmat[0][0]),sizeof(double)*ntri);
      free_block(cmat);
	 
      /*--- Get full dpmata ---*/
      for(i=0;i<num_ir;i++) {
	max = scf_info[i].num_so;
	off = scf_info[i].ideg;
	for(j=0;j<max;j++) {
	  jj = j + off;
	  for(k=0;k<=j;k++) {
	    kk = k + off;
	    dmat[ioff[jj]+kk] = spin_info[0].scf_spin[i].pmat[ioff[j]+k];
	  }
	}
      }
      psio_write_entry(itapDSCF, "Total Alpha Density", (char *) dmat, sizeof(double)*ntri);

      /*--- Get full dpmatb ---*/
      for(i=0;i<num_ir;i++) {
	max = scf_info[i].num_so;
	off = scf_info[i].ideg;
	for(j=0;j<max;j++) {
	  jj = j + off;
	  for(k=0;k<=j;k++) {
	    kk = k + off;
	    dmat[ioff[jj]+kk] = spin_info[1].scf_spin[i].pmat[ioff[j]+k];
	  }
	}
      }
      psio_write_entry(itapDSCF, "Total Beta Density", (char *) dmat, sizeof(double)*ntri);
    }
    free(dmat);
    psio_close(itapDSCF, 1);
  }

  return;
}

}} // namespace psi::cscf
