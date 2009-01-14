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
 * Revision 1.5.4.3  2004/04/10 19:41:32  crawdad
 * Fixed the DIIS code for UHF cases.  The new version uses the Pulay scheme of
 * building the error vector in the AO basis as FDS-SDF, followed by xformation
 * into the orthogonal AO basis.   This code converges faster for test cases
 * like cc8, but fails for linearly dependent basis sets for unknown reasons.
 * -TDC
 *
 * Revision 1.5.4.2  2004/04/09 00:17:37  evaleev
 * Corrected dimensions of the matrix.
 *
 * Revision 1.5.4.1  2004/04/07 03:23:32  crawdad
 * Working to fix UHF-based DIIS.
 * -TDC
 *
 * Revision 1.5  2002/12/06 15:50:32  crawdad
 * Changed all exit values to PSI_RETURN_SUCCESS or PSI_RETURN_FAILURE as
 * necessary.  This is new for the PSI3 execution driver.
 * -TDC
 *
 * Revision 1.4  2000/12/05 19:40:04  sbrown
 * Added Unrestricted Kohn-Sham DFT.
 *
 * Revision 1.3  2000/10/13 19:51:22  evaleev
 * Cleaned up a lot of stuff in order to get CSCF working with the new "Mo-projection-capable" INPUT.
 *
 * Revision 1.2  2000/06/22 22:15:02  evaleev
 * Modifications for KS DFT. Reading in XC Fock matrices and XC energy in formg_direct need to be uncommented (at present those are not produced by CINTS yet).
 *
 * Revision 1.1.1.1  2000/02/04 22:52:34  evaleev
 * Started PSI 3 repository
 *
 * Revision 1.3  1999/11/17 19:40:47  evaleev
 * Made all the adjustments necessary to have direct UHF working. Still doesn't work though..
 *
 * Revision 1.2  1999/11/04 19:24:31  localpsi
 * STB (11/4/99) - Added the orb_mix feature which is equivalent to guess = mix
 * in G94 and also fixed restarting so that if you have different wavefuntions,
 * everything works.  Also if you specify no DOCC and SOCC and restart, if the
 * wavefunctions are different, it will guess again.
 *
 * Revision 1.1  1999/11/02 23:56:00  localpsi
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
 * Revision 1.1.1.1  1999/04/12 16:59:28  evaleev
 * Added a version of CSCF that can work with CINTS.
 * -Ed
 *
 * Revision 1.5  1998/06/30  14:11:12  sbrown
 * *************************************************************
 * *Program Modification                                       *
 * *By: Shawn Brown                                            *
 * *Date: June 30, 1998                                        *
 * *Altered program to make a guess at occupations from the    *
 * *diagonalized core hamiltonian matrix.  Program can now     *
 * *make a guess at the beginning of the calculation or at     *
 * *or at every iteration.  Use the latter at your own risk.   *
 * *See man pages for details on new keywords.                 *
 * *************************************************************
 *
 * Revision 1.4  1995/07/21  17:37:15  psi
 * Made Jan 1st 1995 cscf the current accepted version of cscf.  Some
 * unidentified changes made after that date were causing problems.
 *
 * Revision 1.1  1991/06/15  20:22:40  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: uhf_iter.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

int uhf_iter()
{
  int i,j,l,m,t,ij;
  int nn,num_mo,newci;
  double cimax;
  double **scr;
  double **fock_c;
  double **fock_ct;
  double **ctrans;
  double tol = 1.0e-14;
  struct symm *s;
  struct spin *sp;

  diiser = 0.0;
  scr = block_matrix(nsfmax,nsfmax);
  fock_c = block_matrix(nsfmax,nsfmax);
  fock_ct = block_matrix(nsfmax,nsfmax);
  ctrans = block_matrix(nsfmax,nsfmax);
   
  /* and iterate */

  for (iter=0; iter < itmax ; ) {
    for(t = 0; t<2 ;t++){
      sp = &spin_info[t];

      if(print & 4) {
	for(m=0; m < num_ir ; m++) {
	  if (nn=scf_info[m].num_so) {
	    fprintf(outfile,
		    "\n%s gmat for irrep %s",sp->spinlabel,scf_info[m].irrep_label);
	    print_array(sp->scf_spin[m].gmat,nn,outfile);
	    if (ksdft) {
	      fprintf(outfile,
		      "\n%s xcmat for irrep %s",sp->spinlabel,scf_info[m].irrep_label);
	      print_array(sp->scf_spin[m].xcmat,nn,outfile);
	    }
	  }
	}
      }
	   
      for (m=0; m < num_ir ; m++) {
	s = &scf_info[m];
	       
	if (nn=s->num_so) {

	  /*		   
	  if(!m && !t) {
	    fprintf(outfile, "Fock matrix in top of uhf_iter:\n");
	    print_array(sp->scf_spin[m].fock_pac, nn, outfile);
	  }
	  */

	  /*  form fock matrix = h+g */
	  add_arr(s->hmat,sp->scf_spin[m].gmat,
		  sp->scf_spin[m].fock_pac,ioff[nn]);

	}
      }
    }

    /*----------------------------------------------------
      In KS DFT case, Fock matrix doesn't include Fxc yet
      add them up only after the computation of energy
      ----------------------------------------------------*/
    ecalc(tol);
    if (ksdft) {
      /* now form f = h + g + fxc */
      /* it should be alright to use fock_pac as 2 arguments */
      for(t = 0; t<2 ;t++){
	sp = &spin_info[t];
	for (m=0; m < num_ir ; m++) {
	  s = &(sp->scf_spin[m]);
	  if (nn=scf_info[m].num_so)
	    add_arr(s->fock_pac,s->xcmat,s->fock_pac,ioff[nn]);
	  if(print & 4) {
	    fprintf(outfile,"\n%s fock for irrep %s"
		    ,sp->spinlabel,scf_info[m].irrep_label);
	    print_array(sp->scf_spin[m].fock_pac,nn,outfile);
	  }
	}
      }
    }
       
    /* create new fock matrix in fock_pac or fock_eff */
    if(!diisflg) diis_uhf();
       
    for(t=0;t<2;t++){
      sp = &spin_info[t];
      for (m=0; m < num_ir ; m++) {
	s = &scf_info[m];
	if (nn=s->num_so) {
	  num_mo = s->num_mo;
	
	  /*
	  if(!m && !t) {
	    fprintf(outfile, "\nuhf_iter Fock matrix irrep %d spin %d iter %d\n", m, t, iter);
	    print_array(sp->scf_spin[m].fock_pac, nn, outfile);
	  }
	  */
	   
	  /* transform fock_pac to mo basis */
	  tri_to_sq(sp->scf_spin[m].fock_pac,fock_ct,nn);
	  /*		   mxmb(sp->scf_spin[m].cmat,nn,1
			   ,fock_ct,1,nn,scr,1,nn,nn,nn,nn);
			   mxmb(scr,1,nn,sp->scf_spin[m].cmat,1
			   ,nn,fock_c,1,nn,nn,nn,nn);*/

	  mmult(sp->scf_spin[m].cmat,1,fock_ct,0,scr,0,num_mo,nn,nn,0);
	  mmult(scr,0,sp->scf_spin[m].cmat,0,fock_c,0,num_mo,nn,num_mo,0);

	  /*
	  if(!m && !t) {
	    fprintf(outfile, "\nMO uhf_iter Fock matrix irrep %d spin %d iter %d\n", m, t, iter);
	    print_mat(fock_c, num_mo, num_mo, outfile);
	  }
	  */
		   
	  /*  diagonalize fock_c to get ctrans */
	  sq_rsp(num_mo,num_mo,fock_c,sp->scf_spin[m].fock_evals
		 ,1,ctrans,tol);
		   
	  if(print & 4) {
	    fprintf(outfile,"\n %s eigenvector for irrep %s\n"
		    ,sp->spinlabel,s->irrep_label);
	    eivout(ctrans,sp->scf_spin[m].fock_evals,num_mo,num_mo,outfile);
	  }
		   
	  /*		   mxmb(sp->scf_spin[m].cmat,1,nn,
			   ctrans,1,nn,scr,1,nn,nn,nn,nn);*/
	  mmult(sp->scf_spin[m].cmat,0,ctrans,0,scr,0,nn,num_mo,num_mo,0);
				   
	  if(print & 4) {
	    fprintf(outfile,"\n %s eigenvector after irrep %s\n",
		    sp->spinlabel,s->irrep_label);
	    print_mat(scr,nn,num_mo,outfile);
	  }
		   
	  for (i=0; i < nn; i++)
	    for (j=0; j < num_mo; j++)
	      sp->scf_spin[m].cmat[i][j] = scr[i][j];

	}
      }
	   
      if(converged) {
        free_block(scr);
        free_block(fock_c);
        free_block(fock_ct);
        free_block(ctrans);
        scr = fock_c = fock_ct = ctrans = NULL;
        return 1;
      }
    }
    schmit_uhf(1);
       
    if(print & 4) {
      for(j=0; j < 2; j++){
	for(i=0; i < num_ir ; i++) {
	  s = &scf_info[i];
	  if (nn=s->num_so) {
	    num_mo = s->num_mo;
	    fprintf(outfile,"\northogonalized mos irrep %s\n",
		    s->irrep_label);
	    print_mat(spin_info[j].scf_spin[i].cmat,nn,num_mo,outfile);
	  }
	}
      }
    }
       
    if(mixing && iter ==1)
      orb_mix();

    /* form new density matrix */
    dmatuhf();
       
    /* and form new fock matrix */
    if(iter < itmax) {
      if (!direct_scf)
	formg_open();
      else
	formg_direct();
    }
  }
  return 0;
}

}} // namespace psi::cscf
