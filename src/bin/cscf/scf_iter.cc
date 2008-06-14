/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.11  2007/04/05 15:45:25  crawdad
 * Fixed a few memory leaks identified by valgrind. -TDC
 *
/* Revision 1.10  2004/05/03 04:32:40  crawdad
/* Major mods based on merge with stable psi-3-2-1 release.  Note that this
/* version has not been fully tested and some scf-optn test cases do not run
/* correctly beccause of changes in mid-March 2004 to optking.
/* -TDC
/*
/* Revision 1.9.8.2  2004/04/21 15:45:07  evaleev
/* Modified DIIS algorithm for RHF and ROHF to work in OSO basis rather than in
/* AO basis, to avoid difficulties of transforming between MO and AO bases
/* when linear dependencies are present.
/*
/* Revision 1.9.8.1  2004/04/06 21:29:05  crawdad
/* Corrections to the RHF/ROHF DIIS algorithm, which was simply incorrect.
/* The backtransformation of the DIIS error vectors to the AO basis was not
/* mathematically right.
/* -TDC and EFV
/*
/* Revision 1.9  2002/04/03 02:06:01  janssen
/* Finish changes to use new include paths for libraries.
/*
/* Revision 1.8  2000/07/10 18:03:34  sbrown
/* Enabling cscf to send over just the occupied SCF eigenvector for DFT
/* calculations.  Only done for the RHF case.
/*
/* Revision 1.7  2000/07/06 20:04:02  sbrown
/* Added capabilities to send the eigenvector to cints for DFT
/* calculations.
/*
/* Revision 1.6  2000/07/05 21:47:31  sbrown
/* Enabled the code to export the SCF eigenvector to CINTS when doing DFT.
/*
/* Revision 1.5  2000/06/27 21:12:33  evaleev
/* .
/*
/* Revision 1.4  2000/06/27 21:08:11  evaleev
/* Fixed a minor string manipulation problem in scf_input.c
/*
/* Revision 1.3  2000/06/26 19:04:12  sbrown
/* Added DFT capapbilities to interface with cints using direct scf
/*
/* Revision 1.2  2000/06/22 22:15:02  evaleev
/* Modifications for KS DFT. Reading in XC Fock matrices and XC energy in formg_direct need to be uncommented (at present those are not produced by CINTS yet).
/*
/* Revision 1.1.1.1  2000/02/04 22:52:32  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.6  1999/11/04 19:24:31  localpsi
/* STB (11/4/99) - Added the orb_mix feature which is equivalent to guess = mix
/* in G94 and also fixed restarting so that if you have different wavefuntions,
/* everything works.  Also if you specify no DOCC and SOCC and restart, if the
/* wavefunctions are different, it will guess again.
/*
/* Revision 1.5  1999/11/02 23:56:00  localpsi
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
/* Revision 1.4  1999/10/22 19:47:19  evaleev
/* A direct SCF-enabled version (set DIRECT_SCF=TRUE in input.dat).
/*
/* Revision 1.3  1999/10/11 17:03:18  evaleev
/* Modified the location of nmo in mconst array in file 30.
/*
/* Revision 1.2  1999/08/17 19:04:17  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The program prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes have to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:28  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id: scf_iter.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <libpsio/psio.h>

namespace psi { namespace cscf {

void scf_iter()

{
  int i,j,k,m,ij;
  int max,off,jj,kk;
  int ntri;
  int nn,num_mo,newci;
  double cimax=0;
  double occi, occj, occ0, den;
  double **scr;
  double **cmat;
  double *c1, *c2;
  double **fock_c, **fock_o;
  double **fock_ct;
  double **ctrans;
  double tol = 1.0e-14;
  struct symm *s;

  diiser = 0.0;
  scr = (double **) block_matrix(nsfmax,nsfmax);
  fock_c = (double **) block_matrix(nsfmax,nsfmax);
  fock_ct = (double **) block_matrix(nsfmax,nsfmax);
  ctrans = (double **) block_matrix(nsfmax,nsfmax);
  if(iopen) {
    c1 = (double *) init_array(num_ir);
    c2 = (double *) init_array(num_ir);
    fock_o = (double **) block_matrix(nsfmax,nsfmax);
  }

  for (i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      if (nn=s->num_so) {
	num_mo = s->num_mo;
	s->ucmat = block_matrix(num_mo, num_mo);
	for(j=0; j<num_mo; j++)
	  s->ucmat[j][j] = 1.0;
      }
  }

  /* set up c1 and c2  */
  if(iopen) {
    for (i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      if (nn=s->num_so) {
	num_mo = s->num_mo;
	for (j=0; j < num_mo ; j++) {
	  if(s->occ_num[j]==2.0) c1[i]=2.0;
	  if(s->occ_num[j]==1.0) c2[i]=1.0;
	  if(s->occ_num[j]==0.5) c2[i]=0.5;
	  if(s->occ_num[j]==1.5) c2[i]=1.5;
	}
	den = c1[i]-c2[i];
	if(den) {
	  c1[i] /= den;
	  c2[i] /= den;
	}
      }
    }
  }

  /* and iterate */

  for (iter=0; iter < itmax ; ) {
       
    if(print & 4) {
      for(m=0; m < num_ir ; m++) {
	s = &scf_info[m];
	if (nn=s->num_so) {
	  fprintf(outfile,
		  "\ngmat for irrep %s",s->irrep_label);
	  print_array(s->gmat,nn,outfile);
	  if(iopen) {
	    fprintf(outfile,
		    "\ngmato for irrep %s",s->irrep_label);
	    print_array(s->gmato,nn,outfile);
	  }
	  if (ksdft) {
	    fprintf(outfile,
		    "\nxcmat for irrep %s",s->irrep_label);
	    print_array(s->xcmat,nn,outfile);
	  }
	}
      }
    }
       
       
    for (m=0; m < num_ir ; m++) {
      s = &scf_info[m];
      if (nn=s->num_so) {

	/*  form fock matrix = h+g */
	add_arr(s->hmat,s->gmat,s->fock_pac,ioff[nn]);
	    
	/* for open shell, form fock_open = h+g-q */
	if(iopen) {
	  for (i=0; i < ioff[nn] ; i++)
	    s->fock_open[i]=s->fock_pac[i]-s->gmato[i];
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
      for (m=0; m < num_ir ; m++) {
	s = &scf_info[m];
	if (nn=s->num_so) {
	  add_arr(s->fock_pac,s->xcmat,s->fock_pac,ioff[nn]);
	  if(print & 4) {
	    fprintf(outfile,"\n J+X+C gmat for irrep %d",s->irrep_label);
	    print_array(s->fock_pac,nn,outfile);
	  }
	}
      }
    }
    /* create new fock matrix in fock_pac or fock_eff */
    if(!diisflg) diis(scr,fock_c,fock_ct,c1,c2,cimax,newci);

    if(iopen) {
      for (m=0; m < num_ir ; m++) {
	s = &scf_info[m];
	if (nn=s->num_so) {
	  num_mo = s->num_mo;

	  /* transform fock_pac to mo basis */
	  tri_to_sq(s->fock_pac,fock_ct,nn);
	  /*               mxmb(s->cmat,nn,1,fock_ct,1,nn,scr,1,nn,nn,nn,nn);
			   mxmb(scr,1,nn,s->cmat,1,nn,fock_c,1,nn,nn,nn,nn);*/
	  mmult(s->cmat,1,fock_ct,0,scr,0,num_mo,nn,nn,0);
	  mmult(scr,0,s->cmat,0,fock_c,0,num_mo,nn,num_mo,0);

	  /* transform fock_open to mo basis */
	  tri_to_sq(s->fock_open,fock_ct,nn);
	  /*               mxmb(s->cmat,nn,1,fock_ct,1,nn,scr,1,nn,nn,nn,nn);
			   mxmb(scr,1,nn,s->cmat,1,nn,fock_o,1,nn,nn,nn,nn);*/
	  mmult(s->cmat,1,fock_ct,0,scr,0,num_mo,nn,nn,0);
	  mmult(scr,0,s->cmat,0,fock_o,0,num_mo,nn,num_mo,0);

	  /* form effective fock matrix in mo basis */

	  ij=0;
	  occ0 = s->occ_num[0];
	  for (i=0; i < num_mo; i++ ) {
	    for (j=0; j <= i; j++) {
	      occi = s->occ_num[i];
	      occj = s->occ_num[j];

              /* default: Guest & Saunders general form */
	      if(iter < itmax-1 && !converged && !fock_typ) {
		if(occi == occj) 
		  s->fock_eff[ij] = fock_c[i][j];
		else if(occi)
		  s->fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
		else {
		  if(occj==2.0)
		    s->fock_eff[ij] = fock_c[i][j];
		  else
		    s->fock_eff[ij] = fock_o[i][j];
		}
	      }

	      /* Guest & Saunders' form for high spin */
	      else if(iter < itmax-1 && !converged && fock_typ == 1) {
		if (occi == occj || occi)
		  s->fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
		else if (occj == 2.0) s->fock_eff[ij] = fock_c[i][j];
		else s->fock_eff[ij] = fock_o[i][j];
	      }

	      /* test form (fo fo fo) */
	      else if(iter < itmax-1 && !converged && fock_typ == 2) {
		if (occi == occj) s->fock_eff[ij] = fock_o[i][j];
		else if(occi)
		  s->fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
		else if(occj == 2.0)
		  s->fock_eff[ij] = fock_c[i][j];
		else s->fock_eff[ij] = fock_o[i][j];
	      }

	      /* test form a*(fc fc fc) */
	      else if(iter < itmax-1 && !converged && fock_typ == 3) {
		if (occi == occj) s->fock_eff[ij] = dampd*fock_c[i][j];
		else if(occi)
		  s->fock_eff[ij] = 
		    dampo*(2.0*fock_c[i][j]-fock_o[i][j]);
		else if(occj == 2.0)
		  s->fock_eff[ij] = dampo*fock_c[i][j];
		else s->fock_eff[ij] = dampo*fock_o[i][j];
	      }

	      /* test form a*(fo fo fo) */
	      else if(iter < itmax-1 && !converged && fock_typ == 4) {
		if (occi == occj) s->fock_eff[ij] = dampd*fock_o[i][j];
		else if(occi)
		  s->fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
		else if(occj == 2.0)
		  s->fock_eff[ij] = fock_c[i][j];
		else s->fock_eff[ij] = fock_o[i][j];
	      }

	      /* test form a*(2fc-fo 2fc-fo 2fc-fo) */
	      else if(iter < itmax-1 && !converged && fock_typ == 5) {
		if (occi == occj)
		  s->fock_eff[ij] = 
		    dampd*(2.0*fock_c[i][j]-fock_o[i][j]);
		else if(occi)
		  s->fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
		else if(occj == 2.0)
		  s->fock_eff[ij] = fock_c[i][j];
		else s->fock_eff[ij] = fock_o[i][j];
	      }

	      /* form for converged wavefunction */
	      else {
		if (occi == 2.0) s->fock_eff[ij]=fock_c[i][j];
		else if(occj != 2.0)
		  s->fock_eff[ij]=fock_o[i][j];
		else {
		  if(occi)
		    s->fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
		  else s->fock_eff[ij]=fock_c[i][j];
		}
	      }

	      if(j==i) {
#if 1
		if (occi == occ0 && occi)
		  s->fock_eff[ij] -= lshift;
		else 
		  if (occi) s->fock_eff[ij] -= 0.5*lshift;
#else
		if (occi == 1.0)
		  s->fock_eff[ij] += 0.5*lshift;
		else if (occi == 0.0)
		  s->fock_eff[ij] += 0.5*(lshift+0.1);
#endif
	      }
	      ij++;
	    }
	  }
	}
      }
    }
               
    for (m=0; m < num_ir ; m++) {
      s = &scf_info[m];
      if (nn=s->num_so) {
	num_mo = s->num_mo;
	occ0 = s->occ_num[0];
	if(iopen) {
	  rsp(num_mo,num_mo,ioff[num_mo],s->fock_eff,s->fock_evals,1,ctrans,tol);
	       
	  /*	       mxmb(s->cmat,1,nn,ctrans,1,nn,scr,1,nn,nn,nn,nn);*/
	  mmult(s->cmat,0,ctrans,0,scr,0,nn,num_mo,num_mo,0);
	       
	  for (i=0; i < num_mo; i++) {
	    occi = s->occ_num[i];
#if 1
	    if (occi == occ0 && occi) s->fock_evals[i] += lshift;
	    else if (occi) s->fock_evals[i] += 0.5*lshift;
#else
	    if (occi == 1.0) s->fock_evals[i] -= 0.5*lshift;
	    else if (occi == 0.0) s->fock_evals[i] -= 0.5*(lshift+0.1);
#endif
	  }

	  for(i=0; i < nn; i++)
	    for (j=0; j < num_mo; j++)
	      s->cmat[i][j] = scr[i][j];

	  mmult(s->ucmat,0,ctrans,0,scr,0,num_mo,num_mo,num_mo,0);
	  for(i=0; i < num_mo; i++)
	    for (j=0; j < num_mo; j++)
	      s->ucmat[i][j] = scr[i][j];

	}
	else {

	  /* transform fock_pac to mo basis */
	  tri_to_sq(s->fock_pac,fock_ct,nn);
	  /*               mxmb(s->cmat,nn,1,fock_ct,1,nn,scr,1,nn,nn,nn,nn);
			   mxmb(scr,1,nn,s->cmat,1,nn,fock_c,1,nn,nn,nn,nn);*/
	  mmult(s->cmat,1,fock_ct,0,scr,0,num_mo,nn,nn,0);
	  mmult(scr,0,s->cmat,0,fock_c,0,num_mo,nn,num_mo,0);

	  /*  diagonalize fock_c to get ctrans */
	  sq_rsp(num_mo,num_mo,fock_c,s->fock_evals,1,ctrans,tol);

	  if(print & 4) {
	    fprintf(outfile,"\n eigenvector for irrep %s\n",
		    s->irrep_label);
	    eivout(ctrans,s->fock_evals,num_mo,num_mo,outfile);
	  }

	  /*               mxmb(s->cmat,1,nn,ctrans,1,nn,scr,1,nn,nn,nn,nn);*/
	  mmult(s->cmat,0,ctrans,0,scr,0,nn,num_mo,num_mo,0);

	  if(print & 4) {
	    fprintf(outfile,"\n eigenvector after irrep %s\n",
		    s->irrep_label);
	    print_mat(scr,nn,num_mo,outfile);
	  }

	  for (i=0; i < nn; i++)
	    for (j=0; j < num_mo; j++)
	      s->cmat[i][j] = scr[i][j];

	  mmult(s->ucmat,0,ctrans,0,scr,0,num_mo,num_mo,num_mo,0);
	  for(i=0; i < num_mo; i++)
	    for (j=0; j < num_mo; j++)
	      s->ucmat[i][j] = scr[i][j];

	}
      }
    }

    if(converged) {
      free_block(scr);
      free_block(fock_c);
      free_block(fock_ct);
      free_block(ctrans);
      if(iopen) {
	free_block(fock_o);
	free(c1);
	free(c2);
      }
      //cleanup();
      return;
    }

    schmit(1);

    if(print & 4) {
      for(i=0; i < num_ir ; i++) {
	s = &scf_info[i];
	if (nn=s->num_so) {
	  num_mo = s->num_mo;
	  fprintf(outfile,"\northogonalized mos irrep %s\n",
		  s->irrep_label);
	  print_mat(s->cmat,nn,num_mo,outfile);
	}
      }
    }

    /* reset occupations if needed */
    if (reset_occ) {
      sortev();
      occ_calc();
    }

    /* form new density matrix */
    dmat();

    /* and form new fock matrix */

    if(iter < itmax) {
      if (direct_scf)
	formg_direct();
      else {
	if(iopen) formg_open();
	else formg_closed();
      }
    }
  }

  /* Clean up */
  for(i=0; i < num_ir ; i++) {
    s = &scf_info[i];
    if (nn=s->num_so)
      free_block(s->ucmat);
  }

}

}} // namespace psi::cscf
