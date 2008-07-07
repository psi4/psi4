/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.2  2004/05/03 04:32:40  crawdad
 * Major mods based on merge with stable psi-3-2-1 release.  Note that this
 * version has not been fully tested and some scf-optn test cases do not run
 * correctly beccause of changes in mid-March 2004 to optking.
 * -TDC
 *
/* Revision 1.1.1.1.12.2  2004/04/21 15:45:07  evaleev
/* Modified DIIS algorithm for RHF and ROHF to work in OSO basis rather than in
/* AO basis, to avoid difficulties of transforming between MO and AO bases
/* when linear dependencies are present.
/*
/* Revision 1.1.1.1.12.1  2004/04/06 21:29:05  crawdad
/* Corrections to the RHF/ROHF DIIS algorithm, which was simply incorrect.
/* The backtransformation of the DIIS error vectors to the AO basis was not
/* mathematically right.
/* -TDC and EFV
/*
/* Revision 1.1.1.1  2000/02/04 22:52:29  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.3  1999/11/02 23:55:55  localpsi
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
/* Revision 1.2  1999/08/17 19:04:13  evaleev
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
 * */

static char *rcsid = "$Id: diis.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

double *btemp, **bold, **bmat;
struct diis_mats {
  double ***fock_c;       /* Closed-shell Fock matrix in AO basis */
  double ***fock_o;       /* Open-shell Fock matrix in AO basis */
  double ***error;        /* Error matrix in AO basis */
  } *diism,dtemp;

extern double delta;

void diis_free(void) {
  int i,j;
  free(btemp); btemp = NULL;
  free_block(bold); bold = NULL;
  free_block(bmat);  bmat = NULL;
  for (i=0; i<ndiis; ++i) {
    for(j=0; j<num_ir; ++j) {
      if(scf_info[j].num_so) {
        free_block(diism[i].fock_c[j]);
        diism[i].fock_c[j] = NULL;
        free_block(diism[i].error[j]);
        diism[i].error[j] = NULL;
        if (iopen) {
          free_block(diism[i].fock_o[j]);
          diism[i].fock_o[j] = NULL;
        }
      }
    }
    free(diism[i].fock_c); diism[i].fock_c = NULL;
    if (iopen) { free(diism[i].fock_o); diism[i].fock_o = NULL; }
    free(diism[i].error);  diism[i].error = NULL;
  }
  free(diism); diism = NULL;
}

void diis(double** scr1, double** scr2, double** scr3, double* c1, double* c2, double cim, int newci)
{
  int i,j,k,ij;
  int errcod;
  int m,nn,num_mo;
  double occi, occj;
  int _try = 0;
  int last = iter-1;
  int col = iter+1;
  double etemp, dotp, norm, determ, etempo;
  double scale;
  struct symm *s;
  struct diis_mats *d;
  int diis_print=0;
  double **C;  /* AO->MO transform = P^-1 * C */
  double **X, **S; /* scratch matrix */
  double **tmp1;

  if(diism == NULL) {
    bmat = (double **) block_matrix(ndiis+1,ndiis+1);
    bold = (double **) block_matrix(ndiis,ndiis);
    btemp = (double *) init_array(ndiis+1);

    diism = (struct diis_mats *) malloc(sizeof(struct diis_mats)*ndiis);

    for(i=0; i < ndiis ; i++) {
      diism[i].fock_c = (double ***) malloc(sizeof(double **)*num_ir);
      if(iopen) 
	diism[i].fock_o = (double ***) malloc(sizeof(double **)*num_ir);
      diism[i].error = (double ***) malloc(sizeof(double **)*num_ir);
      for(j=0; j < num_ir ; j++) {
	if(nn=scf_info[j].num_so) {
	  num_mo = scf_info[j].num_mo;
	  diism[i].fock_c[j] = (double **) block_matrix(nn,nn);
	  if(iopen) diism[i].fock_o[j] = (double **) block_matrix(nn,nn);
	  diism[i].error[j] = (double **) block_matrix(num_mo,num_mo);
	}
      }
    }
  }

  scale = 1.0 + dampsv;

  if (iter > ndiis) {
    last = ndiis-1;
    col = ndiis+1;
    dtemp = diism[0];
    for (i=0; i < last ; i++) {
      diism[i] = diism[i+1];
    }
    diism[last] = dtemp;
  }
      
  /* Debugging stuff */
   
  errcod = ip_boolean("DIIS_PRINT",&diis_print,0);
  if(iter == 1 && diis_print)
    ffile(&diis_out,"diis_out.dat",0);
   
  /* save ao fock matrices in fock_save */
   
  d = &diism[last];
  for (m=0; m < num_ir ; m++) {
    s = &scf_info[m];
    if(nn=s->num_so) {
      num_mo = s->num_mo;
      tri_to_sq(s->fock_pac,d->fock_c[m],nn);
      if(iopen) tri_to_sq(s->fock_open,d->fock_o[m],nn);

      /* form error matrix in mo basis */
      mmult(s->cmat,1,d->fock_c[m],0,scr1,0,num_mo,nn,nn,0);
      mmult(scr1,0,s->cmat,0,scr2,0,num_mo,nn,num_mo,0);
      if(iopen) {
	mmult(s->cmat,1,d->fock_o[m],0,scr1,0,num_mo,nn,nn,0);
	mmult(scr1,0,s->cmat,0,scr3,0,num_mo,nn,num_mo,0);
      }

      for (i=0; i < num_mo; i++) {
	occi = s->occ_num[i];
	for (j=0; j <= i ; j++ ) {
	  occj = s->occ_num[j];
	  if (!iopen) {
	    if (occi == 0.0 && occj != 0.0 ) {
	      scr1[i][j]= scr2[i][j];
	      scr1[j][i]= scr2[i][j];
	    }
	    else {
	      scr1[i][j]=scr1[j][i]=0.0;
	    }
	  }
	  else if(!twocon) {
	    if ((occi == 1.0 || occi == 0.5) && occj == 2.0 )
	      etemp = scr2[i][j]-0.5*scr3[i][j];
	    else if(occi == 0.0) {
	      if (occj == 2.0) etemp = scr2[i][j];
	      else if (occj != 0.0) etemp = 0.5*scr3[i][j];
	      else etemp = 0.0;
	    }
	    else etemp = 0.0;
	    scr1[i][j] = scr1[j][i] = etemp;
	  }
	  else {
	    if (occi != 2.0 && occi != 0.0 && occj == 2.0 )
	      etemp = cim*(c1[m]*scr2[i][j]-c2[m]*scr3[i][j]);
	    else if(occi == 0.0) {
	      if (occj == 2.0) etemp = scr2[i][j];
	      else if (occj != 0.0) etemp = cim*scr3[i][j];
	      else etemp = 0.0;
	    }
	    else etemp = 0.0;
	    scr1[i][j] = scr1[j][i] = etemp;
	  }
	}
      }

      /* transform the error matrix into the OSO basis */
      X = block_matrix(num_mo, num_mo);
      C_DGEMM('n', 'n', num_mo, num_mo, num_mo, 1.0, s->ucmat[0], num_mo, scr1[0], nsfmax, 0.0, X[0], num_mo);
      C_DGEMM('n', 't', num_mo, num_mo, num_mo, 1.0, X[0], num_mo, s->ucmat[0], num_mo, 0.0, d->error[m][0], num_mo);
      free_block(X);
      X = NULL;

      for(i=0; i < num_mo ; i++) {
	for(j=0; j <= i ; j++) {
	  etemp=fabs(scr1[i][j]);
	  diiser = MAX0(diiser,etemp);
	}
      }
    }
  }
               
  /* then set up b matrix */

  if (iter > ndiis) {
    for (i=0; i < last ; i++) {
      for (j=0; j <= i ; j++) {
	bold[i][j]=bold[j][i]=bold[i+1][j+1];
      }
    }
  }
  for (i=0; i <= last ; i++) {
    etemp=0.0;
    for (m=0; m < num_ir ; m++) {
      s = &scf_info[m];
      if(nn=s->num_so) {
	num_mo = s->num_mo;
	sdot(diism[i].error[m],diism[last].error[m],num_mo,&dotp);
	etemp += dotp;
      }
    }
    bold[i][last]=bold[last][i] = etemp;
  }

  bmat[0][0] = 0.0;
  btemp[0] = -1.0;
  norm = 1.0/bold[0][0];
  for (i=1; i <= last+1 ; i++) {
    bmat[i][0]=bmat[0][i] = -1.0;
    btemp[i] = 0.0;
    for (j=1; j <= i ; j++) {
      bmat[i][j]=bmat[j][i] = bold[i-1][j-1]*norm;
      if(i==j) bmat[i][j] *= scale;
    }
  }
   
  if(diis_print){
    fprintf(diis_out,"\nBMAT for iter %d",iter);
    print_mat(bmat,col,col,diis_out);
  }
	
  if (iter-1) {
    flin(bmat,btemp,col,1,&determ);

    /* test for poorly conditioned equations */
    while (fabs(determ) < 1.0e-19 && _try < last) {

      _try++;
      col--;

      bmat[0][0] = 0.0;
      btemp[0] = -1.0;
      norm=1.0/bold[_try][_try];
      for (i=1; i <= ndiis-_try ; i++) {
	bmat[i][0]=bmat[0][i] = -1.0;
	for (j=1; j <= i ; j++) {
	  bmat[i][j]=bmat[j][i]=bold[i+_try-1][j+_try-1]*norm;
	  if(i==j) bmat[i][j] *= scale;
	}
	btemp[i] = 0.0;
      }
	 
      if(diis_print){
	fprintf(diis_out,"\nCorrected BMAT for iter %d",iter);
	print_mat(bmat,col,col,diis_out);
      }
         
      flin(bmat,btemp,col,1,&determ);
    }
      
    if(fabs(determ) < 10.0e-20) {
      printf(" try %d no good\n",_try);
      return;
    }
      
    if(iter >= it_diis) {
      for (m=0; m < num_ir ; m++) {
	s = &scf_info[m];
	if(nn=s->num_so) {
	  for (i=ij=0; i < nn ; i++) {
	    for (j=0; j <= i ; j++,ij++) {
	      int kk=1;
	      etemp=0.0;
	      etempo=0.0;
	      for (k=_try; k < last+1 ; k++) {
		if(iopen) etempo += btemp[kk]*diism[k].fock_o[m][i][j];
		etemp += btemp[kk]*diism[k].fock_c[m][i][j];
		kk++;
	      }
	      if(iopen) s->fock_open[ij] = etempo;
	      s->fock_pac[ij] = etemp;
	    }
	  }
	}
      }
    }
  }
}

}} // namespace psi::cscf
