/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.3  2002/12/06 15:50:32  crawdad
 * Changed all exit values to PSI_RETURN_SUCCESS or PSI_RETURN_FAILURE as
 * necessary.  This is new for the PSI3 execution driver.
 * -TDC
 *
/* Revision 1.2  2000/10/13 19:51:22  evaleev
/* Cleaned up a lot of stuff in order to get CSCF working with the new "Mo-projection-capable" INPUT.
/*
/* Revision 1.1.1.1  2000/02/04 22:52:34  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.1  1999/11/02 23:56:00  localpsi
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
/* Revision 1.1.1.1  1999/04/12 16:59:28  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
/*
 * Revision 1.1  1991/06/15  20:22:42  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: schmit_uhf.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void schmit_uhf(int all)
{
    int i,j,ij,nn,num_mo;
    int n,m,t,ncol;
    double *v,**ctmp,vtmp;
    struct symm *s;
    struct spin *sp;
    
    v = (double *) init_array(nsfmax);
    ctmp = (double **) init_matrix(nsfmax,nsfmax);
    
    for(t=0;t < 2; t++){
	sp = &spin_info[t];
	for(n=0; n < num_ir ; n++) {
	    s= &scf_info[n];
	    if(nn=s->num_so) {
		num_mo = s->num_mo;
		for (i=0; i < nn ; i++)
		    for (j=0; j < num_mo ; j++)
			ctmp[j][i] = sp->scf_spin[n].cmat[i][j];
		ncol = sp->scf_spin[n].noccup;
		if(s->nhalf) ncol++;
		if(all) ncol = num_mo;
		if(!ncol) continue;
		for(m=0; m < ncol ; m++) {
		    v[0]=ctmp[m][0]*s->smat[0];
		    for(i=1; i < nn ; i++) {
			for(j=0,vtmp=0.0; j < i ; j++) {
			    ij=ioff[i]+j;
			    vtmp += ctmp[m][j]*s->smat[ij];
			    v[j] += ctmp[m][i]*s->smat[ij];
			}
			v[i] = vtmp+ctmp[m][i]*s->smat[ioff[i]+j];
		    }
		    for(i=0,vtmp=0.0; i < nn ; i++) vtmp += v[i]*ctmp[m][i];
		    if(!vtmp) {
			exit(PSI_RETURN_FAILURE);
		    }
		    if(vtmp < 10.0e-20) vtmp = 10.0e-20;
		    vtmp = 1.0/sqrt(vtmp);
		    
		    for(i=0; i < nn ; i++) {
			v[i] *= vtmp;
			ctmp[m][i] *= vtmp;
		    }
		    
		    if(m < ncol-1) {
			for(i=m+1,vtmp=0.0; i < ncol ; i++) {
			    for(j=0,vtmp=0.0; j<nn ;j++) 
				vtmp += v[j]*ctmp[i][j];
			    for(j=0; j < nn ; j++) 
				ctmp[i][j] -= vtmp*ctmp[m][j];
			}
		    }
		}
		
		for (i=0; i < nn ; i++)
		    for (j=0; j < num_mo ; j++)
			sp->scf_spin[n].cmat[i][j] = ctmp[j][i];
		
	    }
	}
    }

    free(v);
    free_matrix(ctmp,nsfmax);
}

}} // namespace psi::cscf
