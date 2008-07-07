/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.12  2004/05/03 04:32:40  crawdad
 * Major mods based on merge with stable psi-3-2-1 release.  Note that this
 * version has not been fully tested and some scf-optn test cases do not run
 * correctly beccause of changes in mid-March 2004 to optking.
 * -TDC
 *
/* Revision 1.11.4.1  2004/04/06 21:29:05  crawdad
/* Corrections to the RHF/ROHF DIIS algorithm, which was simply incorrect.
/* The backtransformation of the DIIS error vectors to the AO basis was not
/* mathematically right.
/* -TDC and EFV
/*
/* Revision 1.11  2003/08/17 22:57:37  crawdad
/* Removing libfile30 from the repository.  I believe that all code reference
/* to the library have also been properly removed.  The current version
/* passes all test cases on my systems.
/* -TDC
/*
/* Revision 1.10  2003/04/10 20:36:01  crawdad
/* Modifications to cscf to account for *very* large cases.  Mainly converted
/* terms to unsigned ints and more carefully computed pk-block sizes to avoid
/* overflows.
/* -TDC
/*
/* Revision 1.9  2002/11/24 22:52:17  crawdad
/* Merging the gbye-file30 branch into the main trunk.
/* -TDC
/*
/* Revision 1.8.2.1  2002/07/29 23:08:30  evaleev
/* A major set of changes designed to convert all psi modules to use libchkpt.
/*
/* Revision 1.8  2002/05/15 02:29:14  sherrill
/* Read from checkpoint
/*
/* Revision 1.7  2002/04/03 02:06:01  janssen
/* Finish changes to use new include paths for libraries.
/*
/* Revision 1.6  2002/03/24 18:31:19  crawdad
/* NOW it works.
/* -TDC
/*
/* Revision 1.5  2002/03/24 18:30:08  crawdad
/* Beginning mods for libpsio-based file30.  Current version works.
/* -TDC
/*
/* Revision 1.4  2002/03/24 17:28:14  crawdad
/* Minor modifications in preparation for conversion to libpsio-based file30.
/* -TDC
/*
/* Revision 1.3  2000/10/13 19:51:21  evaleev
/* Cleaned up a lot of stuff in order to get CSCF working with the new "Mo-projection-capable" INPUT.
/*
/* Revision 1.2  2000/06/22 22:15:01  evaleev
/* Modifications for KS DFT. Reading in XC Fock matrices and XC energy in formg_direct need to be uncommented (at present those are not produced by CINTS yet).
/*
/* Revision 1.1.1.1  2000/02/04 22:52:31  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.6  1999/11/11 16:00:38  localpsi
/* Fixed minor bug in occupations.  STB (11/11/99)
/*
/* Revision 1.5  1999/11/04 19:24:29  localpsi
/* STB (11/4/99) - Added the orb_mix feature which is equivalent to guess = mix
/* in G94 and also fixed restarting so that if you have different wavefuntions,
/* everything works.  Also if you specify no DOCC and SOCC and restart, if the
/* wavefunctions are different, it will guess again.
/*
/* Revision 1.4  1999/11/02 18:10:13  evaleev
/* Direct SCF improved
/*
/* Revision 1.3  1999/08/17 19:04:15  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.2  1999/08/11 19:24:53  evaleev
/* Unhardwired the size of the ioff array (set it to 1024 for now) and increased MAX_BASIS to 1024.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:27  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id: init_scf.cc 3956 2008-06-09 12:24:49Z rking $";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <libchkpt/chkpt.h>

namespace psi { namespace cscf {

void init_scf()
{
   int i,jj;
   int nn,isadr;
   int nkind,junk;
   int degen[20],*num_so;
   char char_dum[80];
   char **irr_labs;

   ioff[0] = 0;
   for (i = 1; i < MAXIOFF ; i++) {
      ioff[i] = ioff[i-1] + i;
      }

/* EFV 10/24/98 All requests for file30 should be handled with libfile30
   but for now I'll use wreadw */
   num_ir = chkpt_rd_nirreps();
   num_so = chkpt_rd_sopi();
   repnuc = chkpt_rd_enuc();
   irr_labs = chkpt_rd_irr_labs();

/* now initialize scf_info */
   
   n_so_typs=0;
   nsfmax=0;
   nbasis=0;
   
   scf_info = (struct symm *) malloc(sizeof(struct symm)*num_ir);

   /* compute nsfmax */
   for(i=0; i < num_ir; i++) { nn = num_so[i]; if(nn > nsfmax) nsfmax = nn; }

   jj=0;
   for(i=0; i < num_ir ; i++) {
      scf_info[i].num_so = nn = num_so[i];
/* EFV 10/24/98 degeneracy = 1
      scf_info[i].degeneracy = degen[i]; */
      scf_info[i].nclosed = 0;
      scf_info[i].nopen = 0;
      scf_info[i].nhalf = 0;
      scf_info[i].os_num = 0;
      scf_info[i].ideg = 0;

      scf_info[i].irrep_label = irr_labs[i];
/*      scf_info[i].irrep_label[4] = '\0';*/
      jj += 4;

      nbasis += nn;
      if (nn) {
         n_so_typs++;
	 /*         if (nn > nsfmax) nsfmax = nn; */

         scf_info[i].smat = (double *) init_array(ioff[nn]);
         scf_info[i].tmat = (double *) init_array(ioff[nn]);
         scf_info[i].hmat = (double *) init_array(ioff[nn]);
         scf_info[i].pmat = (double *) init_array(ioff[nn]);
         scf_info[i].pmato = (double *) NULL;
         scf_info[i].pmat2 = (double *) NULL;
         scf_info[i].pmato2 = (double *) NULL;
         scf_info[i].dpmat = (double *) init_array(ioff[nn]);
         scf_info[i].dpmato = (double *) NULL;
         scf_info[i].fock_pac = (double *) init_array(ioff[nn]);
         scf_info[i].fock_eff = (double *) NULL;
         scf_info[i].fock_open = (double *) NULL;
         scf_info[i].gmat = (double *) init_array(ioff[nn]);
         scf_info[i].gmato = (double *) NULL;
         scf_info[i].occ_num = (double *) init_array(nn);
         scf_info[i].fock_evals = (double *) init_array(nn);
         scf_info[i].cmat = block_matrix(nn,nn);
	 /* TDC(6/19/96) - Added array for saving original MO vector */
	 scf_info[i].cmat_orig = block_matrix(nn,nn);
         scf_info[i].sahalf = block_matrix(nn,nn);
	 scf_info[i].pinv = block_matrix(nn,nn);
         /* STB(4/1/98) - Added array to store the eigenvalues of the
	                  core hamiltonian for mo guessing*/
         scf_info[i].hevals = (double *) init_array(nn);
	 /* Need separate XC Fock for KS DFT */
	 if (ksdft)
	   scf_info[i].xcmat = init_array(ioff[nn]);
         }
     }
   free(irr_labs);
   free(num_so);
   irr_labs = NULL;
   num_so = NULL;
   /* read in number of atoms and nuclear charges and total number of MO*/
   natom = chkpt_rd_natom();
   zvals = chkpt_rd_zvals();
   nbfso = chkpt_rd_nso();
   
/* Initialize arrays to hold energy and symmetry arrays */
   ener_tot = (double *) init_array(nbfso);
   symm_tot = (int *) init_int_array(nbfso);
   
} 

void init_scf2()
{
   int i,j,k;
   int n,nn,m,mm;
   int junk;
   unsigned int opconst,outbuf,mxcoef3,ntri,mtri;
   struct symm *s;

   opconst = (iopen) ? 3 : 2;

   mtri = ioff[scf_info[0].num_so];
   mxcoef3 = opconst*(mtri*(mtri+1)/2);

   scf_info[0].ideg = 0;

   so2symblk = init_int_array(nbasis);

   junk=j=0;
   for (i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      nn = s->num_so;
      if (i) {
         m = ioff[nn];
         mm = m*(m+1)/2;
         mxcoef3 += opconst*ioff[nn]*mtri;
         mxcoef3 += opconst*mm;
         mtri += m;
         if (nn <= 0) s->ideg = scf_info[i-1].ideg;
         else {
            do {
               n=scf_info[j].num_so;
               j++;
               } while(!n); 
            s->ideg = scf_info[i-1].ideg+n;
            }
         }
       for(k=s->ideg;k<s->ideg+nn;k++)
         so2symblk[k] = i;
      if(s->nopen || s->nhalf) s->os_num = junk++;
      }
   ntri = nbasis*(nbasis+1)/2;

   if (direct_scf) { /* No I/O will be done */
     fprintf(outfile,
	     "\n  direct formation of the Fock matrix is requested\n");
   }
   else {            /* Figure out the size of the PK-buffers */
     readflg = 0;
     maxbuf = mxcoef3/opconst+5;
     outbuf = 884736;
     outbuf = MIN0(outbuf,mxcoef3);
     if(outbuf == 884736) {
       unsigned int pass = mxcoef3/outbuf+1;
       readflg = 1;
       maxbuf = (mxcoef3/pass)/opconst + 2;
       if(iopen) maxbuf /= 2;
#if defined(AIXV3)||defined(SGI)
       maxbuf= (iopen) ? 8192*3 : 8192*5;
       pass = mxcoef3/(maxbuf*opconst)+1;
#endif
       fprintf(outfile,
	       "\n  using buffered io, %u buffers, each %u bytes in size\n",
	       pass,maxbuf*opconst*sizeof(double));
       if(print) fprintf(outfile,"  mxcoef3 = %u maxbuf = %u\n",mxcoef3,maxbuf);
       if(print) fprintf(outfile,"  outbuf = %u\n",outbuf);
     }
     else fprintf(outfile,"\n  keeping integrals in %u bytes of core\n",
		  maxbuf*opconst*sizeof(double));
   }
     
   fflush(outfile);
   return;
   
}

}} // namespace psi::cscf
