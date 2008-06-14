/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.3  2003/08/09 17:39:56  crawdad
 * I added the ability to determine frozen core orbitals for UHF references to
 * cleanup.c.  I also commented out ip_cwk_clear and ip_cwk_add calls in
 * cleanup.c, guess.c, scf_input.c and scf_iter_2.c.  These calls were (1) poor
 * design and (2) interfering with default ip_tree behavior needed to simplify
 * the format of input.dat.
 * -TDC
 *
/* Revision 1.2  2002/03/25 02:17:36  janssen
/* Get rid of tmpl.  Use new naming scheme for libipv1 includes.
/*
/* Revision 1.1.1.1  2000/02/04 22:52:32  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.2  1999/08/17 19:04:17  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:28  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id: scf_iter_2.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <libipv1/ip_lib.h>

namespace psi { namespace cscf {

/* scf procedure for two-configuration scf */

void scf_iter_2()

{
   int i,j,k,m,ij;
   int joff;
   int nn,num_mo;
   int iju=0;
   int newci=1;
   int opblk=0;
   int *optest;
   int errcod;
   double ciconv = 10.0e-9;  // FAE it was 10.0e-9
   double c1i,c1ii;
   double eci1,eci2,eci3;
   double ci1,ci2,ci1old,term,hoff,amax,cn,sn;
   double cimax,cimin,incr;
   double occi, occj, occ0, den;
   double **scr;
   double *c1, *c2;
   double **fock_c, **fock_o;
   double **fock_ct,**fock_eff;
   double **ctrans;
   double tol = 1.0e-16;
   struct symm *s;

   diiser = 0.0;
   /* optest = (int *) init_array(ioff[nbasis]/2); */
   optest = (int *) init_int_array(ioff[nbasis]);
   scr = (double **) init_matrix(nsfmax,nsfmax);
   fock_c = (double **) init_matrix(nsfmax,nsfmax);
   fock_ct = (double **) init_matrix(nsfmax,nsfmax);
   ctrans = (double **) init_matrix(nsfmax,nsfmax);
   c1 = (double *) init_array(num_ir);
   c2 = (double *) init_array(num_ir);
   fock_o = (double **) init_matrix(nsfmax,nsfmax);

   for (i=0; i < num_ir ; i++) {
		 s = &scf_info[i];
		 if (nn=s->num_so) {
		   num_mo = s->num_mo;
		   s->ucmat = block_matrix(num_mo, num_mo);
		   for(j=0; j<num_mo; j++)
		     s->ucmat[j][j] = 1.0;
     }
   }

	 
   /*
   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":SCF");
   */

   incr=0.25;
   errcod = ip_data("INCR","%lf",&incr,0);
   fprintf(outfile,"\n  tcscf increment = %f\n",incr);


/* find open shells */

   for (i=0; i < num_ir ; i++) {
      iju += scf_info[i].num_mo;
      if(scf_info[i].nopen) {
         iju = ioff[iju];
         break;
         }
      }
           
   for (i=0; i < num_ir ; i++) {
      s= &scf_info[i];
      if(nn=s->num_so) {
	 num_mo = s->num_mo;
         for (j=0; j < num_mo ; j++) {
            if (s->occ_num[j] == 1.0) {
               if (!opblk) {
                  opblk++;
                  opshl1 = j;
                  opblk1 = i;
                  }
               else {
                  opshl2 = j;
                  opblk2 = i;
                  }
               }
            }
         }
      }

/* set up array of flags indicating open shells */

   for (k=joff=0; k < num_ir ; k++) {
      s = &scf_info[k];
      if (nn=s->num_so) {
         for (i=0; i < nn ; i++)
            for (j=0; j <= i ; j++)
               if(s->nopen) optest[ioff[i+joff]+j+joff] = 1;
         joff += nn;
         }
      }

/* set up c1 and c2  */
   for (i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      if (nn=s->num_so) {
	 num_mo = s->num_mo;
         for (j=0; j < num_mo ; j++) {
            if(s->occ_num[j]==2.0) c1[i]=2.0;
            if(s->occ_num[j]==1.0) c2[i]=1.0;
            }
         if(i == opblk1) c1i=c1[i];
         if(i == opblk2) c1ii=c1[i];
         den = c1[i]-c2[i];
         if(den != 0.0) {
            c1[i] /= den;
            c2[i] /= den;
            }
         }
      }
   ci1    =  1.0;
   ci1old =  0.0;  // FAE Give it some arbitrary value

/* and iterate */
   
   if (second_root) 
      fprintf(outfile, "Warning: Finding second root of TCSCF\n\n");

   for (iter=0; iter < itmax ; ) {
      if(iter==stop_lshift){
        lshift = 1.0;
        fprintf(outfile,"Setting level shift to %lf\n",lshift);
      }
      schmit(1);

      if(!newci || fabs(ci1old-ci1) < ciconv && iter) goto L2;

/* do some funky stuff to get coeffs ? */

      dmat_2(opblk1);
      dmat_2(opblk2);

      for (i=0; i < num_ir ; i++) {
         double d1,d2;
         s = &scf_info[i];
         if (nn=s->num_so) {
            for (j=0; j < ioff[nn] ; j++) {
               d1 = s->pmat2[j];
               d2 = s->pmato2[j];
               s->pmat2[j] = d1+d2;
               s->pmato2[j] = d1-d2;
               }
            }
         }

      formg_open();

      eci1 = eci2 = eci3 = 0.0;
      for (i=0; i < num_ir ; i++) {
         s = &scf_info[i];
         if (nn=s->num_so) {
            for (j=0; j < ioff[nn] ; j++) {
               eci1 += s->pmato2[j]*s->hmat[j];
               eci2 += s->pmato2[j]*s->dpmat[j];
               eci3 += s->pmato2[j]*s->dpmato[j];
               }
            }
         }

      term = 0.5*(eci1+0.5*eci2);
      hoff = -0.25*eci3;
      amax = term*term + hoff*hoff;
      amax = sqrt(amax);
      if(term < 0.0) amax = -(amax);
      cn = (amax+term)/(2*amax);
      cn = sqrt(cn);
      sn = hoff/(cn*(amax+amax));
      ci1old = ci1;

      if (!second_root) {
         if(term >= 0.0) {
            ci1 = cn;
            ci2 = -sn;
            }
         else {
            ci1=sn;
            ci2=cn;
            }
         }
      else {
         if(term >= 0.0) {
            ci1 = sn;
            ci2 = cn;
            }
         else {
            ci1=cn;
            ci2=-sn;
            }
      }

      save_ci1 = ci1;
      save_ci2 = ci2;

      scf_info[opblk1].occ_num[opshl1] = ci1*ci1*2.0;
      scf_info[opblk2].occ_num[opshl2] = ci2*ci2*2.0;
      den = c1i-scf_info[opblk1].occ_num[opshl1];
      c1[opblk1] = c1i/den;
      c2[opblk1] = scf_info[opblk1].occ_num[opshl1]/den;
      den = c1ii-scf_info[opblk2].occ_num[opshl2];
      c1[opblk2] = c1ii/den;
      c2[opblk2] = scf_info[opblk2].occ_num[opshl2]/den;
      alpha1 = 1.0 - 2.0/scf_info[opblk1].occ_num[opshl1];
      alpha2 = -1.0/(ci1*ci2);
      alpha3 = 1.0 - 2.0/scf_info[opblk2].occ_num[opshl2];

      cimax = MAX0(ci1,ci2);
      cimin = MIN0(ci1,ci2);
      cimax = 1.0/(cimax*cimax+dampsv);
      if(fabs(cimin) >= 0.1) cimax = 1.0; 

      fprintf(outfile,"  ci coeffs %14.10f %14.10f\n",ci1,ci2);

L2:

/* form density matrix */
      
      dmat();

/* form g matrix */

      formg_two(iju,optest);

      for (m=0; m < num_ir ; m++) {
         s = &scf_info[m];
         if (nn = s->num_so) {

          /*  form fock matrix = h+g */
            add_arr(s->hmat,s->gmat,s->fock_pac,ioff[nn]);

          /* form fock_open = h+g-q */

            for (i=0; i < ioff[nn] ; i++)
               s->fock_open[i]=s->fock_pac[i]-s->gmato[i];
            }
         }

      newci = ecalc(incr);

   /* create new fock matrix in fock_eff */
      if(!diisflg) diis(scr,fock_c,fock_ct,c1,c2,cimax,newci);

      for (m=0; m < num_ir ; m++) {
         s = &scf_info[m];
         if (nn=s->num_so) {
	    num_mo = s->num_mo;
         /* transform fock_pac to mo basis */
            tri_to_sq(s->fock_pac,fock_ct,nn);
/*            mxmb(s->cmat,nn,1,fock_ct,1,nn,scr,1,nn,nn,nn,nn);
            mxmb(scr,1,nn,s->cmat,1,nn,fock_c,1,nn,nn,nn,nn);*/
	    mmult(s->cmat,1,fock_ct,0,scr,0,num_mo,nn,nn,0);
	    mmult(scr,0,s->cmat,0,fock_c,0,num_mo,nn,num_mo,0);

         /* transform fock_open to mo basis */
            tri_to_sq(s->fock_open,fock_ct,nn);
/*            mxmb(s->cmat,nn,1,fock_ct,1,nn,scr,1,nn,nn,nn,nn);
            mxmb(scr,1,nn,s->cmat,1,nn,fock_o,1,nn,nn,nn,nn);*/
	    mmult(s->cmat,1,fock_ct,0,scr,0,num_mo,nn,nn,0);
	    mmult(scr,0,s->cmat,0,fock_o,0,num_mo,nn,num_mo,0);

         /* form effective fock matrix in mo basis */

            occ0 = s->occ_num[0];
            for (i=ij=0; i < num_mo; i++) {
               for (j=0; j <= i; j++,ij++) {
                  occi = s->occ_num[i];
                  occj = s->occ_num[j];

              /* default: Guest & Saunders general form */
                  if(iter < itmax-1 && !converged && !fock_typ) {
                     if(occi == occj)
                       s->fock_eff[ij] = fock_c[i][j];
                     else if(occi)
                       s->fock_eff[ij] = 
                         cimax*(c1[m]*fock_c[i][j]-c2[m]*fock_o[i][j]);
                     else {
                       if(occj==2.0)
                         s->fock_eff[ij] = cimax*fock_c[i][j];
                       else
                         s->fock_eff[ij] = cimax*fock_o[i][j];
                       }
                     }

               /* Guest & Saunders' form for high spin */
                  else if(iter < itmax-1 && !converged && fock_typ == 1) {
                     if (occi == occj || occi)
                        s->fock_eff[ij] = c1[m]*fock_c[i][j]-c2[m]*fock_o[i][j];
                     else if (occj == 2.0) s->fock_eff[ij] = fock_c[i][j];
                     else s->fock_eff[ij] = fock_o[i][j];
                     }

               /* test form (fo fo fo) */
                  else if(iter < itmax-1 && !converged && fock_typ == 2) {
                     if (occi == occj) s->fock_eff[ij] = fock_o[i][j];
                     else if(occi)
                        s->fock_eff[ij] = c1[m]*fock_c[i][j]-c2[m]*fock_o[i][j];
                     else if(occj == 2.0)
                        s->fock_eff[ij] = fock_c[i][j];
                     else s->fock_eff[ij] = fock_o[i][j];
                     }

            /* test form a*(fc fc fc) */
                  else if(iter < itmax-1 && !converged && fock_typ == 3) {
                     if (occi == occj) s->fock_eff[ij] = dampd*fock_c[i][j];
                     else if(occi)
                        s->fock_eff[ij] =
                          dampo*(c1[m]*fock_c[i][j]-c2[m]*fock_o[i][j]);
                     else if(occj == 2.0)
                        s->fock_eff[ij] = dampo*fock_c[i][j];
                     else s->fock_eff[ij] = dampo*fock_o[i][j];
                     }

            /* test form a*(fo fo fo) */
                  else if(iter < itmax-1 && !converged && fock_typ == 4) {
                     if (occi == occj) s->fock_eff[ij] = dampd*fock_o[i][j];
                     else if(occi)
                        s->fock_eff[ij] = 
                          dampo*(c1[m]*fock_c[i][j]-c2[m]*fock_o[i][j]);
                     else if(occj == 2.0)
                        s->fock_eff[ij] = dampo*fock_c[i][j];
                     else s->fock_eff[ij] = dampo*fock_o[i][j];
                     }

            /* test form a*(2fc-fo 2fc-fo 2fc-fo) */
                  else if(iter < itmax-1 && !converged && fock_typ == 5) {
                     if (occi == occj)
                        s->fock_eff[ij] =
                          dampd*(c1[m]*fock_c[i][j]-c2[m]*fock_o[i][j]);
                     else if(occi)
                        s->fock_eff[ij] = 
                          dampo*(c1[m]*fock_c[i][j]-c2[m]*fock_o[i][j]);
                     else if(occj == 2.0)
                        s->fock_eff[ij] = dampo*fock_c[i][j];
                     else s->fock_eff[ij] = dampo*fock_o[i][j];
                     }

            /* form for converged wavefunction */
                  else {
                     if (occi == 2.0) s->fock_eff[ij]=fock_c[i][j];
                     else if(occj != 2.0)
                       s->fock_eff[ij]=fock_o[i][j];
                     else {
                       if(occi)
                         s->fock_eff[ij] = 
                           c1[m]*fock_c[i][j]-c2[m]*fock_o[i][j];
                       else s->fock_eff[ij]=fock_c[i][j];
                       }
                     }

                  if(j==i) {
                     if (occi == occ0 && occi) s->fock_eff[ij] -= lshift;
                     else if (occi) s->fock_eff[ij] -= 0.5*lshift;
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

            rsp(num_mo,num_mo,ioff[num_mo],s->fock_eff,s->fock_evals,1,ctrans,tol);

/*            mxmb(s->cmat,1,nn,ctrans,1,nn,scr,1,nn,nn,nn,nn);*/
	    mmult(s->cmat,0,ctrans,0,scr,0,nn,num_mo,num_mo,0);

            for (i=0; i < num_mo; i++) {
               occi = s->occ_num[i];
               if (occi == occ0 && occi) s->fock_evals[i] += lshift;
               else if (occi) s->fock_evals[i] += 0.5*lshift;

               }
	    for(i=0; i < nn; i++)
	      for (j=0; j < num_mo; j++)
		s->cmat[i][j] = scr[i][j];
            }
         }

      if(converged) {
         free_matrix(scr,nsfmax);
         free_matrix(fock_c,nsfmax);
         free_matrix(fock_ct,nsfmax);
         free_matrix(ctrans,nsfmax);
         free_matrix(fock_o,nsfmax);
         free(c1);
         free(c2);
         free(optest);
         /* Clean up */
				 for(i=0; i < num_ir ; i++) {
				   s = &scf_info[i];
				   if (nn=s->num_so)
				   free_block(s->ucmat);
				 }
         //cleanup();
         }
      }
   }

}} // namespace psi::cscf
