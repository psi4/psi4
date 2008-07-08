/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.2  2001/06/29 20:39:28  evaleev
 * Modified cscf to use libpsio to store supermatrix files.
 *
/* Revision 1.1.1.1  2000/02/04 22:52:30  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.2  1999/08/17 19:04:15  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:26  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id: formg2.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void formg_two(int iju, int* optest)
{
   register int i,j,k,joff,nn;
   register int ij,kl;
   int ilast,num;
   int tmpsiz;
   double dotest,tmpval,qtemp;
   struct o_pkints *o_temp;
   struct c_pkints *c_temp;

   tmpsiz=ioff[nbasis];

   if(gtmp2 == NULL) {
      gtmp2 = (double *) init_array(tmpsiz);
      gtmpo2 = (double *) init_array(tmpsiz);
      ptmp2 = (double *) init_array(tmpsiz);
      ptmpo2 = (double *) init_array(tmpsiz);
      }
   else {
      bzero(gtmp2,sizeof(double)*tmpsiz);
      bzero(gtmpo2,sizeof(double)*tmpsiz);
      }
 
   for(k=joff=0; k < num_ir ; k++) {
      if(nn=scf_info[k].num_so) {
         for(i=0; i < nn ; i++)
            for(j=0; j <= i ; j++) {
#if NEW2C
               ptmp2[ioff[i+joff]+j+joff] = scf_info[k].dpmat[ioff[i]+j];
               ptmpo2[ioff[i+joff]+j+joff] = scf_info[k].dpmato[ioff[i]+j];
#else
               ptmp2[ioff[i+joff]+j+joff] = scf_info[k].pmat[ioff[i]+j];
               ptmpo2[ioff[i+joff]+j+joff] = scf_info[k].pmato[ioff[i]+j];
#endif
               }
         joff += nn;
         }
      }

   if(!wherec) {
      /* int_nums_o = (int *) init_array(num_bufs_o+1); */
      int_nums_o = (int *) init_int_array(num_bufs_o+1);
      /* int_nums_c = (int *) init_array(num_bufs_c+1); */
      int_nums_c = (int *) init_int_array(num_bufs_c+1);
      for(i=1; i < num_bufs_o ; i++) int_nums_o[i]=maxbuf;
      for(i=1; i < num_bufs_c ; i++) int_nums_c[i]=maxbuf;
      int_nums_o[num_bufs_o]=lasto;
      int_nums_c[num_bufs_c]=lastc;
      whereo=num_bufs_o;
      wherec=num_bufs_c;
      }

   num=int_nums_o[whereo];
   for (j=0; j < num_bufs_o ; j++) {
      o_temp = o_outbuf;

      for (i=num; i ; i--,o_temp++) {
         ij = (*o_temp).ij;
         kl = (*o_temp).kl;
         tmpval = (*o_temp).pval;
         dotest = (*o_temp).qval;

         gtmp2[ij] += ptmp2[kl]*tmpval;
         gtmp2[kl] += ptmp2[ij]*tmpval;
         if(optest[ij] && optest[kl]) {
            if (ij < iju) {
               gtmpo2[ij] += alpha1*ptmpo2[kl]*tmpval;
               gtmpo2[kl] += alpha1*ptmpo2[ij]*tmpval;
               }
            else if (kl < iju) {
               qtemp = tmpval + alpha2*dotest;
               gtmpo2[ij] += ptmpo2[kl]*qtemp;
               gtmpo2[kl] += ptmpo2[ij]*qtemp;
               }
            else {
               gtmpo2[ij] += alpha3*ptmpo2[kl]*tmpval;
               gtmpo2[kl] += alpha3*ptmpo2[ij]*tmpval;
               }
            }
         }
   
      if (readflgo && j < num_bufs_o-1) {
         if(whereo==num_bufs_o) {
	    PKmat.bufpos = PSIO_ZERO;
            whereo=0;
            }
         whereo++;
         num=int_nums_o[whereo];
	 psio_read(PKmat.unit, PKmat.key, (char *) o_outbuf, sizeof(struct o_pkints)*num,
		   PKmat.bufpos, &(PKmat.bufpos));
         }
      }

   num=int_nums_c[wherec];
   for (j=0; j < num_bufs_c ; j++) {
      c_temp = c_outbuf;

      for (i=num; i ; i--,c_temp++) {
         ij = (*c_temp).ij;
         kl = (*c_temp).kl;
         tmpval = (*c_temp).pval;

         gtmp2[ij] += ptmp2[kl]*tmpval;
         gtmp2[kl] += ptmp2[ij]*tmpval;
         if(optest[ij] && optest[kl]) {
            if (ij < iju) {
               gtmpo2[ij] += alpha1*ptmpo2[kl]*tmpval;
               gtmpo2[kl] += alpha1*ptmpo2[ij]*tmpval;
               }
            else if (kl < iju) {
               gtmpo2[ij] += ptmpo2[kl]*tmpval;
               gtmpo2[kl] += ptmpo2[ij]*tmpval;
               }
            else {
               gtmpo2[ij] += alpha3*ptmpo2[kl]*tmpval;
               gtmpo2[kl] += alpha3*ptmpo2[ij]*tmpval;
               }
            }
         }

      if (readflgc && j < num_bufs_c-1) {
         if(wherec==num_bufs_c) {
	    Pmat.bufpos = PSIO_ZERO;
            wherec=0;
            }
         wherec++;
         num=int_nums_c[wherec];
	 psio_read(Pmat.unit, Pmat.key, (char *) c_outbuf, sizeof(struct c_pkints)*num,
		   Pmat.bufpos, &(Pmat.bufpos));
         }
      }

   for(k=joff=0; k < num_ir ; k++) {
      if(nn=scf_info[k].num_so) {
         for(i=0; i < nn ; i++)
            for(j=0; j <= i ; j++) {
#if NEW2C
               scf_info[k].gmat[ioff[i]+j] += gtmp2[ioff[i+joff]+j+joff];
               scf_info[k].gmato[ioff[i]+j] += gtmpo2[ioff[i+joff]+j+joff];
#else
               scf_info[k].gmat[ioff[i]+j] = gtmp2[ioff[i+joff]+j+joff];
               scf_info[k].gmato[ioff[i]+j] = gtmpo2[ioff[i+joff]+j+joff];
#endif
               }
         joff += nn;
         if(print & 32) {
            fprintf(outfile,"\n gmat for irrep %s\n",scf_info[k].irrep_label);
            print_array(scf_info[k].gmat,nn,outfile);
            fprintf(outfile,"\n gmato for irrep %s\n",scf_info[k].irrep_label);
            print_array(scf_info[k].gmato,nn,outfile);
            }
         }
      }
   }

}} // namespace psi::cscf
