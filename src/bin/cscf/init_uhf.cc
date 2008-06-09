/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/**************************************************

    init_uhf.c - separate code to initialize UHF

    By Shawn Brown

**************************************************/
static char *rcsid = "$Id: init_uhf.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <libchkpt/chkpt.h>

namespace psi { namespace cscf {

void init_uhf()
{
   int i,m,jj;
   int nn,isadr;
   int nkind,junk;
   int degen[20],*num_so;
   char char_dum[80];
   char **irr_labs;
   struct spin *sp;

   ioff[0] = 0;
   for (i = 1; i < 1024 ; i++) {
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
   
/* STB 6/30/99 Initialize structures to hold the information with spin
   spin in mind */
   
   spin_info = (struct spin *) malloc(sizeof(struct spin)*2);

   /* compute nsfmax */
   for(i=0; i < num_ir; i++) { nn = num_so[i]; if(nn > nsfmax) nsfmax = nn; }

   
   for(m=0;m<2;m++){
       spin_info[m].scf_spin = 
	   (struct symm *) malloc(sizeof(struct symm)*num_ir);
   }
   
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
	   
	   scf_info[i].smat = init_array(ioff[nn]);
	   scf_info[i].tmat = init_array(ioff[nn]);
	   scf_info[i].hmat = init_array(ioff[nn]);
	   scf_info[i].sahalf = block_matrix(nn,nn);
	   scf_info[i].pinv = block_matrix(nn,nn);
	   scf_info[i].occ_num = init_array(nn);
	   /* STB (6/30/99) - There is no P matrix in UHF only
	      J and K/2 */
	   
	   scf_info[i].dpmat = init_array(ioff[nn]);
	   scf_info[i].pmat = init_array(ioff[nn]);
	   scf_info[i].cmat = block_matrix(nn,nn);
	   /* STB(4/1/98) - Added array to store the eigenvalues of the
	      core hamiltonian for mo guessing*/
	   scf_info[i].hevals = init_array(nn);
	   for(m=0;m<2;m++){
	       sp = &spin_info[m];
           sp->scf_spin[i].irrep_label = scf_info[i].irrep_label;
	       sp->scf_spin[i].dpmat = init_array(ioff[nn]);
	       sp->scf_spin[i].dpmato = NULL;
	       sp->scf_spin[i].pmat = init_array(ioff[nn]);
	       sp->scf_spin[i].pmato = NULL;
	       sp->scf_spin[i].fock_pac = init_array(ioff[nn]);
	       sp->scf_spin[i].gmat = init_array(ioff[nn]);
	       sp->scf_spin[i].gmato = NULL;
	       sp->scf_spin[i].occ_num = init_array(nn);
	       sp->scf_spin[i].fock_evals = init_array(nn);
	       sp->scf_spin[i].cmat = block_matrix(nn,nn);
	       /* TDC(6/19/96) - Added array for saving original MO vector */
	       sp->scf_spin[i].cmat_orig = block_matrix(nn,nn);
	       /* STB(4/1/98) - Added array to store the eigenvalues of the
		  core hamiltonian for mo guessing*/
	       sp->scf_spin[i].hevals = init_array(nn);
	       /* Need separate XC Fock for KS DFT */
	       if (ksdft)
		 sp->scf_spin[i].xcmat = init_array(ioff[nn]);
	   }
       }
   }
   /* read in number of atoms and nuclear charges and total number of MO*/
   natom = chkpt_rd_natom();
   zvals = chkpt_rd_zvals();
   nbfso = chkpt_rd_nso();
   
   /* Character label for Spin */
   spin_info[0].spinlabel = "Alpha";
   spin_info[1].spinlabel = "Beta";
/* Initialize arrays to hold energy and symmetry arrays */
   ener_tot = init_array(nbfso);
   symm_tot = init_int_array(nbfso);
   
} 

}} // namespace psi::cscf
