/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "common.h"
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>

namespace psi { namespace cscf {

char *determine_functional(void);
char *determine_grid(void);

void scf_input(ip_value_t* ipvalue)
{
   int i,j,k,ijk,m;
   double elast;     
   double **scr_mat;
   char *alabel,*optyp,*wfn,*dertype,*guess,*jobtype;
   char *grid_str;
   char cjunk[80];
   int norder,*iorder,reordr;
   int nc,no,nh,nn,num_mo;
   int ncalcs;
   int optri,ierr,nat;
   int errcod;
   int size;
   int phase_chk;
   int mo_offset, so_offset;
   struct symm *s;
   reftype reftmp;
   int depth;
   int *mopi;

   /*
   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":SCF");
   */
  
   if(ipvalue) ip_print_value(stdout,ipvalue);

   errcod = ip_string("LABEL",&alabel,0);
   if(errcod == IPE_OK) fprintf(outfile,"  label        = %s\n",alabel);

   direct_scf = 0;
   errcod = ip_boolean("DIRECT",&direct_scf,0);
   /* Can do KS DFT direct only */
   if (ksdft) direct_scf=1;

   mixing = 0;
   errcod = ip_boolean("ORB_MIX",&mixing,0);

   /*-----------------------------------------------------
     Which type of guess to use. Sets inflg to:
     0 (AUTO,default) - check if there's an old vector
                        in file30, set inflg to 1 on yes,
		        to 2 otherwise.
     1                - use old vector in file30
     2 (GUESS=CORE)   - use core guess
    -----------------------------------------------------*/
   if(!ip_exist("GUESS",0))  {
     guess = strdup("AUTO");
     inflg=0;
   }
   else {
     errcod = ip_string("GUESS",&guess,0);
     if (!strcmp(guess,"AUTO"))
         inflg=0;
     else if (!strcmp(guess,"CORE"))
         inflg=2;
   }

   hcore_guess = 0;
   if(ip_exist("HCORE_GUESS",0)) {
     char* token;
     errcod = ip_string("HCORE_GUESS",&token,0);
     if (!strcmp(token,"NEW"))
       hcore_guess = 1;
   }
   
   reset_occ = 0;
   errcod = ip_boolean("RESET_OCCUPATIONS",&reset_occ,0);

   reordr = 0;
   norder = 0;
   errcod = ip_boolean("REORDER",&reordr,0);
   if(reordr) {
     norder = 1;

     errcod = ip_count("MOORDER",&size,0);
     errchk(errcod,"MOORDER");
     if(errcod != IPE_OK) {
       fprintf(outfile,"\ncannot find MOORDER. calculation continuing\n");
       norder=0;
       reordr=0;
       }
     else {
       if(size != nbasis) {
         fprintf(outfile,"\n you have not given enough mos to MOORDER\n");
         exit(PSI_RETURN_FAILURE);
         }
       iorder = (int *) malloc(sizeof(int)*size);
       for(i=0; i < size ; i++) {
         errcod = ip_data("MOORDER","%d",&iorder[i],1,i);
         errchk(errcod,"MOORDER");
         }
       }
     }
   
   /* Remove after debugging.  Stop cscf right before going to cints */
   exitflag = 0;
   errcod = ip_boolean("EXIT_CINTS",&exitflag,0);
   
   itmax = 100;
   errcod = ip_data("MAXITER","%d",&itmax,0);

   it_diis = 0;
   errcod = ip_data("DIISSTART","%d",&it_diis,0);
   
   print = 0;
   errcod = ip_data("IPRINT","%d",&print,0);

   fock_typ = 0;
   errcod = ip_data("FOCK_TYPE","%d",&fock_typ,0);

   second_root = 0;
   if (twocon) {
       errcod = ip_boolean("SECOND_ROOT",&second_root,0);
     }
   
   icheck_rot = 1;
   errcod = ip_boolean("CHECK_ROT",&icheck_rot,0);
   
   check_mo_orthonormality = 0;
   errcod = ip_boolean("CHECK_MO_ORTHONORMALITY",&check_mo_orthonormality,0);

   ndiis = (iopen) ? 4 : 6;
   if(twocon) ndiis = 3;
   errcod = ip_data("NDIIS","%d",&ndiis,0);

   if(ipvalue) ip_print_tree(stdout,NULL);

   scf_conv = 7;
   if(ipvalue) ip_print_value(stdout,ipvalue);
   errcod = ip_string("WFN",&wfn,0);
   if(ipvalue) ip_print_value(stdout,ipvalue);
   errcod = ip_string("DERTYPE",&dertype,0);
   if(errcod == IPE_KEY_NOT_FOUND) {
     dertype = (char *) malloc(sizeof(char)*5);
     strcpy(dertype,"NONE");
     }
   if(strcmp(wfn,"SCF")) scf_conv = 10;
   if(!strcmp(dertype,"FIRST")) scf_conv = 10;
   if(!strcmp(dertype,"SECOND")) scf_conv = 12;

   /* for freq finite difference calculations */
   errcod = ip_string("JOBTYPE",&jobtype,0);
   if (errcod == IPE_OK) 
     if (!strcmp(jobtype,"FREQ")) scf_conv = 12;

   errcod = ip_data("CONVERGENCE","%d",&scf_conv,0);

   if (ksdft){
       functional = (char *)determine_functional();
       grid_str = (char *)determine_grid();
   }
   
   if(ipvalue) ip_print_value(stdout,ipvalue);
   fprintf(outfile,"  wfn          = %s\n",wfn);
   fprintf(outfile,"  reference    = %s\n",reference);
   if (ksdft) {
   fprintf(outfile,"  functional   = %s\n",functional);
   fprintf(outfile,"  DFT grid     = %s\n",grid_str);
   }
   fprintf(outfile,"  multiplicity = %d\n",multp);
   fprintf(outfile,"  charge       = %d\n",charge);
   fprintf(outfile,"  direct       = %s\n",(direct_scf) ? "true" : "false");
   if(direct_scf)
   fprintf(outfile,"  dyn_acc      = %s\n",(dyn_acc) ? "true" : "false");
   fprintf(outfile,"  dertype      = %s\n",dertype);
   fprintf(outfile,"  convergence  = %d\n",scf_conv);
   fprintf(outfile,"  maxiter      = %d\n",itmax);
   fprintf(outfile,"  guess        = %s\n",guess); free(guess);
   if(print) fprintf(outfile,"  iprint       = %d\n",print);
   if (second_root)
     fprintf(outfile,"  second_root = TRUE\n");

   diisflg = 1;
   errcod = ip_boolean("DIIS",&diisflg,0);
   /* the convention is to set it to the opposite of it's meaning? */
   diisflg = !diisflg;

   fprintf (outfile,"\n  nuclear repulsion energy %22.13f\n",repnuc);
   fflush(outfile);

   nat    = chkpt_rd_natom();
   ncalcs = chkpt_rd_ncalcs();

/* if inflg is 0 and this isn't the first calc, then get the old vector */
/* from file30.  if inflg is 2, just use core hamiltonian guess */
/* if inflg is 1, get old vector no matter what ncalcs is */

   if ((inflg==0 && ncalcs) || inflg == 1) {
       
       inflg = 1;
       optri = abs(chkpt_rd_iopen());
       reftmp = (reftype) chkpt_rd_ref();
       
       fprintf(outfile,"\n  using old vector from file30 as initial guess\n");
       
/* get old energy from file30 */
       
       elast = chkpt_rd_escf();
       fprintf(outfile,"  energy from old vector: %14.8f\n",elast);
       
       so_offset = 0;
       mo_offset = 0;

       /* Add MO's per/irrep for scf_info */
       mopi = chkpt_rd_orbspi();
       for(k=0; k < num_ir; k++) scf_info[k].num_mo = mopi[k];
       free(mopi);

/* ----------------------------------------------------
** This is the UHF part of the restarting algorithm
** STB (10/29/99)
**
**----------------------------------------------------*/
       
       if(uhf){
	   
	   /* if the reference is not UHF, then just read in the vector for the 
	      restricted calculation into both */
	   if(reftmp != ref_uhf && reftmp != ref_uks){
	       for(k=0; k < num_ir ; k++) {
		   s = &scf_info[k];
		   if(nn=s->num_so) {
		       spin_info[0].scf_spin[k].cmat = chkpt_rd_scf_irrep(k);
		       spin_info[1].scf_spin[k].cmat = chkpt_rd_scf_irrep(k);
		   }
	       }
	   }
	   else{
	       for(k=0; k < num_ir ; k++) {
		   s = &scf_info[k];
		   if(nn=s->num_so) {
		       spin_info[0].scf_spin[k].cmat = chkpt_rd_alpha_scf_irrep(k);
		       spin_info[1].scf_spin[k].cmat = chkpt_rd_beta_scf_irrep(k);
		   }
	       }
	   }
	   
	   for(m=0;m<2;m++){
	       for(k=0; k < num_ir; k++) {
		   s = &scf_info[k];
		   if(nn=s->num_so) {
		     num_mo = s->num_mo;
		       for(i=0; i < nn; i++) 
		           for(j=0; j < num_mo; j++)
			       spin_info[m].scf_spin[k].cmat_orig[i][j] 
				   = spin_info[m].scf_spin[k].cmat[i][j];
		   }
	       }
	   }
	   phase_check = 1;

/* reorder vector if norder = 1 */
	   
	   if (norder) {
	       int loff = 0;
	       int jnew;
	       
	       /* TDC(6/19/96) - If the vector is re-ordered, don't allow
		  phase_checking */
	       phase_check = 0;
	       
	       fprintf(outfile,"\n  mo's will be reordered\n");
	       for (m=0;m<2;m++){
		   for (i=0; i < num_ir; i++) {
		       s = &scf_info[i];
		       if (nn=s->num_so) {
			   num_mo = s->num_mo;
			   scr_mat = (double **) init_matrix(nn,num_mo);
			   for (j=0; j < num_mo; j++) {
			       jnew = iorder[j+loff]-1;
			       for (k=0; k < nn ; k++) {
				   scr_mat[k][j]
				       =spin_info[m].scf_spin[i].cmat[k][jnew];
			       }
			   }
			   for (j=0; j < nn ; j++)
			       for (k=0; k < num_mo ; k++) 
				   spin_info[m].scf_spin[i].cmat[j][k] = scr_mat[j][k];
			   
			   fprintf(outfile,"\n reordered %s mo's for irrep %s\n",
				   spin_info[m].spinlabel,s->irrep_label);
			   print_mat(spin_info[m].scf_spin[i].cmat,nn,num_mo,outfile);
			   loff += num_mo;
			   free_matrix(scr_mat,nn);
		       }
		   }
	       }
	       free(iorder);
	   }
       }
       else{
	   for(k=0; k < num_ir ; k++) {
	       s = &scf_info[k];
	       if(nn=s->num_so) {
		   s->cmat = chkpt_rd_scf_irrep(k);
	       }
	   }

/* TDC(6/19/96) - Make a copy of the vector for later MO phase
   checking and temporarily set the phase_check flag to true */
	   
	   for(k=0; k < num_ir; k++) {
	       s = &scf_info[k];
	       if(nn=s->num_so) {
		   num_mo = s->num_mo;
		   for(i=0; i < nn; i++) 
		       for(j=0; j < num_mo; j++)
                           s->cmat_orig[i][j] = s->cmat[i][j];
	       }
	   }

	   phase_check = 1;
	   
/* reorder vector if norder = 1 */
	   
	   if (norder) {
	       int loff = 0;
	       int jnew;
	       
	       /* TDC(6/19/96) - If the vector is re-ordered, don't allow
		  phase_checking */
	       phase_check = 0;
	       
	       fprintf(outfile,"\n  mo's will be reordered\n");
	       for (i=0; i < num_ir; i++) {
		   s = &scf_info[i];
		   if (nn=s->num_so) {
		       num_mo = s->num_mo;
		       scr_mat = (double **) init_matrix(nn,num_mo);
		       for (j=0; j < num_mo; j++) {
			   jnew = iorder[j+loff]-1;
			   for (k=0; k < nn ; k++) {
			       scr_mat[k][j]=s->cmat[k][jnew];
			   }
		       }
		       for (j=0; j < nn ; j++)
			   for (k=0; k < num_mo ; k++) s->cmat[j][k] = scr_mat[j][k];
		       
		       fprintf(outfile,"\n reordered mo's for irrep %s\n",
			       s->irrep_label);
		       print_mat(s->cmat,nn,num_mo,outfile);
		       loff += num_mo;
		       free_matrix(scr_mat,nn);
		   }
	       }
	       free(iorder);
	   }
       }
   }
   else {
       inflg = 2;
       fprintf(outfile,"  first run, so defaulting to core-hamiltonian guess\n");
       /* TDC(6/19/96) - If not starting from old vector, don't allow
	  phase checking */
       for(k=0; k < num_ir; k++) scf_info[k].num_mo = 0;
       phase_check = 0;
   }

/* TDC(6/20/96) - Check to see if the user will let us do phase
   correction.  The default has already been set above. */
   phase_chk = 1;
   errcod = ip_boolean("PHASE",&phase_chk,0);
   if(phase_check && phase_chk) phase_check = 1;
   else phase_check = 0;

/* read in damping factor and level shift */

   dampsv= (iopen) ? 0.02 : 0.0;
   if(twocon) dampsv = 0.01;
   errcod = ip_data("DIISDAMP","%lf",&dampsv,0);

   lshift=1.0;
   errcod = ip_data("LEVELSHIFT","%lf",&lshift,0);
   if(!iopen && fabs(lshift) > 0.0) lshift = 0.1;
   stop_lshift=10;
   errcod = ip_data("STOP_LEVELSHIFT","%d",&stop_lshift,0);

   dampd=1.0;
   errcod = ip_data("DAMPD","%lf",&dampd,0);

   dampo=1.0;
   errcod = ip_data("DAMPO","%lf",&dampo,0);

   fprintf(outfile,"\n  level shift                      = %f\n",lshift);
   fprintf(outfile,"\n  level shifting will stop after %d cycles\n",stop_lshift);
   if(!diisflg) {
      fprintf(outfile,"  diis scale factor                = %f\n",dampsv+1.0);
      fprintf(outfile,"  iterations before extrapolation  = %d\n",it_diis);
      fprintf(outfile,"  %d error matrices will be kept\n",ndiis);
      }
   else fprintf(outfile,"\n  diis turned off\n");

   switch (fock_typ) {
      case 0:
         break;
      case 1:
         fprintf(outfile,"\n  a fock matrix for high spin will be used\n");
         fprintf(outfile,"  this form may not work well with diis\n");
         break;
      default:
         fprintf(outfile,"\n  an experimental fock matrix will be used\n");
         fprintf(outfile,"  the management will not be held responsible for the results\n");
      }


   /* EFV 10/24/98 Check if delete integrals */
   delete_ints = 0;
   if(!strcmp(wfn,"SCF") && (!strcmp(dertype,"FIRST") || !strcmp(dertype,"NONE")))
     delete_ints = 1;
   errcod = ip_boolean("DELETE_INTS",&delete_ints,0);
     /* These keywords will work only with IWL format */
   delete_1e = delete_ints;
   errcod = ip_boolean("DELETE_1E",&delete_1e,0);
   delete_2e = delete_ints;
   errcod = ip_boolean("DELETE_2E",&delete_2e,0);

   free(dertype);
   
   fflush(outfile);
}

}} // namespace psi::cscf
