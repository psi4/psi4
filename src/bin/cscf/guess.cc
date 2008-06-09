/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/*--------------------------------------------------------

  guess.c: Function that reads the guessing 
  parameters from input and then either
  uses them to form an initial guess or
  calculates an initial guess from
  multiplicity and the charge using a 
  diagonalization of the core guess.

---------------------------------------------------------*/
#define EXTERN
#include "includes.h"
#include "common.h"
#include <libipv1/ip_lib.h>

namespace psi { namespace cscf {

void guess()
{
  struct symm *s;
  int i,j,k,m,l,o;
  int errcod;
  int nc,no,nh,nn;
  int netmp=0;
  int size;
  int optri;
  char *guess_opt,*occ_fix_str;
  reftype reftmp;
  int mguess;

  /*
   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":SCF");
  */

/*--- Set occupations to zero----*/
  for(i=0; i < num_ir; i++){
      scf_info[i].nclosed=0;
      scf_info[i].nopen=0;
    }

/* ---- Get reference type for comparison for restarting ----*/
  if(inflg == 1) reftmp = (reftype) chkpt_rd_ref();

/*---If first calculation either read occupations from DOCC or form an initial guess---*/
  /* check checkpoint flag to see if occupations in ckhpt file should be forced;
    this is needed to override DOCC vectors in finite difference calculations */ 
  if ( chkpt_rd_override_occ() ) {
	occ_read();
  }
  else {
      
	  errcod = ip_exist("DOCC",0);
	  if(errcod) {
	      fprintf(outfile,"\n  Using DOCC and SOCC to \n");
	      fprintf(outfile,"  determine occupations\n\n");
	      
	      errcod = ip_count("DOCC",&size,0);
	      if(errcod == IPE_OK && size != num_ir) {
		  fprintf(outfile,"\n  DOCC array is the wrong size\n");
		  fprintf(outfile,"  is %d, should be %d\n",size,num_ir);
		  exit(PSI_RETURN_FAILURE);
		}
	      if(errcod != IPE_OK && !iopen) {
		  fprintf(outfile,"\n  try adding some electrons buddy!\n");
		  fprintf(outfile,"  need DOCC\n");
		  exit(PSI_RETURN_FAILURE);
	      }
	      
	      if(iopen || uhf) {
		  errcod = ip_count("SOCC",&size,0);
		  if(errcod == IPE_OK && size != num_ir) {
		      fprintf(outfile,"\n  SOCC array is the wrong size\n");
		      fprintf(outfile,"  is %d, should be %d\n",size,num_ir);
		      exit(PSI_RETURN_FAILURE);
		  }
		  
		  errcod = ip_count("HOCC",&size,0);
		  if(errcod == IPE_OK && size != num_ir) {
		      fprintf(outfile,"\n  HOCC array is the wrong size\n");
		      fprintf(outfile,"  is %d, should be %d\n",size,num_ir);
		      exit(PSI_RETURN_FAILURE);
		  }
	      }
	      
	      for (i=0; i < num_ir; i++) {
		  errcod = ip_data("DOCC","%d"
				   ,&scf_info[i].nclosed,1,i);
		  if(iopen || uhf) errcod = 
				ip_data("SOCC","%d"
					,&scf_info[i].nopen,1,i);
		  if(iopen || uhf) errcod = 
				ip_data("HOCC","%d"
					,&scf_info[i].nhalf,1,i);
	      }
	     
	      /* STB (10/29/99) Make sure that DOCC and SOCC have the right
                 number of electrons in them */
	      for(i = 0;i<num_ir;i++)
		  netmp += (2*scf_info[i].nclosed) 
		      + scf_info[i].nopen + scf_info[i].nhalf;
	      if(netmp != nelec){
		 fprintf(outfile,"\n  DOCC and SOCC have the wrong number of");
		 fprintf(outfile,"\n  in them.  Should have %d total electrons"
			 ,nelec);
		 fprintf(outfile,"\n  and there are %d present\n\n\n",netmp);
		 exit(PSI_RETURN_FAILURE);
	      }
	  }
	  else if(inflg == 1 && (reftmp == refnum) )
	      occ_read();
	      
	  else{
	      fprintf(outfile, "  Using core guess to determine occupations\n\n"); 
	      
/* sort orbitals into two arrays, one with the symmetry label and 
   one with the eigenvalues */
	      
	      sortev();
	      
/* calculate the occupations */
	      
	      occ_calc();
	  }
  }
/* STB - 7/10/00 for DFT to ship the eigenvector */
/* calculate the number of total closed and open shells */
	  if(uhf){
	      for(i=0;i<num_ir;i++){
		  b_elec += scf_info[i].nclosed;
		  spin_info[1].scf_spin[i].nclosed = scf_info[i].nclosed;
		  a_elec += scf_info[i].nclosed+scf_info[i].nopen;
		  spin_info[0].scf_spin[i].nclosed = scf_info[i].nclosed
		      + scf_info[i].nopen;
	      }
	  }
	  else{
	      for(i=0; i < num_ir; i++){
		  n_closed += scf_info[i].nclosed;
	      }
	  }
/* output occupations to outfile */
	  occ_out();
	  
/* Setup occupation data for the calculation */

      n_open = 0;

/* Convert occupation data to symmetry info */
      
      for(i=0; i < num_ir; i++) {
	  s=&scf_info[i];
	  
	  if (s->nopen|| s->nhalf) n_open++;
	  if (nn = s->num_so) {
	      nc = s->nclosed;
	      no = s->nopen;
	      nh = s->nhalf;
	      
	      if (  (nn < nc + no + nh)
		    ||(nc < 0)
		    ||(no < 0)
		    ||(nh < 0)) {
		  const char* fmt
		      = "cscf: invalid number of electrons in irrep %d\n";
		  fprintf(stderr,fmt,i);
		  fprintf(outfile,fmt,i);
		  exit(PSI_RETURN_FAILURE);
		} 
	      
	      if(uhf){
		  spin_info[0].scf_spin[i].noccup = nc+no;
		  spin_info[1].scf_spin[i].noccup = nc;
		  if(nh != 0){
		      fprintf(outfile,"\nCannot use HOCC with UHF\n");
		      exit(PSI_RETURN_FAILURE);
		  }
	      }
	      
	      if(iopen) {
		  s->fock_eff = (double *) init_array(ioff[nn]);
		  s->fock_open = (double *) init_array(ioff[nn]);
		  s->pmato = (double *) init_array(ioff[nn]);
		  s->dpmato = (double *) init_array(ioff[nn]);
		  s->gmato = (double *) init_array(ioff[nn]);
	      }
	      if(twocon) {
		  s->pmat2 = (double *) init_array(ioff[nn]);
		  s->pmato2 = (double *) init_array(ioff[nn]);
	      }
	      if(uhf){
		  for(j=0; j < nc+no; j++)
		      spin_info[0].scf_spin[i].occ_num[j] = 1.0;
		  
		  for(j=0;j < nc; j++)
		      spin_info[1].scf_spin[i].occ_num[j] = 1.0;

	      }
	      else{
		  for (j=0; j < nc ; j++) {
		      s->occ_num[j] = 2.0;
		  }
		  for (j=nc; j < nc+no ; j++) {
		      s->occ_num[j] = 1.0;
		  }
		  if(nh) s->occ_num[nc+no] = nh*0.5;
	      }
	  }
      }
      
      optri = n_open*(n_open+1)/2;
      
      if (iopen) {
	  int mm1=1,mm2=1;
	  if(optri == 0){
	      fprintf(outfile,"\nNot an open shell molecule.\n");
	      fprintf(outfile,"Re-check opentype!!!!\n\n");
	      exit(PSI_RETURN_FAILURE);
	  }
	  
	  if(iter == 0 || print & 1){
	      fprintf(outfile,"\n  open-shell energy coeffs\n");
	      fprintf(outfile,"  open shell pair    alpha         beta\n");}
	  optri = (n_open*(n_open+1))/2;
	  
	  alpha = (double *) init_array(ioff[optri]);
	  beta = (double *) init_array(ioff[optri]);
	  
	  if (twocon) {
	      if(n_open == 2) {
		  alpha[0] = 0.0;
		  alpha[1] = 0.0;
		  alpha[2] = 0.0;
		  beta[0] = 0.0;
		  beta[1] = 1.0;
		  beta[2] = 0.0;
	      }
	      else {
		  fprintf(outfile,
			  "this program cannot handle same symmetry\n");
		  fprintf(outfile," tcscf. try SCFX\n");
		  exit(PSI_RETURN_FAILURE);
	      }
	  }
	  else if(singlet) {
	      if(n_open == 2) {
		  alpha[0] = 0.0;
		  alpha[0] = 0.0;
		  alpha[1] = 0.0;
		  alpha[2] = 0.0;
		  beta[0] = 1.0;
		  beta[1] = -3.0;
		  beta[2] = 1.0;
	      }
	      else {
		  fprintf(outfile,
			  "this program cannot handle same symmetry\n");
		  fprintf(outfile," singlets. try SCFX\n");
		  exit(PSI_RETURN_FAILURE);
	      }
	  }
	  else if(hsos) {
	      for(i=0; i < optri ; i++) {
		  alpha[i]=0.0;
		  beta[i]=1.0;
	      }
	  }
	  else {
	      for (i=0; i < optri ; i++) {
		  errcod = ip_data("ALPHA","%lf",&alpha[i],1,i);
		  errcod = ip_data("BETA","%lf",&beta[i],1,i);
		  beta[i] = -beta[i];
	      }
	  }
	  for (i=0; i < optri; i++) {
	      if(iter == 0 || print & 1)
		  fprintf(outfile,"        %d  %d       %f     %f\n",mm1,mm2,
			  alpha[i],-beta[i]);
	      mm2++;
	      if (mm2 > mm1) {
		  mm1++;
		  mm2 = 1;
	      }
	  }
      }
}       

}} // namespace psi::cscf
