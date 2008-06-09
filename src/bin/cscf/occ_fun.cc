/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/*****************************************************/
/*                                                   */
/*      occ_init - This function will take care of   */
/*                   all of the input parameters     */
/*                   concerning occupations          */
/*      occ_calc - Determines the occupations from    */
/*                multiplicity and charge            */
/*                                                   */
/*      By Shawn Brown - 6/23/99                     */
/*****************************************************/

#define EXTERN
#include "includes.h"
#include "common.h"
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>

namespace psi { namespace cscf {

/* Read the Opentype */

void occ_init(void){
    int i;
    int errcod; 
    int mguess=0;
    int open;

    iopen=0;
    mflag=0;
    uhf=special=twocon=hsos=singlet=0;
    
    num_ir = chkpt_rd_nirreps();
    
    multp=1;
    if(ip_exist("MULTIPLICITY", 0)) errcod = ip_data("MULTIPLICITY","%d",&multp,0);
    else if(ip_exist("MULTI", 0)) errcod = ip_data("MULTI","%d",&multp,0);
    else if(ip_exist("MULTP", 0)) errcod = ip_data("MULTP","%d",&multp,0);
    else {
	/*open = (int *) init_array(num_ir);*/
        if(ip_exist("SOCC",0)){
	    for(i=0;i<num_ir;i++){
		errcod = ip_data("SOCC","%d",&open,1,i);
		mguess += open;
	    }
	    
	    if(mguess == 1) multp = 2;
	    else if(mguess == 0) multp = 1;
	    else if(mguess == 2){
		fprintf(outfile,"\n  You must specify the MULTP keyword.");
		fprintf(outfile,"\n  I have no way of discerning between");
		fprintf(outfile,"\n  a triplet and an open shell singlet.\n");
		exit(PSI_RETURN_FAILURE);
	    }
	    else{
		fprintf(outfile,"\n The multiplicity will be ``highspin''.\n");
		multp = mguess + 1;
	    }
	}
	fprintf(outfile,"\n  I think the multiplicity is %d.\n",multp);
	fprintf(outfile,"  If this is wrong, please specify the MULTP keyword\n\n");
    }
	    
    
    if(!ip_exist("REFERENCE",0)) {
      reference = strdup("RHF");
      refnum = ref_rhf; /* RHF for file30 flag */
      ksdft = 0;  /* do Kohn-Sham DFT? Default - no */
    }
    else {
      errcod = ip_string("REFERENCE",&reference,0);
      if(!strcmp(reference,"ROHF")){
  	refnum = ref_rohf; /* ROHF for file30 flag */
  	iopen = 1;
  	if(multp == 1)
  	    singlet = 1;
  	else if(multp > 1)
  	    hsos = 1;
  	
      }
      else if(!strcmp(reference,"TWOCON")){
  	refnum = ref_tcscf; /*TCSCF for file30 flag */
  	iopen = 2;
  	twocon = 1;
      }
      else if(!strcmp(reference,"SPECIAL")){
  	/* NO FILE30 flag for Special */
  	iopen = 1;
  	special = 1;
      }
      else if(!strcmp(reference,"UHF")){
  	refnum = ref_uhf; /* UHF for file30 flag */
  	uhf = 1;
      }
      else if(!strcmp(reference,"RKS")){
  	refnum = ref_rks; /* flag for spin-restricted Kohn-Sham DFT */
  	ksdft = 1;
      }
      else if(!strcmp(reference,"UKS")){
  	refnum = ref_uks; /* flag for spin-unrestricted Kohn-Sham DFT */
  	uhf = 1;
  	ksdft = 1;
      }
      else{
  	if(multp != 1){
  	    fprintf(outfile,
  		    "\n Please specify an open shell reference\n");
  	    fprintf(outfile," with multpicity > 1\n");
  	    fprintf(outfile," multiplicity = %d\n",multp);
  	    fprintf(outfile," reference    = %s\n",reference);
  	    exit(PSI_RETURN_FAILURE);
  	}
      }
    }
	
/* Read Charge same as above */
    
    charge=0;
    errcod = ip_data("CHARGE","%d",&charge,0);

/* read in number of atoms and nuclear charges and total number of MO*/
    natom = chkpt_rd_natom();
    zvals = chkpt_rd_zvals();
    nbfso = chkpt_rd_nso();
    
/* Let's make sure that the molecule
   can have the specified multiplicity */
    
    nelec = 0;
    for(i=0; i < natom;i++)
	{
	    nelec = nelec + static_cast <int> (zvals[i]);
	}
    nelec = nelec - charge;
    
    if(multp == 1 && iopen == 0) {
	if(nelec%2==1) {
	    fprintf(outfile,"\n Impossible multiplicity with charge");
	    fprintf(outfile," and # of electrons specified\n\n");
	    fprintf(outfile,"\nMultiplicity = %d\nCharge = %d",multp,charge);
	    fprintf(outfile,"\nNumber of Electrons = %d\n",nelec);
	    exit(PSI_RETURN_FAILURE); }
    }
    else if(multp == 3 || (multp == 1 && singlet == 1)) {
	if(nelec%2==1) {
	    fprintf(outfile,"\n Impossible multiplicity with charge");
	    fprintf(outfile," and # of electrons specified\n\n");
	    fprintf(outfile,"\nMultiplicity = %d\nCharge = %d",multp,charge);
	    fprintf(outfile,"\nNumber of Electrons = %d\n",nelec);
	    exit(PSI_RETURN_FAILURE); }
    }
    else if(multp%2 == 1 && nelec%2 ==1){
	fprintf(outfile,"\nImpossible multiplicity with charge");
	fprintf(outfile," and # of electrons specified\n");
	fprintf(outfile,"\nMultiplicity = %d\nCharge = %d",multp,charge);
	fprintf(outfile,"\nNumber of Electrons = %d\n",nelec);
	exit(PSI_RETURN_FAILURE); 
    }
    else if(multp%2 == 0 && nelec%2 ==0){
	fprintf(outfile,"\nImpossible multiplicity with charge");
	fprintf(outfile," and # of electrons specified\n");
	fprintf(outfile,"\nMultiplicity = %d\nCharge = %d",multp,charge);
	fprintf(outfile,"\nNumber of Electrons = %d\n",nelec);
	exit(PSI_RETURN_FAILURE); 
    }
    else {
	fprintf(outfile,
		"\nCannot check consistency of the multiplicity\n");
	fprintf(outfile,"\nand number of electrons, double check\n");
	fprintf(outfile,"your occupations\n\n");
    }
}

void occ_calc(void){
    int i,j;
    int a,b;
    
    /* zero out occupations */
    for(i=0; i<num_ir; ++i) {
      scf_info[i].nclosed = 0;
      scf_info[i].nopen = 0;
      scf_info[i].nhalf = 0;
    } 
        
    if(multp == 1 && iopen == 0) {
	for(i = 0;i < (nelec/2);i++) 
		  scf_info[symm_tot[i]].nclosed++; 
    }
    
    else if(multp == 3 || (multp == 1 && singlet == 1)) {
	for(i = 0;i < ((nelec/2)-1);i++)
	    scf_info[symm_tot[i]].nclosed++;
	for(j = (nelec/2)-1;j < (nelec/2)+1;j++)
	    scf_info[symm_tot[j]].nopen++; 
    }
    
    else if(multp == 2) {
	for(i = 0;i < (nelec-1)/2;i++)
	    scf_info[symm_tot[i]].nclosed++;
	scf_info[symm_tot[((nelec-1)/2)]].nopen++;
    }
    
    else if(multp == 4) {
	if(nelec < 3){
	    fprintf(outfile,"\nNot enough electrons for a quartet\n\n");
	    exit(PSI_RETURN_FAILURE);
	}
	fprintf(outfile,
		"\n*****CSCF3.0 can only guess at highspin quartets*****\n");
	for(i=0;i < (nelec/2)-1;i++)
	    scf_info[symm_tot[i]].nclosed++;
	for(i=(nelec/2)-1;i < (nelec/2)+2;i++)
	    scf_info[symm_tot[i]].nopen++;
	    }
    /* STB (10/29/99) Can define all other highspin multiplicities by
       whether it is odd or even */
    
    else if(multp%2 == 0){
	if (nelec < multp - 1){
	    fprintf(outfile,
		    "\n  Not enough electrons for a Multiplicity %d \n\n"
		    ,multp);
	    exit(PSI_RETURN_FAILURE);
	}
	fprintf(outfile,
		"\n*****CSCF3.0 can only guess at highspin \n");
	fprintf(outfile,"     for Multplicity %d  *****\n\n",multp);
	
	a = (nelec/2)-((multp)-3);
	b = a+(multp-1);
		
	for(i=0;i<a;i++)
	    scf_info[symm_tot[i]].nclosed++;
	for(i=a;i<b;i++)
	    scf_info[symm_tot[i]].nopen++;
	    }
    
    else if(multp%2 == 1){
	if (nelec < multp - 1){
	    fprintf(outfile,
		    "\n  Not enough electrons for a Multiplicity %d \n\n"
		    ,multp);
	    exit(PSI_RETURN_FAILURE);
	}
	
	a = (nelec/2)-((multp)-3);
	b = a+(multp-1);
			
	for(i=0;i<a;i++)
	    scf_info[symm_tot[i]].nclosed++;
	for(i=a;i<b;i++)
	    scf_info[symm_tot[i]].nopen++;
    }

    /* determine alpha and beta occupations, if UHF */
    if(uhf){
      for(i=0;i<num_ir;i++){
        spin_info[1].scf_spin[i].nclosed = scf_info[i].nclosed;
        spin_info[0].scf_spin[i].nclosed = scf_info[i].nclosed
                                         + scf_info[i].nopen;
      }
    }

    // now assign occupation vectors
    for (i=0; i<num_ir; ++i) {
        const int nc = scf_info[i].nclosed;
        const int no = scf_info[i].nopen;
        const int nh = scf_info[i].nhalf;
        const int nmo = scf_info[i].num_mo;

        if (uhf) {
          int j;
          for (j=0; j < nc+no; j++)
            spin_info[0].scf_spin[i].occ_num[j] = 1.0;
          for ( ; j < nmo; j++)
            spin_info[0].scf_spin[i].occ_num[j] = 0.0;
          for (j=0; j < nc; j++)
            spin_info[1].scf_spin[i].occ_num[j] = 1.0;
          for ( ; j < nmo; j++)
            spin_info[1].scf_spin[i].occ_num[j] = 0.0;
        }
        else {
          int j;
          for (j=0; j < nc; j++) {
            scf_info[i].occ_num[j] = 2.0;
          }
          for ( ; j < nc+no; j++) {
            scf_info[i].occ_num[j] = 1.0;
          }
          for ( ; j < nc+no+nh; j++) {
            scf_info[i].occ_num[j] = 0.5;
          }
          for ( ; j < nmo; j++) {
            scf_info[i].occ_num[j] = 0.0;
          }
        }
      }

}

void occ_read(){
    int i;	
    
    fprintf(outfile,"\n  Reading Occupations from checkpoint file.\n");
    int* cldpi = chkpt_rd_clsdpi();
    int* openpi = chkpt_rd_openpi();
    
    for(i = 0;i < num_ir;i++){
	scf_info[i].nclosed=cldpi[i];
	scf_info[i].nopen=openpi[i];
    }
    free(cldpi);
    free(openpi);
}

void occ_out(void){
    
    int i;
    
    fprintf(outfile,"\n  Symmetry block:  ");
    for(i=0;i<num_ir;i++) {
	fprintf(outfile,"%4s  ",scf_info[i].irrep_label);
    }
    fprintf(outfile,"\n  DOCC:            ");
    for(i=0;i<num_ir;i++) {
	fprintf(outfile,"%3d   ",scf_info[i].nclosed);
    }
    fprintf(outfile,"\n  SOCC:            ");
    for(i=0;i<num_ir;i++) {
	fprintf(outfile,"%3d   ",scf_info[i].nopen);
    }
    fprintf(outfile,"\n\n");
    fflush(outfile);
}
    
}} // namespace psi::cscf
